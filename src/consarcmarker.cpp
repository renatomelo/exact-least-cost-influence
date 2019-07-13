#include "consarcmarker.h"

struct SCIP_ConsData
{
    Digraph *graph;
    Arc arc;
    int nPropagatedVars;         /* number of variables that existed, the last time,
                                  *  the related node was propagated, used to determine
                                  *  whether the constraint should be repropagated */
    CONSTYPE type;               // stores whether the arc has to be in the solution or not
    int nPropagations;           // stores the number propagations runs of this constraint
    unsigned int propagated : 1; // is constraint already propagated?
    SCIP_NODE *node;             // the node in the B&B-tree at which the cons is sticking
};

/**
 * fixes variables
 */
SCIP_RETCODE checkVariable(
    SCIP *scip,
    SCIP_CONSDATA *consdata, // constraint data
    SCIP_VAR *var,           // variables to check
    int *nFixedVars,         // pointer to store the number of fixed variables
    SCIP_Bool *cuttoff       // pointer to store if a cutoff was detected
)
{
    SCIP_VARDATA *vardata;
    CONSTYPE type;
    SCIP_Bool fixed;
    SCIP_Bool infeasible;

    assert(consdata != NULL);
    assert(var != NULL);
    assert(nFixedVars != NULL);
    assert(cuttoff != NULL);

    // if variables is locally fixed to zero continue
    if (SCIPvarGetUbLocal(var) < 0.5)
        return SCIP_OKAY;

    //check if the arc which corresponds to the varialbe is feasible for this constraint
    vardata = SCIPvarGetData(var);

    if (type == WITHOUT)
    {
        SCIP_CALL(SCIPfixVar(scip, var, 0.0, &infeasible, &fixed));
        if (infeasible)
        {
            assert(SCIPvarGetLbLocal(var) > 0.5);
            SCIPinfoMessage(scip, NULL, "-> cutoff\n");
            (*cuttoff) = TRUE;
        }
        else
        {
            assert(fixed);
            (*nFixedVars)++;
        }
    }
    // otherwise the variable is in the solution
    else
    {
        SCIP_CALL(SCIPfixVar(scip, var, 1.0, &infeasible, &fixed));
        if (infeasible)
        {
            assert(SCIPvarGetUbLocal(var) < 0.5);
            SCIPinfoMessage(scip, NULL, "-> cutoff\n");
            (*cuttoff) = TRUE;
        }
        else
        {
            assert(fixed);
            (*nFixedVars)++;
        }
    }

    return SCIP_OKAY;
}

/**
 * fixes variables to zero if the corresponding arcs are not valid for this node 
 * due to branching
 */
SCIP_RETCODE fixVariables(
    SCIP *scip,
    SCIP_CONSDATA *consdata,
    ArcSCIPVarMap &z,
    SCIP_RESULT *result)
{
    int nFixedVars = 0;
    SCIP_Bool cutoff = FALSE;

    SCIPinfoMessage(scip, NULL, "check variables %d to %d\n", consdata->nPropagatedVars, nVars);

    for (int i = consdata->nPropagatedVars; i < nVars; i++)
    {
        SCIP_CALL( checkVariable(scip, consdata, vars[i], &nFixedVars, &cutoff) );
    }
    SCIPinfoMessage(scip, NULL, "fixed %d variables locally\n", nFixedVars);

    if (cutoff)  
        *result = SCIP_CUTOFF;
    else if (nFixedVars > 0)  
        *result = SCIP_REDUCEDDOM;

    return SCIP_OKAY;
}

/**
 * frees constraint data
 */
SCIP_DECL_CONSDELETE(ConshdlrArcMarker::scip_delete)
{
    assert(consdata != NULL);

    SCIPfreeBlockMemory(scip, consdata);

    return SCIP_OKAY;
}

/**
 * transforms constraint data into data belonging to the trasformed problem
 */
SCIP_DECL_CONSTRANS(ConshdlrArcMarker::scip_trans)
{
    SCIP_CONSDATA *sourcedata;
    SCIP_CONSDATA *targetdata;

    sourcedata = SCIPconsGetData(sourcecons);
    assert(sourcedata != NULL);

    //create constraint data for the target constraint
    SCIP_CALL(
        createConsData(scip, &targetdata, sourcedata->arc, sourcedata->type, sourcedata->node)
    );

    //create target constraint
    SCIP_CALL(
        SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr,
                         targetdata, SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
                         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
                         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
                         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons))
    );
    return SCIP_OKAY;
}

/**
 * domain propagation method of constraint handler
 */
SCIP_DECL_CONSPROP(ConshdlrArcMarker::scip_prop)
{
    SCIP_CONSDATA *consdata;
    SCIP_VAR **vars;
    int nVars;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(result != NULL);

    SCIPinfoMessage(scip, NULL, "propagation constraints of constraint handler <" CONSHDLR_NAME ">\n");

    // get vars and nvars if needed

    *result = SCIP_DIDNOTFIND;

    for (int i = 0; i < nconss; i++)
    {
        consdata = SCIPconsGetData(conss[i]);

        // chech if all previously generated variables are valid for this constraint
        //assert(checkConsData(instance, consdata, TRUE));

        // remember of verify about the NDEBUG flag of SCIP
        if (!consdata->propagated)
        {
            SCIPinfoMessage(scip, NULL, "propagate constraint <%s> ", SCIPconsGetName(conss[i]));
            printConsData(instance, consdata);

            SCIP_CALL(fixVariables(scip, consdata, z, result));
            consdata->nPropagations++;

            if (*result != SCIP_CUTOFF)
            {
                consdata->propagated = TRUE;
                consdata->nPropagatedVars = nVars;
            }
            else
                break;
        }

        //check if constraint is completely propagated
        //assert(checkConsData(instance, consdata, FALSE));
    }

    return SCIP_OKAY;
}

/**
 * constraint activation notification method of constraint handler
 */
SCIP_DECL_CONSACTIVE(ConshdlrArcMarker::scip_active)
{
    SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->nPropagatedVars <= instance.m);

    SCIPinfoMessage(scip, NULL,
                    "activate constraint <%s> at node <%" SCIP_LONGINT_FORMAT "> in depth <%d>: ",
                    SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node),
                    SCIPnodeGetDepth(consdata->node));
    //SCIPinfoMessage(consdataPrint(scip, consdata, NULL));

    if (consdata->nPropagatedVars != instance.m)
    {
        SCIPinfoMessage(scip, NULL, "-> mark constraint to be repropagated\n");
        consdata->propagated = FALSE;
        SCIP_CALL(SCIPrepropagateNode(scip, consdata->node));
    }

    return SCIP_OKAY;
}

SCIP_DECL_CONSDEACTIVE(ConshdlrArcMarker::scip_deactive)
{
    SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

    SCIPinfoMessage(scip, NULL,
                    "deactivate constraint <%s> at node <%" SCIP_LONGINT_FORMAT "> in depth <%d>: ",
                    SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node),
                    SCIPnodeGetDepth(consdata->node));

    // set the number of propagated varables to current number of ARC variables in SCIP
    consdata->nPropagatedVars = instance.m;

    return SCIP_OKAY;
}

/**
 * variable deletion method of constraint handler
 */
SCIP_DECL_CONSDELVARS(ConshdlrArcMarker::scip_delvars)
{
    return SCIP_OKAY;
}

SCIP_DECL_CONSPRINT(ConshdlrArcMarker::scip_print)
{
    SCIP_CONSDATA *consdata;
    Arc a;

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    a = consdata->arc;
    DNode u = consdata->graph->source(a);
    DNode v = consdata->graph->target(a);

    SCIPinfoMessage(scip, NULL, "branching on the arc (%s,%s)", instance.nodeName[u], instance.nodeName[v]);
}

SCIP_RETCODE createConsData(
    SCIP *scip,
    SCIP_CONSDATA **consdata, //pointer to store the constraint data
    Arc arc,
    CONSTYPE type,   //stores whether the arc have to be in the solution or not
    SCIP_NODE *node) //the node in the B&B-tree at which the cons is sticking
{
    assert(scip != NULL);
    assert(consdata != NULL);
    assert(arc != NULL);
    assert(type == WITH || type == WITHOUT);

    SCIP_CALL(SCIPallocBlockMemory(scip, consdata));

    (*consdata)->arc = arc;
    (*consdata)->type = type;
    (*consdata)->nPropagatedVars = 0;
    (*consdata)->nPropagations = 0;
    (*consdata)->propagated = FALSE;
    (*consdata)->node = node;

    return SCIP_OKAY;
}

void printConsData(GLCIPInstance &instance, SCIP_CONSDATA *consdata)
{
    DNode u = instance.g.source(consdata->arc);
    DNode v = instance.g.target(consdata->arc);
    std::cout << "arc (" << instance.nodeName[u] << "," << instance.nodeName[v]
              << ") at node " << SCIPnodeGetNumber(consdata->node) << "is "
              << consdata->type ==
            WITH
        ? ""
        : "NOT"
              << " on the solution" << std::endl;
}

/**
 * creates and captures the arc marker constraint
 */
SCIP_RETCODE ConshdlrArcMarker::createConsArcMarker(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name,
    Arc arc,
    CONSTYPE type,
    SCIP_NODE *node,
    SCIP_Bool local)
{
    SCIP_CONSHDLR *conshdlr;
    SCIP_CONSDATA *consdata;

    // find the arc marker constraint handler
    conshdlr = SCIPfindConshdlr(scip, "arcmarker");
    if (conshdlr == NULL)
    {
        SCIPerrorMessage("arc marker constraint handler not found\n");
        return SCIP_PLUGINNOTFOUND;
    }

    // create constraint data
    SCIP_CALL(createConsData(scip, &consdata, arc, type, node));

    //create constraint
    SCIP_CALL(
        SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE,
                       TRUE, local, FALSE, FALSE, FALSE, TRUE));

    SCIPinfoMessage(scip, NULL, "created constraint: ");
    printConsData(instance, consdata);

    return SCIP_OKAY;
}