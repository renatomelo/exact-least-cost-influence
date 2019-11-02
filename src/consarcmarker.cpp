#include "consarcmarker.h"
#include "pricer_glcip.h"

using namespace arcmarker;

struct SCIP_ConsData
{
    Digraph *graph;
    ArcIt arc;
    SCIP_VAR *arcVar;
    int nPropagatedVars;         /* number of variables that existed, the last time,
                                  *  the related node was propagated, used to determine
                                  *  whether the constraint should be repropagated */
    CONSTYPE type;               // stores whether the arc has to be in the solution or not
    int nPropagations;           // stores the number propagations runs of this constraint
    unsigned int propagated : 1; // is constraint already propagated?
    SCIP_NODE *node;             // the node in the B&B-tree at which the cons is sticking
};

SCIP_RETCODE createConsData(
    SCIP *scip,
    SCIP_CONSDATA **consdata, //pointer to store the constraint data
    ArcIt& arc,
    SCIP_VAR *arcVar,
    CONSTYPE type,   //stores whether the arc have to be in the solution or not
    SCIP_NODE *node) //the node in the B&B-tree at which the cons is sticking
{
    assert(scip != NULL);
    assert(consdata != NULL);
    //assert(arc != NULL);
    assert(type == WITH || type == WITHOUT);

    SCIP_CALL(SCIPallocBlockMemory(scip, consdata));

    (*consdata)->arc = arc;
    (*consdata)->arcVar = arcVar;
    (*consdata)->type = type;
    (*consdata)->nPropagatedVars = 0;
    (*consdata)->nPropagations = 0;
    (*consdata)->propagated = FALSE;
    (*consdata)->node = node;

    return SCIP_OKAY;
}

void printConsData(SCIP_CONSDATA *consdata)
{
    /* DNode u = instance.g.source(consdata->arc);
    DNode v = instance.g.target(consdata->arc);
    std::cout << "arc (" << instance.nodeName[u] << "," << instance.nodeName[v]
              << ") at node " << SCIPnodeGetNumber(consdata->node) << "is "
              << (consdata->type == WITH ? " " : "NOT") << " on the solution" << std::endl; */

    std::cout << "arc variable  (" << SCIPvarGetName(consdata->arcVar) << ") at node "
              << SCIPnodeGetNumber(consdata->node) << " is "
              << (consdata->type == WITH ? "" : "NOT") << " on the solution" << std::endl;
}

/**
 * fixes variables
 * TODO: probably in this function I need to create a copy of the graph and remove edges
 */
SCIP_RETCODE checkVariable(
    SCIP *scip,
    GLCIPInstance& instance,
    DNodeSCIPVarMap& x,
    SCIP_CONSDATA *consdata, // constraint data
    SCIP_VAR *var,           // variables to check
    int *nFixedVars,         // pointer to store the number of fixed variables
    SCIP_Bool *cuttoff       // pointer to store if a cutoff was detected
)
{
    //std::cout << "checkVariable() function\n";
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
    if (consdata->type == WITHOUT)
    {
        SCIP_CALL(SCIPfixVar(scip, var, 0.0, &infeasible, &fixed));
        //std::cout << "both bounds of variable " << SCIPvarGetName(var) << " are set to 0\n";

        //SCIP_PRICER* pricer = SCIPfindPricer(scip, "GLCIP_pricer");
        ObjPricerGLCIP* pricer_glcip = (ObjPricerGLCIP*)  SCIPfindObjPricer(scip, "GLCIP_pricer");
        //pricer_glcip->isOnSolution[consdata->arc] = -1;
        //SCIP_PRICERDATA* pricerdata =  SCIPpricerGetData(scip, pricer);
        //SCIPpricerSetData()

        if (infeasible)
        {
            std::cout << "IS INFEASIBLE\n";
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
        // fix the variable x_i to 1 meaning that the source(arc) is active
        //SCIP_CALL(SCIPfixVar(scip, x[instance.g.source(consdata->arc)], 1.0, &infeasible, &fixed));
        
        //std::cout << "both bounds of variable " << SCIPvarGetName(var) << " are set to 1\n";

        ObjPricerGLCIP* pricer_glcip = (ObjPricerGLCIP*)  SCIPfindObjPricer(scip, "GLCIP_pricer");
        //pricer_glcip->isOnSolution[consdata->arc] = 1;

        if (infeasible)
        {
            std::cout << "IS INFEASIBLE\n";
            assert(SCIPvarGetLbLocal(var) < 0.5);
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
    GLCIPInstance& instance,
    DNodeSCIPVarMap& x,
    SCIP_CONSDATA *consdata,
    SCIP_VAR *arc_var,
    SCIP_RESULT *result)
{
    int nFixedVars = 0;
    SCIP_Bool cutoff = FALSE;

    SCIP_CALL(checkVariable(scip, instance, x, consdata, arc_var, &nFixedVars, &cutoff));

    //SCIPinfoMessage(scip, NULL, "variable locally fixed\n");

    if (cutoff)
        *result = SCIP_CUTOFF;
    else if (nFixedVars > 0)
        *result = SCIP_REDUCEDDOM;

    return SCIP_OKAY;
}

/**
 * Frees all data of a constraint
 */
SCIP_DECL_CONSDELETE(ConshdlrArcMarker::scip_delete)
{
    assert(consdata != NULL);

    SCIPfreeBlockMemory(scip, consdata);

    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOLP(ConshdlrArcMarker::scip_enfolp)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOPS(ConshdlrArcMarker::scip_enfops)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_DECL_CONSLOCK(ConshdlrArcMarker::scip_lock)
{
    return SCIP_OKAY;
}

SCIP_DECL_CONSCHECK(ConshdlrArcMarker::scip_check)
{
    *result = SCIP_FEASIBLE;
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
        createConsData(scip, &targetdata, sourcedata->arc, sourcedata->arcVar, sourcedata->type, sourcedata->node));

    //create target constraint
    SCIP_CALL(
        SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr,
                       targetdata, SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
                       SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
                       SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
                       SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
                       SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)));
    return SCIP_OKAY;
}

/**
 * Propagation callback of the constraint handler, needed to fix variables
 * for enforcing the branching decision.
 * TODO read about the SCIPrepropagateNode() method
 */
SCIP_DECL_CONSPROP(ConshdlrArcMarker::scip_prop)
{
    SCIP_CONSDATA *consdata;
    //int nVars = instance.m;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(result != NULL);

    //std::cout << "propagation of constraint handler <" << CONSHDLR_NAME << ">\n";

    *result = SCIP_DIDNOTFIND;

    for (int i = 0; i < nconss; i++)
    {
        consdata = SCIPconsGetData(conss[i]);
        /* std::cout << "number of propagated variables: " << consdata->nPropagatedVars
                  << std::endl; */

        // chech if all previously generated variables are valid for this constraint
        //assert(checkConsData(instance, consdata, TRUE));

        // remember of verify about the NDEBUG flag of SCIP
        if (!consdata->propagated)
        {
            /* std::cout << "propagate constraint <" << SCIPconsGetName(conss[i]) << ">:";
            printConsData(consdata); */

            SCIP_CALL(fixVariables(scip, instance, x, consdata, consdata->arcVar, result));
            consdata->nPropagations++;

            if (*result != SCIP_CUTOFF)
                consdata->propagated = TRUE;
            else
                break;
        }

        //check if constraint is completely propagated
        //assert(checkConsData(instance, consdata, FALSE));
    }

    return SCIP_OKAY;
}

/**
 * Constraint activation notification method of constraint handler.
 * This method is always called when a node is entered on which the constraint
 * has been added. Here, we apply the changes to the pricing data structures. 
 */
SCIP_DECL_CONSACTIVE(ConshdlrArcMarker::scip_active)
{
   /*  SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL); */
    //assert(consdata->nPropagatedVars <= instance.m);

    /* std::cout << "activate constraint <" << SCIPconsGetName(cons) << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata); */

    //std::cout << "number of propagated variables: " << consdata->nPropagatedVars << std::endl;
    /*  if (consdata->nPropagatedVars != instance.m)
    { */
   /*  SCIPinfoMessage(scip, NULL, "-> mark constraint to be repropagated\n");
    consdata->propagated = FALSE;
    SCIP_CALL(SCIPrepropagateNode(scip, consdata->node)); */
    //}

    return SCIP_OKAY;
}

/**
 * The CONSDEACTIVE method will be called if the node is left again. Since the 
 * CONSACTIVE and CONSDEACTIVE methods of different constraints are always called in a
 * stack-like fashion, this should be exactly what you need.
 */
SCIP_DECL_CONSDEACTIVE(ConshdlrArcMarker::scip_deactive)
{
    SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

    /* std::cout << "deactivate constraint <" << SCIPconsGetName(cons) << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata); */

    ObjPricerGLCIP* pricer_glcip = (ObjPricerGLCIP*)  SCIPfindObjPricer(scip, "GLCIP_pricer");
    //pricer_glcip->isOnSolution[consdata->arc] = 0;

    return SCIP_OKAY;
}

/**
 * variable deletion method of constraint handler
 */
SCIP_DECL_CONSDELVARS(ConshdlrArcMarker::scip_delvars)
{
    return SCIP_OKAY;
}

/* SCIP_DECL_CONSPRINT(ConshdlrArcMarker::scip_print)
{
    SCIP_CONSDATA *consdata;
    Arc a;

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    a = consdata->arc;
    DNode u = consdata->graph->source(a);
    DNode v = consdata->graph->target(a);

    SCIPinfoMessage(scip, NULL, "branching on the arc (%s,%s)", instance.nodeName[u], instance.nodeName[v]);
} */

/**
 * creates and captures the arc marker constraint
 */
SCIP_RETCODE arcmarker::createConsArcMarker(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name,
    SCIP_VAR *arcVar,
    ArcIt& arc,
    CONSTYPE type,
    SCIP_NODE *node)
{

    //std::cout << "Entenring in createConsArcMarker()\n";
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
    SCIP_CALL(createConsData(scip, &consdata, arc, arcVar, type, node));

    //create constraint
    SCIP_CALL(
        SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE,
                       TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

    /* std::cout << "created constraint: ";
    printConsData(consdata);
 */
    return SCIP_OKAY;
}