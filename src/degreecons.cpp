#include "degreecons.h"
#include "pricer_glcip.h"

using namespace degreecons;

struct SCIP_ConsData
{
    DNode vertex;
    int bound;                   // the right-hand side of the constraint \sum_{a \in N} z[a] <= d
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
    DNode &vertex,
    int bound,
    CONSTYPE type,   //stores whether the arc have to be in the solution or not
    SCIP_NODE *node) //the node in the B&B-tree at which the cons is sticking
{
    assert(scip != NULL);
    assert(consdata != NULL);
    assert(vertex != NULL);
    assert(type == LOWERBOUND || type == UPPERBOUND || type = FIXEDSIZE);

    SCIP_CALL(SCIPallocBlockMemory(scip, consdata));
    (*consdata)->vertex = vertex;
    (*consdata)->bound = bound;
    (*consdata)->type = type;
    (*consdata)->nPropagatedVars = 0;
    (*consdata)->nPropagations = 0;
    (*consdata)->propagated = FALSE;
    (*consdata)->node = node;

    return SCIP_OKAY;
}

void printConsData(SCIP_CONSDATA *consdata)
{
    string type = "";
    if (consdata->type == UPPERBOUND)
        type = "UPPERBOUND";
    else if (consdata->type == LOWERBOUND)
        type = "LOWERBOUND";
    else
        type = "FIXEDSIZE";

    /*  std::cout << "branching on vertex  (" << consdata->vertexName << ") at node "
               << SCIPnodeGetNumber(consdata->node) << ", constraint of type "
               << type << std::endl; */
    std::cout << "branching at node " << SCIPnodeGetNumber(consdata->node) << ", constraint of type "
              << type << std::endl;
}

/** 
 * fixes a variable to zero if the corresponding influencing-sets are not valid for this 
 * constraint/node (due to branching) 
 */
SCIP_RETCODE checkVariable(
    SCIP *scip,
    GLCIPInstance &instance,
    SCIP_CONSDATA *consdata, // constraint data
    InfluencingSet &ifs,     // variables to check
    int *nFixedVars,         // pointer to store the number of fixed variables
    SCIP_Bool *cuttoff       // pointer to store if a cutoff was detected
)
{
    std::cout << "checkVariable() function\n";
    SCIP_Bool fixed;
    SCIP_Bool infeasible;

    assert(scip != NULL);
    assert(consdata != NULL);
    assert(ifs != NULL);
    assert(nFixedVars != NULL);
    assert(cuttoff != NULL);

    //cout << "checking variable " << SCIPvarGetName(ifs.getVar()) << endl;
    // if variables is locally fixed to zero continue
    if (SCIPvarGetUbLocal(ifs.getVar()) < 0.5)
        return SCIP_OKAY;

    //check if the influencing-set which corresponds to the varialbe is feasible for this constraint
    int size = ifs.getNodes().size();
    int bound = consdata->bound;
    CONSTYPE type = consdata->type;

    if ((type == LOWERBOUND && size < bound) ||
        (type == UPPERBOUND && size > bound) ||
        (type == FIXEDSIZE && size != bound))
    {
        SCIP_CALL(SCIPfixVar(scip, ifs.getVar(), 0.0, &infeasible, &fixed));
        //std::cout << "variable " << SCIPvarGetName(ifs.getVar()) << " fixed in 0\n";

        if (infeasible)
        {
            std::cout << "IS INFEASIBLE\n";
            assert(SCIPvarGetLbLocal(ifs.getVar()) > 0.5);
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
    GLCIPInstance &instance,
    DNodeInfSetsMap &infSet,
    SCIP_CONSDATA *consdata,
    int nVars,
    SCIP_RESULT *result)
{
    int nFixedVars = 0;
    SCIP_Bool cutoff = FALSE;

    DNode v = consdata->vertex;

    printf("check variables %d to %d\n", consdata->nPropagatedVars, nVars);

    for (int i = consdata->nPropagatedVars; i < nVars; i++)
        SCIP_CALL(checkVariable(scip, instance, consdata, infSet[v][i], &nFixedVars, &cutoff));

    printf("fixed %d variables locally\n", nFixedVars);

    if (cutoff)
        *result = SCIP_CUTOFF;
    else if (nFixedVars > 0)
        *result = SCIP_REDUCEDDOM;

    return SCIP_OKAY;
}

/** 
 * check if all variables are valid for the given consdata 
 * (it need to be done only in case of debug) 
 */
SCIP_Bool checkConsData(
    SCIP *scip,              /**< SCIP data structure */
    GLCIPInstance &instance, /**< problem data */
    SCIP_CONSDATA *consdata, /**< constraint data */
    SCIP_Bool beforeprop     /**< is this check performed before propagation? */
)
{
    // TODO get the influencing set variables of the vertex in the consdata (?)
    // iterate over the influencing set variables
    // if variables is locally fixed to zero continue
    // check if the size of influencing-sets which corresponds to the variables is
    // feasible for this constraints
    return TRUE;
}

/**
 * Frees all data of a constraint
 */
SCIP_DECL_CONSDELETE(DegreeCons::scip_delete)
{
    assert(consdata != NULL);

    SCIPfreeBlockMemory(scip, consdata);

    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOLP(DegreeCons::scip_enfolp)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOPS(DegreeCons::scip_enfops)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_DECL_CONSLOCK(DegreeCons::scip_lock)
{
    return SCIP_OKAY;
}

SCIP_DECL_CONSCHECK(DegreeCons::scip_check)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

/**
 * transforms constraint data into data belonging to the trasformed problem
 */
SCIP_DECL_CONSTRANS(DegreeCons::scip_trans)
{
    cout << "SCIP_DECL_CONSTRANS" << endl;
    SCIP_CONSDATA *sourcedata;
    SCIP_CONSDATA *targetdata;

    sourcedata = SCIPconsGetData(sourcecons);
    assert(sourcedata != NULL);

    //create constraint data for the target constraint
    SCIP_CALL(
        createConsData(scip, &targetdata, sourcedata->vertex, sourcedata->bound,
                       sourcedata->type, sourcedata->node));

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
SCIP_DECL_CONSPROP(DegreeCons::scip_prop)
{
    cout << "SCIP_DECL_CONSPROP" << endl;
    SCIP_CONSDATA *consdata;
    int nVars = 0;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(result != NULL);

    std::cout << "propagation of constraint handler <" << CONSHDLR_NAME << ">\n";

    //get vars and number of vars

    *result = SCIP_DIDNOTFIND;

    for (int i = 0; i < nconss; i++)
    {
        consdata = SCIPconsGetData(conss[i]);
        std::cout << "number of propagated variables: " << consdata->nPropagatedVars
                  << std::endl;

        // chech if all previously generated variables are valid for this constraint
        assert(checkConsData(scip, instance, consdata, TRUE));

        // remember of verify about the NDEBUG flag of SCIP
        if (!consdata->propagated)
        {
            std::cout << "propagate constraint <" << SCIPconsGetName(conss[i]) << ">:";
            printConsData(consdata);

            nVars = infSet[consdata->vertex].size();

            SCIP_CALL(fixVariables(scip, instance, infSet, consdata, nVars, result));
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
        assert(checkConsData(scip, instance, consdata, FALSE));
    }

    return SCIP_OKAY;
}

/**
 * Constraint activation notification method of constraint handler.
 * This method is always called when a node is entered on which the constraint
 * has been added. Here, we apply the changes to the pricing data structures.
 */
SCIP_DECL_CONSACTIVE(DegreeCons::scip_active)
{
    cout << "SCIP_DECL_CONSACTIVE" << endl;
    SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->nPropagatedVars <= infSet[consdata->vertex].size());

    std::cout << "activate constraint <" << SCIPconsGetName(cons) << "> on vertex <"
              << instance.nodeName[consdata->vertex] << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata);

    //std::cout << "number of propagated variables: " << consdata->nPropagatedVars << std::endl;
    if (consdata->nPropagatedVars != (int) infSet[consdata->vertex].size())
    {
        SCIPinfoMessage(scip, NULL, "-> mark constraint to be repropagated\n");
        consdata->propagated = FALSE;
        SCIP_CALL(SCIPrepropagateNode(scip, consdata->node));
    }

    return SCIP_OKAY;
}

/**
 * The CONSDEACTIVE method will be called if the node is left again. Since the
 * CONSACTIVE and CONSDEACTIVE methods of different constraints are always called in a
 * stack-like fashion, this should be exactly what you need.
 */
SCIP_DECL_CONSDEACTIVE(DegreeCons::scip_deactive)
{
    cout << "SCIP_DECL_CONSDEACTIVE" << endl;
    SCIP_CONSDATA *consdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

    std::cout << "deactivate constraint <" << SCIPconsGetName(cons) << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata);

    /* set the number of propagated variables to current number of variables in SCIP */
    consdata->nPropagatedVars = infSet[consdata->vertex].size();

    return SCIP_OKAY;
}

/**
 * variable deletion method of constraint handler
 */
SCIP_DECL_CONSDELVARS(DegreeCons::scip_delvars)
{
    //cout << "SCIP_DECL_CONSDELVARS" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_CONSPRINT(DegreeCons::scip_print)
{
    SCIP_CONSDATA *consdata;

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    printConsData(consdata);

    return SCIP_OKAY;
}

/**
 * creates and captures the degree constraint
 */
SCIP_RETCODE degreecons::createDegreeCons(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name,
    DNode &vertex,
    int bound,
    CONSTYPE type,
    SCIP_NODE *node)
{
    std::cout << "Entenring in createDegreeCons()\n";

    SCIP_CONSHDLR *conshdlr;
    SCIP_CONSDATA *consdata;

    // find the degree constraint handler
    conshdlr = SCIPfindConshdlr(scip, "degree-constraint-handler");
    if (conshdlr == NULL)
    {
        SCIPerrorMessage("degree constraint handler not found\n");
        return SCIP_PLUGINNOTFOUND;
    }

    // create constraint data
    SCIP_CALL(createConsData(scip, &consdata, vertex, bound, type, node));

    //create constraint
    SCIP_CALL(
        SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE,
                       TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

    //std::cout << "created constraint: ";
    //printConsData(consdata);

    return SCIP_OKAY;
}