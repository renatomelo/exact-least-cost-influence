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
    SCIP_Bool added;
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
    (*consdata)->added = FALSE;

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
    //std::cout << "checkVariable() function\n";
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

/* Local methods */
/**
 * get the set of candidate vertices viable for multinode branching based on the constraint
 */
SCIP_RETCODE getBranchCands(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    DNode *branchCands, //the address of branching candidates
    SCIP_Real *branchCandsFrac,
    int *nCands //number of branching candidates
)
{
    std::cout << "getBranchCands() FUNCTION \n";
    // all the vertices in which the sum of incoming arcs' weights is fractional in the LP solution are viable candidates
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //skip the vertex variables fixed in zero
        if (SCIPisEQ(scip, SCIPvarGetUbLocal(x[v]), 0))
            continue;

        double sum = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            sum += SCIPvarGetLPSol(z[a]);
        }

        if (!SCIPisFeasIntegral(scip, sum))
        {
            //cout << "sum of " << instance.nodeName[v] << ": " << sum << endl;
            (branchCands)[*nCands] = v;
            (branchCandsFrac)[*nCands] = MAX(1 - sum, sum);
            (*nCands)++;
        }
    }

    return SCIP_OKAY;
}

/**
 * branch on a selected vertex based on its incoming arcs
 */
SCIP_RETCODE branchOnVertexDegree(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    DNode *candidates,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;
    SCIP_CONS *consUpper;
    SCIP_CONS *consLower;

    std::cout << "branchOnVertexDegree()\n";
    //cout << "number of candidates: " << nCands << endl;

    // search the best candidate
    // the vertex with greatest sum of incoming arcs' weights among the candidates
    double maximum = branchCandsFrac[0];
    int bestCand = 0;
    for (int i = 1; i < nCands; i++)
    {
        if (maximum < branchCandsFrac[i])
        {
            maximum = branchCandsFrac[i];
            bestCand = i;
        }
    }

    //cout << "best candidate: " << instance.nodeName[candidates[bestCand]] << ", value "<< branchCandsFrac[bestCand] << endl;

    // chosen vertex
    DNode v = candidates[bestCand];

    /* cout << "branching on vertex " << instance.nodeName[candidates[bestCand]]
         << " with fractional degree sum = " << branchCandsFrac[bestCand] << endl; */

    double d = branchCandsFrac[bestCand];
    //cout << "d = " << d << endl;

    // perform the branching (If the sum is fractional, create two child nodes. Otherwise, create three child nodes)
    // create the branch-and-bound tree child nodes of the current node
    SCIP_CALL(SCIPcreateChild(scip, &leftChild, 0, SCIPgetLocalTransEstimate(scip)));
    SCIP_CALL(SCIPcreateChild(scip, &rightChild, 0, SCIPgetLocalTransEstimate(scip)));

    cout << "add the constraints here\n";

    //std::cout << "creating the constraint handlers \n";
    // create corresponding constraints
    /* SCIP_CALL(createDegreeCons(scip, &consUpper, "upper", candidates[bestCand], 
                                floor(d), UPPERBOUND, rightChild)); */
    SCIP_CALL(SCIPcreateConsLinear(scip, &consUpper, "upper-degree-cons", 0, NULL, NULL, 0, floor(d),
                                   TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE));

    //add the variables to the new constraint
    int inDegree = 0; // count the number of incoming arcs
    for (InArcIt a(instance.g, v); a != INVALID; ++a)
    {
        SCIPaddCoefLinear(scip, consUpper, z[a], 1);
        inDegree++;
    }
    //SCIPaddCons(scip, consUpper);

    SCIPprintCons(scip, consUpper, NULL);
    SCIPinfoMessage(scip, NULL, "\n");

    /* SCIP_CALL(createDegreeCons(scip, &consLower, "lower", candidates[bestCand], 
                                ceil(d), LOWERBOUND, leftChild)); */

    SCIP_CALL(SCIPcreateConsLinear(scip, &consLower, "lower-degree-cons", 0, NULL, NULL, ceil(d),
                                   inDegree, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE));

    //add the variables to the new constraint
    for (InArcIt a(instance.g, v); a != INVALID; ++a)
    {
        SCIPaddCoefLinear(scip, consLower, z[a], 1);
    }
    //SCIPaddCons(scip, consLower);

    SCIPprintCons(scip, consLower, NULL);
    SCIPinfoMessage(scip, NULL, "\n");
    //exit(0);

    // add constraints to nodes
    SCIP_CALL(SCIPaddConsNode(scip, rightChild, consUpper, NULL));
    SCIP_CALL(SCIPaddConsNode(scip, leftChild, consLower, NULL));

    // release constraints
    SCIP_CALL(SCIPreleaseCons(scip, &consUpper));
    SCIP_CALL(SCIPreleaseCons(scip, &consLower));

    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}

SCIP_RETCODE printDegreeRows(SCIP *scip)
{
    /* get data */
    SCIP_ROW **rows;
    int nrows;

    SCIP_CALL(SCIPgetLPRowsData(scip, &rows, &nrows));
    assert(nrows > 0);
    assert(rows != NULL);
    cout << "number of rows at node " << SCIPnodeGetNumber(SCIPgetFocusNode(scip)) << ": " << nrows << endl;

    string ub = "upper-degree-cons";
    string lb = "lower-degree-cons";
    for (int i = 0; i < nrows; i++)
    {
        if (SCIProwGetName(rows[i]) == ub || SCIProwGetName(rows[i]) == lb)
        {
            SCIPprintRow(scip, rows[i], NULL);
        }
    }

    return SCIP_OKAY;
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

/* SCIP_DECL_CONSENFOLP(DegreeCons::scip_enfolp)
{
    std::cout << "---------- BRANCHING ON CONSENFOLP ----------\n";

    printf("Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", 
            SCIPnodeGetNumber(SCIPgetFocusNode(scip)), SCIPgetFocusDepth(scip));

    DNode *candidates;            // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array

    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.n));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.n));

    // get branching candidates
    SCIP_CALL(getBranchCands(scip, instance, x, z, candidates, fractionalitiesCands, &nCands));
    assert(nCands > 0);

    // *result = SCIP_DIDNOTRUN;
    *result = SCIP_FEASIBLE;

    if (nCands > 0)
        // perform the branching
        SCIP_CALL(branchOnVertexDegree(scip, instance, z, candidates, fractionalitiesCands, nCands, result));

    // free memory
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    printDegreeRows(scip);
    std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
} */

SCIP_DECL_CONSENFOLP(DegreeCons::scip_enfolp)
{
    std::cout << "---------- CONSENFOLP ----------\n";

    cout << "Focus node: " << SCIPnodeGetNumber(SCIPgetFocusNode(scip)) << endl;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(result != NULL);

    SCIP_CONS *cons;
    double ub = SCIPinfinity(scip);
    double lb = 0;

    cout << "nconss = " << nconss << endl;

    SCIP_CONSDATA *consdata;
    consdata = SCIPconsGetData(conss[nconss - 1]);

    if (consdata->added == TRUE)
    {
        cout << "constraint already added!!" << endl;

        //printDegreeRows(scip);

        *result = SCIP_FEASIBLE;
        return SCIP_OKAY;
    }

    string name = "";
    if (consdata->type == UPPERBOUND)
    {
        name = "upper-degree-cons";
        ub = consdata->bound;
        //cout << name << endl;
    }
    else if (consdata->type == LOWERBOUND)
    {
        name = "lower-degree-cons";
        lb = consdata->bound;
        //cout << name << endl;
    }
    else
        name = "fixedsize-degree-cons";

    SCIP_CALL(SCIPcreateConsLinear(scip, &cons, name.c_str(), 0, NULL, NULL, lb, ub,
                                   TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE));

    DNode v = consdata->vertex;
    //add the variables to the new constraint
    for (InArcIt a(instance.g, v); a != INVALID; ++a)
        SCIPaddCoefLinear(scip, cons, z[a], 1);

    SCIPprintCons(scip, cons, NULL);
    SCIPinfoMessage(scip, NULL, "\n");

    consdata->added = TRUE;

    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);

    *result = SCIP_CONSADDED;

    printDegreeRows(scip);
    std::cout << "---------- CONSTRAINT ADDED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOPS(DegreeCons::scip_enfops)
{
    cout << "SCIP_DECL_CONSENFOPS" << endl;
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_DECL_CONSLOCK(DegreeCons::scip_lock)
{
    return SCIP_OKAY;
}

SCIP_DECL_CONSCHECK(DegreeCons::scip_check)
{
    cout << "SCIP_DECL_CONSCHECK" << endl;
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

    //std::cout << "propagation of constraint handler <" << CONSHDLR_NAME << ">\n";

    //get vars and number of vars

    *result = SCIP_DIDNOTFIND;

    for (int i = 0; i < nconss; i++)
    {
        consdata = SCIPconsGetData(conss[i]);
        //cout << "number of propagated variables: " << consdata->nPropagatedVars << endl;

        // chech if all previously generated variables are valid for this constraint
        assert(checkConsData(scip, instance, consdata, TRUE));

        if (!consdata->propagated)
        {
            //std::cout << "propagate constraint <" << SCIPconsGetName(conss[i]) << ">:";
            //printConsData(consdata);

            nVars = infSet[consdata->vertex].size();

            SCIP_CALL(fixVariables(scip, instance, infSet, consdata, nVars, result));
            consdata->nPropagations++;

            if (*result != SCIP_CUTOFF)
            {
                cout << "SCIP_CUTOFF" << endl;
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

    /* std::cout << "activate constraint <" << SCIPconsGetName(cons) << "> on vertex <"
              << instance.nodeName[consdata->vertex] << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata); */

    //std::cout << "number of propagated variables: " << consdata->nPropagatedVars << std::endl;
    if (consdata->nPropagatedVars != (int)infSet[consdata->vertex].size())
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

    /* std::cout << "deactivate constraint <" << SCIPconsGetName(cons) << "> at node <"
              << SCIPnodeGetNumber(consdata->node) << "> in depth <"
              << SCIPnodeGetDepth(consdata->node) << ">: ";
    printConsData(consdata); */

    //set the number of propagated variables to current number of variables in SCIP
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
        SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, TRUE, FALSE,
                       TRUE, TRUE, FALSE, FALSE, FALSE, TRUE));

    //std::cout << "created constraint: ";
    //printConsData(consdata);

    return SCIP_OKAY;
}

SCIP_RETCODE degreecons::createDegreeCons2(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name)
{
    std::cout << "Entenring in createDegreeCons2()\n";

    SCIP_CONSHDLR *conshdlr;
    SCIP_CONSDATA *consdata;

    // find the degree constraint handler
    conshdlr = SCIPfindConshdlr(scip, "degree-constraint-handler");
    if (conshdlr == NULL)
    {
        SCIPerrorMessage("degree constraint handler not found\n");
        return SCIP_PLUGINNOTFOUND;
    }

    //create constraint
    SCIP_CALL(
        SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, TRUE, FALSE,
                       FALSE, TRUE, FALSE, FALSE, FALSE, FALSE));

    //std::cout << "created constraint: ";
    //printConsData(consdata);

    return SCIP_OKAY;
}