#include "branch_glcip.h"

ObjBranchruleGLCIP::ObjBranchruleGLCIP(
    SCIP *scip,
    const char *p_name,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z,
    DNodeInfSetsMap &p_inf_set) : ObjBranchrule(scip, p_name, "Defines the branching rule", 50000, -1, 1.0),
                                  instance(p_instance),
                                  x(p_x),
                                  z(p_z),
                                  infSet(p_inf_set)
{
}

ObjBranchruleGLCIP::~ObjBranchruleGLCIP() {}

/* Local methods */
/**
 * get the branching candidates viable for multinode branching
 */
SCIP_RETCODE getBranchCands(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    Arc **arcCands,
    SCIP_VAR **branchCands,     //the address of branching candidates
    SCIP_Real *branchCandsFrac, //pointer to fractionalities of the candidates
    int *nCands                 //number of branching candidates
)
{
    std::cout << "getBranchCands() FUNCTION \n";
    // all arc variables that are in the LP, and have fractional values are viable candidates
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        //SCIP_Bool isFeasibleIntegral = SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a]));
        //std::cout << "var status = " << SCIPvarGetStatus(z[a]) <<std::endl;

        //if (SCIPvarGetStatus(z[a]) == SCIP_VARSTATUS_COLUMN && !isFeasibleIntegral)
        if (!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a])))
        {
            //std::cout << *nCands << " if condition  \n";
            (branchCands)[*nCands] = z[a];
            (branchCandsFrac)[*nCands] = MAX(1 - SCIPvarGetLPSol(z[a]), SCIPvarGetLPSol(z[a]));
            (*nCands)++;
        }
    }

    return SCIP_OKAY;
}

/**
 * branch on a selected binary arc variable
 * TODO decide about the child nodes
 */
SCIP_RETCODE branchOnArcVar(
    SCIP *scip,
    SCIP_VAR **candidates,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    /* SCIP_Real priority;
    SCIP_Real estimate;
    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;

    // check all candidates
    for (int i = 0; i < nCands; i++)
    {
        if (SCIPvarGetStatus(candidates[i]) == SCIP_VARSTATUS_LOOSE ||
            SCIPvarGetStatus(candidates[i]) == SCIP_VARSTATUS_COLUMN &&
                !SCIPisZero(scip, SCIPvarGetLPSol(z[arc])) && !SCIPisEQ(scip, SCIPvarGetSol(z[arc], 1.0)))
        {
            priority = SCIPcalcNodeselPriority(scip, candidates[i], SCIP_BRANCHDIR_UPWARDS, 1.0);
            estimate = SCIPcalcChildEstimate(scip, candidates[i], 1.0);
            tmp = SCIPcalcChildEstimate(scip, candidates[i], 0);
            minestzero = MIN(tmp, minestzero);

            // branch all viable candidates
            SCIP_CALL(SCIPcreateChild(scip, &node, priority, estimate));
            SCIP_CALL(SCIPchgVarLbNode(scip, node, candidates[i], 1.0));

            branched = TRUE;
            *result = SCIP_BRANCHED;
        }
    } */

    std::cout << "branchOnArcVar FUNCTION \n";

    // variables for finding the most fractional column
    double fractionality;
    double bestFractionality;

    int bestCand = -1;

    // search the least fractional candidate
    bestFractionality = 1;
    for (int i = 0; i < nCands; i++)
    {
        assert(candidates[i] != NULL);
        fractionality = branchCandsFrac[i];
        fractionality = min(fractionality, 1 - fractionality);

        if (fractionality < bestFractionality)
        {
            bestFractionality = fractionality;
            bestCand = i;
        }
    }

    assert(bestCand >= 0);
    assert(SCIPisFeasPositive(scip, bestFractionality));

    SCIPinfoMessage(scip, NULL, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, factor=%g)\n",
                    nCands, bestCand, SCIPvarGetName(candidates[bestCand]), branchCandsFrac[bestCand],
                    SCIPvarGetBranchFactor(candidates[bestCand]));

    // perform the branching
    SCIP_CALL(SCIPbranchVar(scip, candidates[bestCand], NULL, NULL, NULL));
    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}

SCIP_RETCODE branchOnArcVar2(
    SCIP *scip,
    SCIP_VAR **candidates,
    Arc **arcCands,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;
    SCIP_CONS *consWith;
    SCIP_CONS *consWithout;

    Arc arc;

    std::cout << "branchOnArcVar2 FUNCTION \n";

    // variables for finding the most fractional column
    double fractionality;
    double bestFractionality;

    int bestCand = -1;

    // search the least fractional candidate
    bestFractionality = 1;
    for (int i = 0; i < nCands; i++)
    {
        assert(candidates[i] != NULL);
        fractionality = branchCandsFrac[i];
        //fractionality = min(fractionality, 1 - fractionality);

        if (fractionality < bestFractionality)
        {
            bestFractionality = fractionality;
            bestCand = i;
        }
    }

    assert(bestCand >= 0);
    assert(SCIPisFeasPositive(scip, bestFractionality));

   /*  SCIPinfoMessage(scip, NULL, " -> %d candidates, 
                selected candidate %d: variable <%s> (frac=%g, factor=%g)\n",
                nCands, bestCand, SCIPvarGetName(candidates[bestCand]),
                branchCandsFrac[bestCand], SCIPvarGetBranchFactor(candidates[bestCand])); */

    // perform the branching
    //SCIP_CALL(SCIPbranchVar(scip, candidates[bestCand], NULL, NULL, NULL));
    // create the branch-and-bound tree child nodes of the current node
    SCIP_CALL( SCIPcreateChild(scip, &leftChild, 0, SCIPgetLocalTransEstimate(scip)) );
    SCIP_CALL( SCIPcreateChild(scip, &rightChild, 0, SCIPgetLocalTransEstimate(scip)) );

    // create corresponding constraints (createConsArcMarker())
    
    SCIP_CALL( ConshdlrArcMarker::createConsArcMarker(scip, &consWithout, "without",
                arcCands[bestCand],  WITHOUT,  leftChild, TRUE) );
    SCIP_CALL( ConshdlrArcMarker::createConsArcMarker(scip, &consWith, "with",
                arc[bestCand],  WITH, rightChild, TRUE) );

    // add constraints to nodes
    SCIP_CALL( SCIPaddConsNode(scip, leftChild, consWithout, NULL) );
    SCIP_CALL( SCIPaddConsNode(scip, rightChild, consWith, NULL) );    

    // release constraints
    SCIP_CALL(SCIPreleaseCons(scip, &consWithout));
    SCIP_CALL(SCIPreleaseCons(scip, &consWith));

    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}

/* SCIP_DECL_BRANCHINIT(ObjBranchruleGLCIP::scip_init)
{
    return SCIP_OKAY;
} */
/**
 * branching execution method for fractional LP solutions
 */
SCIP_DECL_BRANCHEXECLP(ObjBranchruleGLCIP::scip_execlp)
{
    //std::cout << "---------- BRANCHING ----------\n";
    SCIPinfoMessage(scip, NULL, "Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));
    SCIP_VAR **candidates;        // candidates for branching
    Arc **arcCands;
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array

    SCIP_CALL(SCIPallocClearBufferArray(scip, &arcCands, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.m));
    // get branching candidates
    SCIP_CALL(getBranchCands(scip, instance, z, arcCands, candidates, fractionalitiesCands, &nCands));
    assert(nCands > 0);
    *result = SCIP_DIDNOTRUN;

    // perform the branching
    SCIP_CALL(branchOnArcVar2(scip, candidates, fractionalitiesCands, nCands, result));

    // free memory
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}