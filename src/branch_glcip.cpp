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
    SCIP_VAR **branchCands,     //the address of branching candidates
                                //    SCIP_Real *branchCandsSol,  //pointer to solution values of the candidates
    SCIP_Real *branchCandsFrac, //pointer to fractionalities of the candidates
    int *nCands                 //number of branching candidates
)
{
    std::cout << "getBranchCands() FUNCIOTN \n";
    //SCIP_Real *branchCandsSol;
    // all arc variables that are in the LP, and have fractional values are viable candidates
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        SCIP_Bool isFeasibleIntegral = SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a]));

        std::cout << "var status = " << SCIPvarGetStatus(z[a]) <<std::endl;

        if (SCIPvarGetStatus(z[a]) == SCIP_VARSTATUS_COLUMN && !isFeasibleIntegral)
        {
            std::cout << *nCands << " if condition  \n";
            (branchCands)[*nCands] = z[a];
            std::cout << "fgdhjfdsafghdsa " << std::endl;
            (branchCandsFrac)[*nCands] = MAX(1 - SCIPvarGetLPSol(z[a]), SCIPvarGetLPSol(z[a]));
            /* (branchCandsSol)[*nCands] = SCIPvarGetLPSol(z[a]);
            (branchCandsFrac)[*nCands] = MAX(1 - (branchCandsSol)[*nCands], (branchCandsSol)[*nCands]); */
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
    /* SCIP_Bool branched;
    SCIP_Real priority;
    SCIP_Real estimate;
    SCIP_Real tmp;
    SCIP_Real minestzero;

    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;
    int nCands = 0;

    SCIP_CALL(SCIPallocClearBufferArray(scip, &branched, 2));

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

    std::cout << "branchOnArcVar FUNCIOTN \n";

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

    /* SCIPinfoMessage(scip, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, obj=%g, factor=%g, score=%g)\n",
                    nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]), lpcandsfrac[bestcand], bestobj,
                    SCIPvarGetBranchFactor(lpcands[bestcand]), bestscore); */

    // perform the branching
    SCIP_CALL(SCIPbranchVar(scip, candidates[bestCand], NULL, NULL, NULL));
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
    std::cout << "---------- BRANCHING ----------\n\n";
    SCIP_VAR **candidates;        // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array

    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.m));
    // get branching candidates
    SCIP_CALL(getBranchCands(scip, instance, z, candidates, fractionalitiesCands, &nCands));
    assert(nCands > 0);
    *result = SCIP_DIDNOTRUN;

    // SCIP_CALL(SCIPgetLPBranchCands(scip,
    //                                &candidates,           /* array of LP branching candidates, or NULL */
    //                                NULL,                  /* array of LP candidate solution values, or NULL */
    //                                &fractionalitiesCands, /* array of LP candidate fractionalities, or NULL */
    //                                NULL,                  /* number of LP branching candidates, or NULL */
    //                                &numLPCands,           /* number of candidates with maximal priority, or NULL */
    //                                NULL));

    // perform the branching
    SCIP_CALL(branchOnArcVar(scip, candidates, fractionalitiesCands, nCands, result));

    // free memory
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}