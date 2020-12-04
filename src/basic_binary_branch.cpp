#include "basic_binary_branch.h"
//#include "consarcmarker.h"

//using namespace arcmarker;

BasicBinaryBranch::BasicBinaryBranch(
    SCIP *scip,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z
    //DNodeInfSetsMap &p_inf_set
    ) : ObjBranchrule(scip, "basic-binary-branching", "Defines the basic 0-1 branching rule", 50000, -1, 1.0),
        instance(p_instance),
        x(p_x),
        z(p_z)
//infSet(p_inf_set)
{
}

BasicBinaryBranch::~BasicBinaryBranch() {}

// search the least fractional candidate
int BasicBinaryBranch::leastFractional(
    SCIP *scip,
    SCIP_VAR **candidates,
    SCIP_Real *branchCandsFrac,
    int nCands)
{
    double fractionality;
    double bestFractionality;

    int bestCand = -1;

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

    return bestCand;
}

/* SCIP_DECL_BRANCHINIT(ObjBranchruleGLCIP::scip_init)
{
    return SCIP_OKAY;
} */
/**
 * branching execution method for fractional LP solutions
 */
SCIP_DECL_BRANCHEXECLP(BasicBinaryBranch::scip_execlp)
{
    //std::cout << "---------- BRANCHING ----------\n";
    //SCIPinfoMessage(scip, NULL, "Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));

    SCIP_VAR **candidates;   // candidates for branching
    double *fractionalities; // fractionalities of candidates
    int nCands = 0;          // length of array
    //ArcIt *arcCands;
    double nvars = SCIPgetNBinVars(scip);

    *result = SCIP_DIDNOTRUN;

    SCIP_VAR *var;

    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, nvars));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalities, nvars));

    // get branching candidates
    SCIPgetLPBranchCands(scip, &candidates, NULL, &fractionalities, &nCands, NULL, NULL);
    assert(nCands > 0);

    if (!nCands)
        return SCIP_OKAY;

    //if we don't find a arc or vertex cut we choose the least fractional variable
    int best = leastFractional(scip, candidates, fractionalities, nCands);
    var = candidates[best];

    // free memory
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalities);

    // perform the branching
    SCIP_CALL(SCIPbranchVar(scip, var, NULL, NULL, NULL));
    *result = SCIP_BRANCHED;

    //std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}