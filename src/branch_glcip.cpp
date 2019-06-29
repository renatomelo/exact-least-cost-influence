#include "branch_glcip.h"

ObjBranchruleGLCIP::ObjBranchruleGLCIP(
    SCIP*                   scip,
    const char*             p_name,
    GLCIPInstance&          p_instance,
    DNodeInfSetsMap&        p_inf_set
    ) : ObjBranchrule(scip, p_name, "Defines the branching rule", 50000, -1, 1.0),
        instance(p_instance),
        infSet(p_inf_set)
{
}

ObjBranchruleGLCIP::~ObjBranchruleGLCIP() {}

/* SCIP_DECL_BRANCHINIT(ObjBranchruleGLCIP::scip_init)
{
    return SCIP_OKAY;
} */
/**
 * branching execution method for fractional LP solutions
 */
SCIP_DECL_BRANCHEXECLP(ObjBranchruleGLCIP::scip_execlp)
{
    //std::cout << "---------- BRANCHING ----------\n\n";
    SCIP_VAR **candidates;        // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int numLPCands;               // length of array

    // variables for finding the most fractional column
    double fractionality;
    double bestFractionality;

    int bestCand;

    // array of variables
    SCIP_VAR** vars;
    int nvars;

    // the 2 influencing

    // nodes of branch-and-bound-tree
    //SCIP_NODE lft, rgt;

    *result = SCIP_DIDNOTRUN;

    // get branching candidates
    SCIP_CALL(SCIPgetLPBranchCands(scip,
                                   &candidates,           /* array of LP branching candidates, or NULL */
                                   NULL,                  /* array of LP candidate solution values, or NULL */
                                   &fractionalitiesCands, /* array of LP candidate fractionalities, or NULL */
                                   NULL,                  /* number of LP branching candidates, or NULL */
                                   &numLPCands,           /* number of candidates with maximal priority, or NULL */
                                   NULL));
    assert(numLPCands > 0);

    bestCand = -1;

    // search the least fractional candidate
    bestFractionality = 1;
    for (int i = 0; i < numLPCands; i++)
    {
        assert(candidates[i] != NULL);
        fractionality = fractionalitiesCands[i];
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

    //std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}