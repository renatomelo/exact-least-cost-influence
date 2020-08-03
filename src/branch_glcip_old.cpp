#include "branch_glcip.h"
#include "consarcmarker.h"

using namespace arcmarker;

ObjBranchruleGLCIP::ObjBranchruleGLCIP(
    SCIP *scip,
    const char *p_name,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z,
    DNodeInfSetsMap &p_inf_set,
    ArcBoolMap &p_isAble) : ObjBranchrule(scip, p_name, "Defines the branching rule", 50000, -1, 1.0),
                            instance(p_instance),
                            x(p_x),
                            z(p_z),
                            infSet(p_inf_set),
                            isAble(p_isAble)
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
    ArcIt *arcCands,
    SCIP_VAR **branchCands,     //the address of branching candidates
    SCIP_Real *branchCandsFrac, //pointer to fractionalities of the candidates
    int *nCands,                //number of branching candidates
    ArcBoolMap &isAble)
{
    //std::cout << "getBranchCands() FUNCTION \n";
    // all arc variables that are in the LP, and have fractional values are viable candidates
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        if (isAble[a] && !SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a])))
        //if (!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a])))
        {
            (branchCands)[*nCands] = z[a];
            (arcCands)[*nCands] = a;
            (branchCandsFrac)[*nCands] = MAX(1 - SCIPvarGetLPSol(z[a]), SCIPvarGetLPSol(z[a]));
            (*nCands)++;
        }
    }

    return SCIP_OKAY;
}

/**
 * branch on a selected binary arc variable
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

    /* SCIPinfoMessage(scip, NULL, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, factor=%g)\n",
                    nCands, bestCand, SCIPvarGetName(candidates[bestCand]), branchCandsFrac[bestCand],
                    SCIPvarGetBranchFactor(candidates[bestCand])); */

    // perform the branching
    SCIP_CALL(SCIPbranchVar(scip, candidates[bestCand], NULL, NULL, NULL));
    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}

SCIP_RETCODE branchOnArcVar2(
    SCIP *scip,
    SCIP_VAR **candidates,
    ArcIt *arcCands,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;
    SCIP_CONS *consWith;
    SCIP_CONS *consWithout;

    //std::cout << "branchOnArcVar2 FUNCTION \n";

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

        if (fractionality < bestFractionality)
        {
            bestFractionality = fractionality;
            bestCand = i;
        }
    }

    assert(bestCand >= 0);
    assert(SCIPisFeasPositive(scip, bestFractionality));

    ArcIt arc = arcCands[bestCand];

    /* SCIPinfoMessage(scip, NULL, "-> %d candidates, selected candidate: variable <%s> (frac=%g, factor=%g)\n",
                    nCands, SCIPvarGetName(candidates[bestCand]),
                    branchCandsFrac[bestCand], SCIPvarGetBranchFactor(candidates[bestCand])); */

    // perform the branching
    // create the branch-and-bound tree child nodes of the current node
    SCIP_CALL(SCIPcreateChild(scip, &leftChild, 0, SCIPgetLocalTransEstimate(scip)));
    SCIP_CALL(SCIPcreateChild(scip, &rightChild, 0, SCIPgetLocalTransEstimate(scip)));

    //std::cout << "creating the constraint handlers \n";
    // create corresponding constraints
    SCIP_CALL(createConsArcMarker(scip, &consWith, "with", candidates[bestCand],
                                  arc, WITH, rightChild));
    SCIP_CALL(createConsArcMarker(scip, &consWithout, "without", candidates[bestCand],
                                  arc, WITHOUT, leftChild));

    // add constraints to nodes

    SCIP_CALL(SCIPaddConsNode(scip, rightChild, consWith, NULL));
    SCIP_CALL(SCIPaddConsNode(scip, leftChild, consWithout, NULL));

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
/*     cout << "vars able to be branched" << endl;
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        if (isAble[a])
        {
            DNode u = instance.g.source(a);
            DNode v = instance.g.target(a);
            cout << "(" << instance.nodeName[u] << "," << instance.nodeName[v] << ")\n";
        }
    } */
    //std::cout << "---------- BRANCHING ----------\n";
    //SCIPinfoMessage(scip, NULL, "Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));
    SCIP_VAR **candidates;        // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array
    ArcIt *arcCands;

    SCIP_CALL(SCIPallocClearBufferArray(scip, &arcCands, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.m));
    // get branching candidates
    SCIP_CALL(getBranchCands(scip, instance, z, arcCands, candidates, fractionalitiesCands, &nCands, isAble));
    assert(nCands > 0);
    *result = SCIP_DIDNOTRUN;

    // perform the branching
    //SCIP_CALL(branchOnArcVar2(scip, candidates, arcCands, fractionalitiesCands, nCands, result));
    SCIP_CALL(branchOnArcVar(scip, candidates, fractionalitiesCands, nCands, result));


    // free memory
    SCIPfreeBufferArray(scip, &arcCands);
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    //std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}