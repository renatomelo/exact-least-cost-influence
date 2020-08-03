#include "branch_glcip.h"
#include "consarcmarker.h"

using namespace arcmarker;

ObjBranchruleGLCIP::ObjBranchruleGLCIP(
    SCIP *scip,
    const char *p_name,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z) : ObjBranchrule(scip,
                                          p_name,                      
                                          "Defines the branching rule", // description
                                          50000,    // priority
                                          -1,       // maximal depth level (-1 for no limit)
                                          1.0),     //maximal relative distance from current node's dual bound to primal bound
                            instance(p_instance),
                            x(p_x),
                            z(p_z)
{
}

ObjBranchruleGLCIP::~ObjBranchruleGLCIP() {}

/* Local methods */
/**
 * get the set of candidates vertices viable for multinode branching based on the constraint
 */
SCIP_RETCODE getBranchCands(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    DNode **branchCands,     //the address of branching candidates
    int *nCands                 //number of branching candidates
    )
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
 * branch on a selected vertex based on its incoming arcs
 */
SCIP_RETCODE branchOnVertexDegree(
    SCIP *scip,
    DNode **candidates,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    SCIP_NODE *leftChild;
    SCIP_NODE *rightChild;
    SCIP_CONS *consWith;
    SCIP_CONS *consWithout;

    //std::cout << "branchOnVertexDegree()\n";

    // search the best candidate (criterion to be defined)

    // chosen vertex
    DNode best;

    /* SCIPinfoMessage(scip, NULL, "-> %d candidates, selected candidate: variable <%s> (frac=%g, factor=%g)\n",
                    nCands, SCIPvarGetName(candidates[bestCand]),
                    branchCandsFrac[bestCand], SCIPvarGetBranchFactor(candidates[bestCand])); */

    // perform the branching (If the sum is fractional, create two child nodes. Otherwise, create three child nodes)
    // create the branch-and-bound tree child nodes of the current node
    SCIP_CALL(SCIPcreateChild(scip, &leftChild, 0, SCIPgetLocalTransEstimate(scip)));
    SCIP_CALL(SCIPcreateChild(scip, &rightChild, 0, SCIPgetLocalTransEstimate(scip)));

    //std::cout << "creating the constraint handlers \n";
    // create corresponding constraints
    SCIP_CALL(createConsArcMarker(scip, &consWith, "upper", candidates[bestCand],
                                  best, UPPERBOUND, rightChild));
    SCIP_CALL(createConsArcMarker(scip, &consWithout, "lower", candidates[bestCand],
                                  best, LOWERBOUND, leftChild));

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