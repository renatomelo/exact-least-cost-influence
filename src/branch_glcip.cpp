#include "branch_glcip.h"
#include "degreecons.h"

using namespace degreecons;

ObjBranchruleGLCIP::ObjBranchruleGLCIP(
    SCIP *scip,
    const char *p_name,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z) : ObjBranchrule(scip,
                                          p_name,                      
                                          "Defines the branching rule", // description
                                          50000,    // priority
                                          3,       // maximal depth level (-1 for no limit)
                                          1.0),     //maximal relative distance from current node's dual bound to primal bound
                            instance(p_instance),
                            x(p_x),
                            z(p_z)
{
}

ObjBranchruleGLCIP::~ObjBranchruleGLCIP() {}

/* Local methods */
/**
 * get the set of candidate vertices viable for multinode branching based on the constraint
 */
SCIP_RETCODE getBranchCands2(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap& x,
    ArcSCIPVarMap &z,
    DNode *branchCands,     //the address of branching candidates
    SCIP_Real *branchCandsFrac,
    int *nCands                 //number of branching candidates
    )
{
    //std::cout << "getBranchCands() FUNCTION \n";
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
SCIP_RETCODE branchOnVertexDegree2(
    SCIP *scip,
    GLCIPInstance &instance,
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

    //std::cout << "branchOnVertexDegree()\n";
    //cout << "number of candidates: " << nCands << endl;

    // search the best candidate
    // the vertex with greatest sum of incoming arcs' weights among the candidates
    double maximum = branchCandsFrac[0];
    int bestCand = 0;
    for (int i = 1; i < nCands; i++)
    {
        if(maximum < branchCandsFrac[i])
        {
            maximum = branchCandsFrac[i];
            bestCand = i;
        }
    }

   //cout << "best candidate: " << instance.nodeName[candidates[bestCand]] << ", value "<< branchCandsFrac[bestCand] << endl;

    // chosen vertex
    //DNode vertex = candidates[bestCand];

    cout << "branching on vertex " << instance.nodeName[candidates[bestCand]] 
         << " with fractional degree sum = " << branchCandsFrac[bestCand] << endl;
    
    double d = branchCandsFrac[bestCand];

    // perform the branching (If the sum is fractional, create two child nodes. Otherwise, create three child nodes)
    // create the branch-and-bound tree child nodes of the current node
    SCIP_CALL(SCIPcreateChild(scip, &leftChild, 0, SCIPgetLocalTransEstimate(scip)));
    SCIP_CALL(SCIPcreateChild(scip, &rightChild, 0, SCIPgetLocalTransEstimate(scip)));

    //std::cout << "creating the constraint handlers \n";
    // create corresponding constraints
    SCIP_CALL(createDegreeCons(scip, &consUpper, "upper", candidates[bestCand], 
                                floor(d), UPPERBOUND, rightChild));
    SCIP_CALL(createDegreeCons(scip, &consLower, "lower", candidates[bestCand], 
                                ceil(d), LOWERBOUND, leftChild));

    // add constraints to nodes

    SCIP_CALL(SCIPaddConsNode(scip, rightChild, consUpper, NULL));
    SCIP_CALL(SCIPaddConsNode(scip, leftChild, consLower, NULL));

    // release constraints
    SCIP_CALL(SCIPreleaseCons(scip, &consUpper));
    SCIP_CALL(SCIPreleaseCons(scip, &consLower));

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
    std::cout << "---------- BRANCHING ----------\n";

    printf("Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", 
            SCIPnodeGetNumber(SCIPgetFocusNode(scip)), SCIPgetFocusDepth(scip));
            
    DNode *candidates;        // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array

    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.n));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.n));

    // get branching candidates
    SCIP_CALL(getBranchCands2(scip, instance, x, z, candidates, fractionalitiesCands, &nCands));
    assert(nCands > 0);
    
    *result = SCIP_DIDNOTRUN;

    if (nCands > 0)
        // perform the branching
        SCIP_CALL(branchOnVertexDegree2(scip, instance, candidates, fractionalitiesCands, nCands, result));

    // free memory
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    /* TODO: use the branching rule and the degreecons, add the constraint in fact at the enfolp of degreecons (or just a row) */
    return SCIP_OKAY;
}