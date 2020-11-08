#include "binary_branch.h"
//#include "consarcmarker.h"

//using namespace arcmarker;

BinaryBranch::BinaryBranch(
    SCIP *scip,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z
    //DNodeInfSetsMap &p_inf_set
    ) : ObjBranchrule(scip, "branching-rule", "Defines the branching rule", 50000, -1, 1.0),
        instance(p_instance),
        x(p_x),
        z(p_z)
//infSet(p_inf_set)
{
}

BinaryBranch::~BinaryBranch() {}

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
    int *nCands                 //number of branching candidates
)
{
    //std::cout << "getBranchCands() FUNCTION \n";
    // all arc variables that are in the LP, and have fractional values are viable candidates
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        if (!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[a])))
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

// search the least fractional candidate
int leastFractional(
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

void getSubGraphLocally(
    SCIP *scip,
    GLCIPInstance &instance,
    Digraph &graph,
    DNodeDNodeMap &nodeRef,
    ArcArcMap &arcRef,
    ArcSCIPVarMap &z)
{
    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        if (SCIPisEQ(scip, SCIPvarGetUbLocal(z[arcRef[a]]), 0))
        {
            //removing arc variables fixed in zero
            graph.erase(a);
        }
    }

    //GraphViewer::ViewGLCIPSupportGraph(instance, graph, "Sub-graph", nodeRef);
}

/**
 * Based in a fractional solution, look for some unique arc or vertice that separates the graph in more strongly connected components
 */
SCIP_VAR *thereIsUniqueSeparator(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    DNodeSCIPVarMap &x)
{
    //bool thereIs = false;
    int n = 0; //save the maximum number of strongly connected components
    SCIP_VAR *var = NULL;
    vector<ArcIt> arcs;     //saves fractional arc-variables
    vector<DNode> vertices; //saves fractional vertex-variables

    //TODO look for cut arcs
    //TODO for each arc, create a subgraph excluding the fixed variables and the arc a

    //get the support graph of the current feasible solution
    Digraph graph;
    DNodeDNodeMap nodeRef(graph); //save the reference to the original node
    ArcArcMap arcRef(graph);      //save the reference to the original arc
    //getSubGraphLocally(scip, instance, graph, nodeRef, arcRef, z);

    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        //removing arc variables fixed in zero
        if (SCIPisEQ(scip, SCIPvarGetUbLocal(z[arcRef[a]]), 0))
            graph.erase(a);
        //getting fractional arc-variables
        else if (!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(z[arcRef[a]])))
        {
            //cout << SCIPvarGetName(z[arcRef[a]]) << " has value = " << SCIPvarGetLPSol(z[arcRef[a]]) << endl;
            arcs.push_back(a);
        }
    }

    n = countStronglyConnectedComponents(graph);

    for (ArcIt a : arcs)
    {
        DNode u = graph.source(a);
        DNode v = graph.target(a);

        //remove fractional arc a
        graph.erase(a);

        //count the number of strongly connected components in this subgraph without a
        int tmp = countStronglyConnectedComponents(graph);
        if (tmp > n)
        {
            cout << "n = " << n << ", number of new components (arc) = " << tmp << endl;

            //save the variable associated to the most separator arc
            n = tmp;
            var = z[arcRef[a]];
            //thereIs = true;
        }

        //add again the removed arc
        graph.addArc(u, v);
    }

    //look for cut vertices
    if (instance.alpha < 1)
    {
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            //removing vertex variables fixed in zero
            if (SCIPisEQ(scip, SCIPvarGetUbLocal(x[nodeRef[v]]), 0))
                graph.erase(v);
            //getting fractional arc-variables
            else if (!SCIPisFeasIntegral(scip, SCIPvarGetLPSol(x[nodeRef[v]])))
                vertices.push_back(v);
        }

        for (DNode v : vertices)
        {
            //save the incoming and outgoing neighbors to add the arcs again later
            vector<DNode> incoming, outgoing;
            //cout << "incoming neighbors:";
            for (InArcIt a(graph, v); a != INVALID; ++a)
            {
                incoming.push_back(graph.source(a));
                //cout << " " << instance.nodeName[graph.source(a)];
            }
            //cout << endl;

            //cout << "outgoing neighbors:";
            for (OutArcIt a(graph, v); a != INVALID; ++a)
            {
                outgoing.push_back(graph.target(a));
                //cout << " " << instance.nodeName[graph.target(a)];
            }
            //cout << endl;

            graph.erase(v);

            int tmp = countStronglyConnectedComponents(graph);
            if (tmp > n)
            {
                cout << "n = " << n << ", number of new components (vertex) = " << tmp << endl;

                //save the variable associated to the most separator arc
                n = tmp;
                var = x[nodeRef[v]];
            }

            //add the node and the arcs again to the graph
            DNode v2 = graph.addNode();
            for (DNode u : incoming)
                graph.addArc(u, v2);
            for (DNode w : outgoing)
                graph.addArc(v2, w);
        }
    }

    return var;
}

/**
 * branch on a selected binary arc variable
 */
SCIP_RETCODE branchOnArcVar(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    DNodeSCIPVarMap &x,
    ArcIt *arcCands,
    SCIP_VAR **candidates,
    SCIP_Real *branchCandsFrac,
    int nCands,
    SCIP_RESULT *result //pointer to store result of branching
)
{
    SCIP_VAR *var = thereIsUniqueSeparator(scip, instance, z, x);
    if (!var)
    {
        //cout << "there is no unique separator\n";
        if (!nCands)
            return SCIP_OKAY;
        //if we don't find a arc or vertex cut we choose the least fractional variable
        int best = leastFractional(scip, candidates, branchCandsFrac, nCands);
        var = candidates[best];
    }
    else
        cout << SCIPvarGetName(var) << endl;
    /* SCIPinfoMessage(scip, NULL, " -> %d candidates, selected candidate %d: variable <%s> (frac=%g, factor=%g)\n",
                    nCands, bestCand, SCIPvarGetName(candidates[bestCand]), branchCandsFrac[bestCand],
                    SCIPvarGetBranchFactor(candidates[bestCand])); */

    // perform the branching
    SCIP_CALL(SCIPbranchVar(scip, var, NULL, NULL, NULL));
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
SCIP_DECL_BRANCHEXECLP(BinaryBranch::scip_execlp)
{
    std::cout << "---------- BRANCHING ----------\n";
    //SCIPinfoMessage(scip, NULL, "Start branching at node %" SCIP_LONGINT_FORMAT ", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));
    SCIP_VAR **candidates;        // candidates for branching
    double *fractionalitiesCands; // fractionalities of candidates
    int nCands = 0;               // length of array
    ArcIt *arcCands;

    SCIP_CALL(SCIPallocClearBufferArray(scip, &arcCands, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &candidates, instance.m));
    SCIP_CALL(SCIPallocClearBufferArray(scip, &fractionalitiesCands, instance.m));
    // get branching candidates
    SCIP_CALL(getBranchCands(scip, instance, z, arcCands, candidates, fractionalitiesCands, &nCands));
    assert(nCands > 0);
    *result = SCIP_DIDNOTRUN;

    // perform the branching
    //SCIP_CALL(branchOnArcVar2(scip, candidates, arcCands, fractionalitiesCands, nCands, result));
    SCIP_CALL(branchOnArcVar(scip, instance, z, x, arcCands, candidates, fractionalitiesCands, nCands, result));

    // free memory
    SCIPfreeBufferArray(scip, &arcCands);
    SCIPfreeBufferArray(scip, &candidates);
    SCIPfreeBufferArray(scip, &fractionalitiesCands);

    //std::cout << "---------- BRANCHED SUCCESFULLY ----------\n\n";
    return SCIP_OKAY;
}