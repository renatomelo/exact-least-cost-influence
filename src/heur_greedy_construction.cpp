#include "heur_greedy_construction.h"
#include <lemon/min_cost_arborescence.h>
#include <queue>

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurGreedyConstruction::scip_free)
{
    return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurGreedyConstruction::scip_init)
{
    //cout << "SCIP_DECL_HEURINIT\n";

    /* create heuristic data */
    SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
    return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurGreedyConstruction::scip_exit)
{
    //cout << "SCIP_DECL_HEUREXIT\n";

    /* free everything which was created in scip_init */
    SCIP_CALL(SCIPfreeSol(scip, &sol));

    return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process 
 * is about to begin)
 *
 * This method is called when the presolving was finished and the branch and bound process is 
 * about to begin. The primal heuristic may use this call to initialize its branch and bound 
 * specific data.
 */
SCIP_DECL_HEURINITSOL(HeurGreedyConstruction::scip_initsol)
{
    //cout << "SCIP_DECL_HEURINITSOL\n";

    return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurGreedyConstruction::scip_exitsol)
{
    return SCIP_OKAY;
}

double getCostOfSolution(
    GLCIPInstance &instance,
    Digraph &graph,
    DNodeDNodeMap &nodeRef,
    ArcArcMap &arcRef)
{
    double totalCost = 0;

    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        double sum = 0;
        for (InArcIt a(graph, v); a != INVALID; ++a)
            sum += instance.influence[arcRef[a]];

        double solCost = 0;
        // assuming that the incentives are sorted in an increasing order
        // uses the first incentive that overcomes the threshold of v
        for (size_t i = 0; i < instance.incentives[nodeRef[v]].size(); i++)
        {
            if (sum + instance.incentives[nodeRef[v]][i] >= instance.threshold[nodeRef[v]])
            {
                solCost = instance.incentives[nodeRef[v]][i];
                break;
            }
        }
        totalCost += solCost;
    }

    return totalCost;
}

//saves the set of edges not in arborescence and add arbitrarely while there is no cycles
void addExternalArcs(
    GLCIPInstance &instance,
    ArcValueMap &w,
    DNode &root)
{
    ArcBoolMap arborescence(instance.g);

    minCostArborescence(instance.g, w, root, arborescence);

    Digraph graph;
    DNodeDNodeMap nodeRef(graph);
    ArcArcMap arcRef(graph);

    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    vector<Arc> externalArcs;

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        if (!arborescence[arcRef[a]])
        {
            //removing arcs that are not in the arborescence
            externalArcs.push_back(a);
            graph.erase(a);
        }
    }
    //add new arcs to the arborescence while does not form cycle
    while (!externalArcs.empty())
    {
        Arc a = externalArcs.back();
        externalArcs.pop_back();

        Arc e = graph.addArc(graph.source(a), graph.target(a));
        if (!dag(graph))
        {
            //cout << "cycle created\n";
            graph.erase(e);
        }
        else
        {
            //cout << "arc added to the arborescence\n";
            arborescence[arcRef[a]] = TRUE;
        }
    }

    double totalCost = getCostOfSolution(instance, graph, nodeRef, arcRef);
    cout << "cost of addExternalArcs() = " << totalCost << endl;
}

//save the set of edges not in arborescence and sort them by the influence weights in decreasing order
void addSortedArcs(
    GLCIPInstance &instance,
    ArcValueMap &w,
    DNode &root)
{
    ArcBoolMap arborescence(instance.g);

    minCostArborescence(instance.g, w, root, arborescence);

    Digraph graph;
    DNodeDNodeMap nodeRef(graph);
    ArcArcMap arcRef(graph);

    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    vector<Arc> externalArcs;

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        if (!arborescence[arcRef[a]])
        {
            //removing arcs that are not in the arborescence
            externalArcs.push_back(a);
            graph.erase(a);
        }
    }

    // Using lambda to compare arcs by the influence weights
    auto compare_arcs = [&](Arc &a, Arc &b) {
        return instance.influence[arcRef[a]] < instance.influence[arcRef[b]];
    };

    sort(externalArcs.begin(), externalArcs.end(), compare_arcs);

    //add new arcs to the arborescence while does not form cycle
    while (!externalArcs.empty())
    {
        Arc a = externalArcs.back();
        externalArcs.pop_back();

        Arc e = graph.addArc(graph.source(a), graph.target(a));
        if (!dag(graph))
        {
            //cout << "cycle created\n";
            graph.erase(e);
        }
        else
        {
            //cout << "arc added to the arborescence\n";
            arborescence[arcRef[a]] = TRUE;
        }
    }

    double totalCost = getCostOfSolution(instance, graph, nodeRef, arcRef);
    cout << "cost of addSortedArcs() = " << totalCost << endl;
}

typedef Digraph::NodeMap<vector<Arc>> DNodeArcsMap;

//sort the vertices by the cost of given incentives in decreasing order
//while there is no cycles add new arcs to the most expensive vertices
//for each vertex, the incoming arcs are sorted in decreasing order by the influence weights
void addArcsBySortedVertices(
    GLCIPInstance &instance,
    ArcValueMap &w,
    DNode &root)
{
    ArcBoolMap arborescence(instance.g);

    minCostArborescence(instance.g, w, root, arborescence);

    Digraph graph;
    DNodeDNodeMap nodeRef(graph);
    ArcArcMap arcRef(graph);

    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    DNodeArcsMap incomingArcs(graph);

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        if (!arborescence[arcRef[a]])
        {
            //removing arcs that are not in the arborescence
            incomingArcs[graph.target(a)].push_back(a);
            graph.erase(a);
        }
    }

    // Using lambda to compare arcs by the influence weights
    /* auto compare_arcs = [&](Arc &a, Arc &b) {
        return instance.influence[arcRef[a]] < instance.influence[arcRef[b]];
    }; */

    /* for (DNodeIt v(graph); v != INVALID; ++v)
    {
        sort(incomingArcs[v].begin(), incomingArcs[v].end(), compare_arcs);

         cout << "weight of external incoming arcs of " << instance.nodeName[nodeRef[v]] << ": ";
        while (!incomingArcs[v].empty())
        {
            Arc a = incomingArcs[v].back();
            incomingArcs[v].pop_back();
            cout << instance.influence[arcRef[a]] << " ";
        }
        
        cout << endl;
    } */

    // Using lambda to compare vertices by the cost of incentives to be given
    auto cmp = [&](DNode &u, DNode &v) {
        double sum1 = 0;
        double sum2 = 0;

        for (InArcIt a(graph, u); a != INVALID; ++a)
            sum1 += instance.influence[arcRef[a]];

        for (InArcIt a(graph, v); a != INVALID; ++a)
            sum2 += instance.influence[arcRef[a]];

        double dif1 = instance.threshold[nodeRef[u]] - sum1;
        double dif2 = instance.threshold[nodeRef[v]] - sum2;

        return dif1 < dif2;
    };

    // list of vertices sorted by the cost of incentive offered
    priority_queue<DNode, vector<DNode>, decltype(cmp)> ordering(cmp);

    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        double sum = 0;
        for (InArcIt a(graph, v); a != INVALID; ++a)
            sum += instance.influence[arcRef[a]];
        double dif = instance.threshold[nodeRef[v]] - sum;

        if (dif > 0 && nodeRef[v] != root)
        {
            //cout << "inserting: " << instance.nodeName[nodeRef[v]] << endl;
            ordering.push(v);
        }
    }

    while (!ordering.empty())
    {
        double sum = 0;
        double diff = 0;
        DNode v = ordering.top();

        //cout << "removing: " << instance.nodeName[nodeRef[v]] << ": " << dif << endl;
        ordering.pop();

        //add a new edge to the most expensive if does not form cycles
        while (!incomingArcs[v].empty())
        {
            Arc a = incomingArcs[v].back();

            incomingArcs[v].pop_back();

            Arc e = graph.addArc(graph.source(a), graph.target(a));
            if (dag(graph))
            {
                //cout << "arc added to the arborescence\n";
                arborescence[arcRef[a]] = TRUE;

                for (InArcIt b(graph, v); b != INVALID; ++b)
                    sum += instance.influence[arcRef[b]];

                diff = instance.threshold[nodeRef[v]] - sum;

                if (diff - instance.influence[arcRef[a]] > 0)
                {
                    ordering.push(v);
                    //cout << "re inserting " << instance.nodeName[nodeRef[v]] << ": " << diff - instance.influence[arcRef[a]] << endl;
                }
                break;
            }
            else
            {
                //cout << "cycle created\n";
                graph.erase(e);
            }
        }
    }

    double totalCost = getCostOfSolution(instance, graph, nodeRef, arcRef);
    cout << "cost of addArcsBySortedVertices() = " << totalCost << endl;
}

void heurOrdering(GLCIPInstance &instance)
{
    // Using lambda to compare vertices
    auto cmp = [&](DNode &u, DNode &v) {
        double sum1 = 0;
        double sum2 = 0;

        for (OutArcIt a(instance.g, u); a != INVALID; ++a)
            sum1 += instance.influence[a];

        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
            sum2 += instance.influence[a];

        double dif1 = sum1 - instance.threshold[u];
        double dif2 = sum2 - instance.threshold[v];

        return dif1 < dif2;
    };

    priority_queue<DNode, vector<DNode>, decltype(cmp)> ordering(cmp);

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //cout << "inserting: " << instance.nodeName[v] << endl;
        ordering.push(v);
    }

    set<DNode> actives; // start with a empty solution
    double totalCost = 0;

    while (!ordering.empty())
    {
        DNode v = ordering.top();
        actives.insert(v);
        ordering.pop();

        assert(v != INVALID);

        double sum = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            DNode u = instance.g.source(a);
            if (actives.count(u))
                sum += instance.influence[a];
        }

        double cost = 0;
        // assuming that the incentives are sorted in an increasing order
        // uses the first incentive that overcomes the threshold of v
        for (size_t i = 0; i < instance.incentives[v].size(); i++)
        {
            if (sum + instance.incentives[v][i] >= instance.threshold[v])
            {
                cost = instance.incentives[v][i];
                break;
            }
        }
        totalCost += cost;
    }
    cout << "cost of heurOrdering() = " << totalCost << endl;
}

double getMinIncentiveNodeLocal(
    GLCIPInstance &instance,
    set<DNode> actives,
    DNode &node)
{
    double minCost = 1e+20;
    set<DNode> activeNeigbors;

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        activeNeigbors.clear();

        if (!actives.count(v))
        {
            for (InArcIt a(instance.g, v); a != INVALID; ++a)
            {
                DNode u = instance.g.source(a);

                if (actives.count(u))
                {
                    activeNeigbors.insert(u);

                    // if the cost is zero it cannot be improved, then stop
                    if (GLCIPBase::costInfluencingSet(instance, v, activeNeigbors) == 0)
                        break;
                }
            }

            double cost = GLCIPBase::costInfluencingSet(instance, v, activeNeigbors);
            if (cost < minCost)
            {
                minCost = cost;
                node = v;

                // if the cost is zero it cannot be improved, then stop
                if (cost == 0)
                    break;
            }
        }
    }

    return minCost;
}

void heurMinIncentive(GLCIPInstance &instance)
{
    set<DNode> actives; // start with a empty solution
    double totalCost = 0;

    //find the vertex of minimum threshold
    double minThr = 1e+20;
    DNode node = INVALID;

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //choose the minimum cost
        if (instance.threshold[v] < minThr)
        {
            minThr = instance.threshold[v];
            node = v;
            //stops if the threshold is 1 because there is no threshold less than 1
            if (minThr == 1)
                break;
        }
    }

    actives.insert(node);
    totalCost += GLCIPBase::cheapestIncentive(instance, node, 0);

    while (actives.size() < instance.alpha * instance.n)
    {
        DNode v = INVALID;
        set<DNode> nodes;
        //double cost;
        totalCost += getMinIncentiveNodeLocal(instance, actives, v);

        actives.insert(v);
    }
    cout << "cost of heurMinIncentive() = " << totalCost << endl;
}

SCIP_RETCODE HeurGreedyConstruction::constructNewSol(
    SCIP *scip,
    SCIP_SOL *newsol)
{
    set<DNode> actives; // start with a empty solution

    //negated weights
    ArcValueMap w(instance.g);
    for (ArcIt a(instance.g); a != INVALID; ++a)
        w[a] = -instance.influence[a];

    //find the maximum spanning arborecence here. As we negated the weiths on the arcs,
    //we get the maximum spanning arborecence instead of the minimum
    double minCost = 0;
    DNode root = INVALID;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        double cost = 0;
        ArcBoolMap arborescence(instance.g);

        cost = minCostArborescence(instance.g, w, v, arborescence);

        if (cost < minCost)
        {
            minCost = cost;
            root = v;
        }
    }

    //cout << "cost = " << minCost << ", root: " << instance.nodeName[root] << endl;

    addArcsBySortedVertices(instance, w, root);

    addExternalArcs(instance, w, root);

    addSortedArcs(instance, w, root);

    heurOrdering(instance);

    heurMinIncentive(instance);

    //TODO set the new solution and get the cost

    return SCIP_OKAY;
}

SCIP_DECL_HEUREXEC(HeurGreedyConstruction::scip_exec)
{
    //cout << "SCIP_DECL_HEUREXEC\n";

    assert(heur != NULL);

    SCIP_Bool success = FALSE;

    SCIP_SOL *newsol = SCIPgetBestSol(scip);

    assert(result != NULL);

    /* since the timing is SCIP_HEURTIMING_AFTERLPNODE, the current node should have an LP */
    assert(SCIPhasCurrentNodeLP(scip));

    *result = SCIP_DIDNOTRUN;

    /* only call heuristic, if an optimal LP solution is at hand */
    if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
        return SCIP_OKAY;

    /* copy the current LP solution to the working solution */
    SCIP_CALL(SCIPlinkLPSol(scip, sol));

    /* allocate local memory */
    SCIP_CALL(SCIPcreateSol(scip, &newsol, heur));

    //call here the greedy construction
    constructNewSol(scip, newsol);
    exit(0);

    // due to construction we already know, that the solution will be feasible
    SCIP_CALL(SCIPtrySol(scip, newsol, TRUE, TRUE, FALSE, FALSE, FALSE, &success));
    if (success)
    {
        //cout << "heur solution feasible !\n";
        *result = SCIP_FOUNDSOL;
    }
    else
    { // the solution total cost is worst than already existent solutions
        *result = SCIP_DIDNOTFIND;
        //cout << "the solution total cost is worst than already existent solutions\n";
    }

    /* free all local memory */
    SCIP_CALL(SCIPfreeSol(scip, &newsol));

    return SCIP_OKAY;
}