#include "heur_greedy_construction.h"
#include <lemon/min_cost_arborescence.h>

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

SCIP_RETCODE HeurGreedyConstruction::constructNewSol(
    SCIP *scip,
    SCIP_SOL *newsol)
{
    double totalCost = 0;
    set<DNode> actives; // start with a empty solution

    //negated weights
    ArcValueMap w(instance.g);
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        w[a] = -instance.influence[a];
    }

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

    cout << "cost = " << minCost << ", root: " << instance.nodeName[root] << endl;

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

    //GraphViewer::ViewGLCIPSupportGraph(instance, graph, "Arborescence", nodeRef);

    totalCost = 0;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        double sum = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            if (arborescence[a])
            {
                sum += instance.influence[a];
            }
        }

        double solCost = 0;
        // assuming that the incentives are sorted in an increasing order
        // uses the first incentive that overcomes the threshold of v
        for (size_t i = 0; i < instance.incentives[v].size(); i++)
        {
            if (sum + instance.incentives[v][i] >= instance.threshold[v])
            {
                solCost = instance.incentives[v][i];
                break;
            }
        }
        totalCost += solCost;
    }

    cout << "Cost of the arborescence: " << totalCost << endl;
    delete[] sorting;

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

    totalCost = 0;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        double sum = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            if (arborescence[a])
            {
                sum += instance.influence[a];
            }
        }

        double solCost = 0;
        // assuming that the incentives are sorted in an increasing order
        // uses the first incentive that overcomes the threshold of v
        for (size_t i = 0; i < instance.incentives[v].size(); i++)
        {
            if (sum + instance.incentives[v][i] >= instance.threshold[v])
            {
                solCost = instance.incentives[v][i];
                break;
            }
        }
        totalCost += solCost;
    }

    cout << "New total cost: " << totalCost << endl;
    //set the new solution and get the cost
    /* while (actives.size() < instance.alpha * instance.n)
    {
        DNode v = ordering.top();
        actives.insert(v);
        ordering.pop();

        assert(v != INVALID);
        SCIP_CALL(SCIPsetSolVal(scip, newsol, x[v], 1.0));

        double sum = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            DNode u = instance.g.source(a);
            if (actives.count(u) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
            {
                sum += instance.influence[a];

                SCIP_CALL(SCIPsetSolVal(scip, newsol, z[a], 1.0));
            }
        }

        double cost = 0;
        // assuming that the incentives are sorted in an increasing order
        // uses the first incentive that overcomes the threshold of v
        for (size_t i = 0; i < instance.incentives[v].size(); i++)
        {
            if (sum + instance.incentives[v][i] >= instance.threshold[v])
            {
                cost = instance.incentives[v][i];
                SCIP_CALL(SCIPsetSolVal(scip, newsol, xip[v][i], 1.0));
                break;
            }
        }
        //cout << "cost of " << instance.nodeName[v] << " is " << cost << endl;
        totalCost += cost;
    } */

    cout << "Total cost: " << totalCost << endl;

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