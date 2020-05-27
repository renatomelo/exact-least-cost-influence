#include "heur_ordering.h"
#include <queue>
#include <functional> // std::greater

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurOrdering::scip_free)
{
    return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurOrdering::scip_init)
{
    //cout << "SCIP_DECL_HEURINIT\n";

    /* create heuristic data */
    SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
    return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurOrdering::scip_exit)
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
SCIP_DECL_HEURINITSOL(HeurOrdering::scip_initsol)
{
    //cout << "SCIP_DECL_HEURINITSOL\n";

    return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurOrdering::scip_exitsol)
{
    return SCIP_OKAY;
}

void sortVertices(GLCIPInstance &instance)
{
    // Using lambda to compare vertices
    auto cmp = [&](DNode &u, DNode &v) {
        double sum1 = 0;
        for (OutArcIt a(instance.g, u); a != INVALID; ++a)
        {
            sum1 += instance.influence[a];
        }

        double sum2 = 0;
        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
        {
            sum2 += instance.influence[a];
        }
        double dif1 = sum1 - instance.threshold[u];
        double dif2 = sum2 - instance.threshold[v];

        return dif1 < dif2;
    };

    priority_queue<DNode, vector<DNode>, decltype(cmp)> ordering(cmp);

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        cout << "inserting: " << instance.nodeName[v] << endl;
        ordering.push(v);
    }

    while (!ordering.empty())
    {
        DNode v = ordering.top();

        double sum1 = 0;
        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
        {
            sum1 += instance.influence[a];
        }
        double score = sum1 - instance.threshold[v];
        cout << "removing: " << instance.nodeName[v] << ", thr = " << instance.threshold[v]
             << ", influence = " << sum1 << ", score = " << score << endl;
        ordering.pop();
    }
}

SCIP_RETCODE HeurOrdering::constructNewSol(
    SCIP *scip,
    SCIP_SOL *newsol)
{
    double totalCost = 0;
    set<DNode> actives; // start with a empty solution

    // Using lambda to compare vertices
    auto cmp = [&](DNode &u, DNode &v) {
        double sum1 = 0;
        for (OutArcIt a(instance.g, u); a != INVALID; ++a)
        {
            if (SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
            {
                sum1 += instance.influence[a];
            }
        }

        double sum2 = 0;
        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
        {
            if (SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
            {
                sum2 += instance.influence[a];
            }
        }
        double dif1 = sum1 - instance.threshold[u];
        double dif2 = sum2 - instance.threshold[v];

        return dif1 < dif2;
    };

    // list of vertices sorted by the score defined above
    priority_queue<DNode, vector<DNode>, decltype(cmp)> ordering(cmp);

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //cout << "inserting: " << instance.nodeName[v] << endl;
        ordering.push(v);
    }

    while (actives.size() < instance.alpha * instance.n)
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
    }

    //cout << "Total cost: " << totalCost << endl;

    return SCIP_OKAY;
}

/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(HeurOrdering::scip_exec)
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

    //TODO call here the ordering
    constructNewSol(scip, newsol);
    //exit(0);

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