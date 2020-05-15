#include "heur_ordering.h"
#include <queue>
#include <functional>     // std::greater

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

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
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

struct compare
{
    GLCIPInstance &instance;
    compare(GLCIPInstance &_instance) : instance(_instance) {}

    bool operator()(const DNode &u, const DNode &v)
    {
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
    }
};

void sortVerticesTest(GLCIPInstance &instance)
{
   compare cmp(instance);
   priority_queue<DNode, vector<DNode>, cmp> q;
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cout << "inserting: " << instance.nodeName[v] << endl;
      q.push(v);
   }

   while (!q.empty())
   {
      cout << "removing: " << instance.nodeName[q.top()] << endl;
      q.pop();
   }
}

/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(HeurOrdering::scip_exec)
{
    sortVerticesTest(instance);
    exit(0);
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

    //call here the ordering
    greedyConstruction(scip, newsol);

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
        //cout << "heur solution not found !\n";
    }

    /* free all local memory */
    SCIP_CALL(SCIPfreeSol(scip, &newsol));
    //cout << "free all local memory \n";

    return SCIP_OKAY;
}