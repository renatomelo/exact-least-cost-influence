#include "heur_minincentive.h"
#include <queue>

/**
 * The node with minimal incentive to activate it is chosen
 * An LP solution is taken into account by reducing the MinIncentive greedy value
 * by the incentive given to the node in the LP solution. Ties are broken by 
 * selecting the node with maximal number of non-active neighbors on outgoing arcs.
 */
double HeurMinIncentive::getMinIncentiveNode(
    SCIP *scip,
    set<DNode> actives,
    DNode &node)
{
   double minCost = 1e+20;
   set<DNode> activeNeigbors;

   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      activeNeigbors.clear();

      if (!actives.count(v) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, x[v])))
      {
         for (InArcIt a(instance.g, v); a != INVALID; ++a)
         {
            DNode u = instance.g.source(a);

            if (actives.count(u) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
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

/** 
 * Greedy construction heuristic to obtain feasible solutions to warm-start column 
 * generation with an initial set of influensing-set variables:
 * 
 * At each iteration we activate a not yet active node by paying the minimum available 
 * incentive to reach its hurdle, taking into account the current influence coming 
 * from already active neighbors.
 */
SCIP_RETCODE HeurMinIncentive::greedyConstruction(
    SCIP *scip,
    SCIP_SOL *newsol)
{
   set<DNode> actives; // start with a empty solution
   double totalCost = 0;

   //find the vertex of minimum threshold
   double minThr = 1e+20;
   DNode node = INVALID;

   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //choose the minimum cost
      if (instance.threshold[v] < minThr && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, x[v])))
      {
         minThr = instance.threshold[v];
         node = v;
         //stops if the threshold is 1 because there is no threshold less than 1
         if (minThr == 1)
            break;
      }
   }

   int idx = GLCIPBase::getIndexOfChepeastIncentive(instance, node);

   //cout << "setting var: " << SCIPvarGetName(x[node]) << endl;
   SCIP_CALL(SCIPsetSolVal(scip, newsol, x[node], 1.0));

   //fix the incentive variable for the chosen node to start the propagation
   //cout << "setting vars: " << SCIPvarGetName(xip[node][idx]) << endl;
   SCIP_CALL(SCIPsetSolVal(scip, newsol, xip[node][idx], 1.0));
   totalCost += instance.incentives[node][idx];
   
   actives.insert(node);

   while (actives.size() < instance.alpha * instance.n)
   {
      DNode v = INVALID;
      set<DNode> nodes;
      //double cost;
      totalCost += getMinIncentiveNode(scip, actives, v);

      assert(v != INVALID);
      //cout << "setting var: " << SCIPvarGetName(x[v]) << endl;
      SCIP_CALL(SCIPsetSolVal(scip, newsol, x[v], 1.0));

      double sum = 0;
      // save v's influencing-set and activate it
      for (InArcIt a(instance.g, v); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         if (actives.count(u) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
         {
            sum += instance.influence[a];
            //cout << "setting var: " << SCIPvarGetName(z[a]) << endl;
            SCIP_CALL(SCIPsetSolVal(scip, newsol, z[a], 1.0));
         }
      }

      idx = GLCIPBase::getIndexOfChepeastIncentive(instance, v, sum);
      //cout << "setting var: " << SCIPvarGetName(xip[v][idx]) << endl;
      SCIP_CALL(SCIPsetSolVal(scip, newsol, xip[v][idx], 1.0));

      actives.insert(v);
   }

   //cout << "cost of heurMinIncentive() = " << totalCost << endl;

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurMinIncentive::scip_free)
{
   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurMinIncentive::scip_init)
{
   //cout << "SCIP_DECL_HEURINIT\n";

   /* create heuristic data */
   SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurMinIncentive::scip_exit)
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
SCIP_DECL_HEURINITSOL(HeurMinIncentive::scip_initsol)
{
   //cout << "SCIP_DECL_HEURINITSOL\n";
   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurMinIncentive::scip_exitsol)
{
   //cout << "SCIP_DECL_HEUREXITSOL\n";
   return SCIP_OKAY;
}

/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(HeurMinIncentive::scip_exec)
{
   //sortVerticesTest(instance);
   //cout << "SCIP_DECL_HEUREXEC\n";
   assert(heur != NULL);

   SCIP_Bool success = FALSE;

   //SCIP_SOL *newsol = SCIPgetBestSol(scip);
    SCIP_SOL *newsol;

   assert(result != NULL);
   /* since the timing is SCIP_HEURTIMING_AFTERLPNODE, the current node should have an LP */
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
      return SCIP_OKAY;

   /* copy the current LP solution to the working solution */
   SCIP_CALL(SCIPlinkLPSol(scip, sol));

   /* creates a primal solution living in the original problem space, initialized to zero;
   a solution in original space allows to set original variables to values that would be 
   invalid in thetransformed problem due to preprocessing fixings or aggregations */
   SCIP_CALL(SCIPcreateOrigSol(scip, &newsol, heur));

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
      /* cout << "heur solution worst than already existent solutions!\n";
      cout << "NEW SOL\n";
      double totalCost = 0;
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         for (unsigned int i = 0; i < xip[v].size(); i++)
         {
            if (SCIPgetSolVal(scip, newsol, xip[v][i]) > 0)
            {
               cout << SCIPvarGetName(xip[v][i]) << " = " << SCIPgetSolVal(scip, newsol, xip[v][i]) << endl;
               totalCost += instance.incentives[v][i];
            }
         }
      }

      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         if (SCIPgetSolVal(scip, newsol, z[a]) > 0)
            cout << SCIPvarGetName(z[a]) << " = " << SCIPgetSolVal(scip, newsol, z[a]) << endl;
      }

      // get the cost of this solution
      cout << "solution total cost: " << totalCost << endl;
      exit(0); */
   }

   /* free all local memory */
   SCIP_CALL(SCIPfreeSol(scip, &newsol));

   return SCIP_OKAY;
}
