#include "heur_mininfluence.h"

/**
 * The node with minimal incentive to activate it is chosen
 * An LP solution is taken into account by reducing the MinIncentive greedy value
 * by the incentive given to the node in the LP solution. Ties are broken by 
 * selecting the node with maximal number of non-active neighbors on outgoing arcs.
 */
double HeurMinInfluence::getMinIncentiveNode(
    SCIP *scip,
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
//set<DNode> HeurMinInfluence::greedyConstruction(
SCIP_RETCODE HeurMinInfluence::greedyConstruction(
    SCIP *scip,
    SCIP_SOL *newsol)
{
   set<DNode> actives; // start with a empty solution

   //find the cheapest vertex with Lambda_v_empty in the LP solution
   /* double minC = 1e+20;
   DNode node = INVALID;
   unsigned int idx = -1;
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (infSet[v][i].nodes.empty() && SCIPgetSolVal(scip, sol, infSet[v][i].var) > 0)
         {
            //choose the minimum cost
            if (infSet[v][i].cost < minC)
            {
               minC = infSet[v][i].cost;
               node = v;
               idx = i;
            }
            //since we fond the variable Lambda_v_empty we don't need to visit the rest
            break; //break the inner loop
         }
      }
   }

   assert(idx != -1);
   //fix the variable Lambda_v_empty for the chosen node to start the propagation
   SCIP_CALL(SCIPsetSolVal(scip, newsol, infSet[node][idx].var, 1.0));
   SCIP_CALL(SCIPsetSolVal(scip, newsol, x[node], 1.0));
   actives.insert(node);
   cout << "vertex " << instance.nodeName[node] << " chosen to start the propagation\n"; */

   /* while (actives.size() < instance.alpha * instance.n)
   {
      DNode v;
      double minCost = getMinIncentiveNode(scip, actives, v);
      // save v's influencing-set and activate it
      InfluencingSet ifs;
      for (InArcIt a(instance.g, v); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         if (actives.count(u))
         {
            ifs.nodes.insert(u);
            // if the minCost was achieved we don't need to increase the influencing-set more
            if (minCost >= GLCIPBase::costInfluencingSet(instance, v, ifs.nodes))
               break;
         }
      }
      ifs.cost = minCost;
      infSet[v].push_back(ifs);
      actives.insert(v);

      //std::cout << instance.nodeName[v] << " was activated " << std::endl;
   } */

   while (actives.size() < instance.alpha * instance.n)
   {
      DNode v = INVALID;
      set<DNode> nodes;
      //double cost;
      double minCost = getMinIncentiveNode(scip, actives, v);

      assert(v != INVALID);
      SCIP_CALL(SCIPsetSolVal(scip, newsol, x[v], 1.0));
      //cout << "setting variable " << SCIPvarGetName(x[v]) << endl;
      double sum = 0;
      // save v's influencing-set and activate it
      for (InArcIt a(instance.g, v); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         if (actives.count(u) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
         {
            sum += instance.influence[a];
            //cout << "influence :" << instance.influence[a] << endl;
            nodes.insert(u);
            SCIP_CALL(SCIPsetSolVal(scip, newsol, z[a], 1.0));
            //cout << "setting variable " << SCIPvarGetName(z[a]) << endl;

            // if the minCost was achieved we don't need to increase the influencing-set more
            if (minCost >= GLCIPBase::costInfluencingSet(instance, v, nodes))
               break;
         }
      }
      //cost = minCost;

      double slack = (sum + minCost) - instance.threshold[v];

      //if influencing-set 'nodes' isn't minimal remove surplus elements
      for (auto u = nodes.begin(); u != nodes.end();)
      {
         Arc a = findArc(instance.g, *u, v);
         assert(a != INVALID);

         if (instance.influence[a] < slack)
         {
            //cout << "set isn't minimal, removing vertex " << instance.nodeName[*u] << endl;
            u = nodes.erase(u);
            slack += instance.influence[a];
            SCIP_CALL(SCIPsetSolVal(scip, newsol, z[a], 0.0));
         }
         else
            ++u;
      }

      //find the variable related to the influencing-set found
      bool varFound = FALSE;
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (nodes == infSet[v][i].nodes)
         {
            //cout << "var found!" << endl;
            SCIP_CALL(SCIPsetSolVal(scip, newsol, infSet[v][i].var, 1.0));
            //cout << "setting variable " << SCIPvarGetName(infSet[v][i].var) << endl;
            varFound = TRUE;
            break;
         }
      }

      //in this case nodes is not minimal
      if (!varFound)
      {
         cout << "var not found\n";
         cout << "sum of influence = " << sum << ", threshold = "
              << instance.threshold[v] << endl;
      }

      actives.insert(v);

      //std::cout << instance.nodeName[v] << " was activated " << std::endl;
   }

   //showActivatedNodes(instance, actives, infSet);
   //return actives;
   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurMinInfluence::scip_free)
{
   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurMinInfluence::scip_init)
{
   cout << "SCIP_DECL_HEURINIT\n";

   /* create heuristic data */
   SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurMinInfluence::scip_exit)
{
   cout << "SCIP_DECL_HEUREXIT\n";

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
SCIP_DECL_HEURINITSOL(HeurMinInfluence::scip_initsol)
{
   cout << "SCIP_DECL_HEURINITSOL\n";
   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurMinInfluence::scip_exitsol)
{
   cout << "SCIP_DECL_HEUREXITSOL\n";
   return SCIP_OKAY;
}

/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(HeurMinInfluence::scip_exec)
{
   cout << "SCIP_DECL_HEUREXEC\n";
   assert(heur != NULL);

   SCIP_Bool success = FALSE;

   SCIP_SOL *newsol = SCIPgetBestSol(scip);
   /* SCIP_Bool *visited;
   SCIP_Bool success; */

   assert(result != NULL);
   /* since the timing is SCIP_HEURTIMING_AFTERLPNODE, the current node should have an LP */
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
      return SCIP_OKAY;

   /* copy the current LP solution to the working solution */
   SCIP_CALL(SCIPlinkLPSol(scip, sol));

   /* for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      cout << SCIPgetSolVal(scip, sol, z[a]) << endl;
   } */
   /* for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (SCIPgetSolVal(scip, sol, infSet[v][i].var) > 0)
         {
            cout << SCIPvarGetName(infSet[v][i].var) << " = "
                 << SCIPgetSolVal(scip, sol, infSet[v][i].var) << endl;
         }
      }
   } */
   // TODO verify which case the heuristic can fail, maybe the size of active set
   *result = SCIP_DIDNOTFIND;

   /* allocate local memory */
   SCIP_CALL(SCIPcreateSol(scip, &newsol, heur));

   greedyConstruction(scip, newsol);

   // due to construction we already know, that the solution will be feasible
   SCIP_CALL(SCIPtrySol(scip, newsol, FALSE, FALSE, FALSE, FALSE, FALSE, &success));
   if (success)
   {  
      cout << "solution sussesful!\n";
      *result = SCIP_FOUNDSOL;
   }
   else
      cout << "solution failed!\n";
   
   /* free all local memory */
   SCIP_CALL(SCIPfreeSol(scip, &newsol));

   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
/* SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurMinInfluence::clone)
{
   return new HeurMinInfluence(scip);
} */
