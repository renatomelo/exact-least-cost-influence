#include "heur_mininfluence.h"
#include "pricer_glcip.h"
#include <queue>

/** add influencing-set variable to problem */
SCIP_RETCODE addInfluencingSetVar(
    SCIP *scip,
    GLCIPInstance &instance,
    DNode &v,
    set<DNode> &nodes,
    DNodeInfSetsMap &infSet,
    ArcConsMap *arcCons,    /**< map of arc constraints */
    DNodeConsMap *vertCons, /**< map of partitioning constraints */
    vector<Phi> *gpcRows)
{
   // data structure to save the variables and associated nodes
   InfluencingSet ifs(instance, v, nodes);
   ifs.setCost(GLCIPBase::costInfluencingSet(instance, v, nodes));

   string name;
   if (nodes.size() > 0)
   {
      std::stringstream stream;
      const char *separator = "";
      for (DNode u : nodes)
      {
         stream << separator << instance.nodeName[u];
         separator = ",";
      }

      name = "PrLambda_" + instance.nodeName[v] + "_{" + stream.str() + "}";
   }
   else
      name = "PrLambda_" + instance.nodeName[v] + "_empty";
   ifs.setName(name);
   //double cost = GLCIPBase::costInfluencingSet(instance, v, nodes);
   SCIP_VAR *var;
   SCIP_CALL(SCIPcreateVar(scip, &var,
                           ifs.getName().c_str(),   // var name
                           0.0,                     // lower bound
                           SCIPinfinity(scip),      // upper bound
                           ifs.getCost(),           // coeficient in the objective function
                           SCIP_VARTYPE_CONTINUOUS, // continuous variable
                           FALSE,                   // initial variable
                           FALSE,                   // removable variable
                           NULL, NULL, NULL, NULL, NULL));
   // add new variable to the list of variables to price into LP
   SCIP_CALL(SCIPaddVar(scip, var));

   SCIP_CALL(SCIPaddCoefLinear(scip, (*vertCons)[v], var, 1.0));

   ifs.setVar(var);

   if (nodes.size() != 0)
   {
      for (DNode u : nodes)
      {
         Arc a = findArc(instance.g, u, v);
         assert(a != INVALID);
         SCIP_CALL(SCIPaddCoefLinear(scip, (*arcCons)[a], var, 1.0));
         //in.nodes.insert(u);
      }
   }

   //update the GPC rows by adding the new var on each constraint in which the set X contains v
   for (unsigned int i = 0; i < (*gpcRows).size(); i++)
   {
      if ((*gpcRows)[i].generalizedSet.count(v))
      {
         if (!GLCIPBase::intersects((*gpcRows)[i].generalizedSet, ifs.getNodes()))
         {
            SCIPaddVarToRow(scip, (*gpcRows)[i].row, ifs.getVar(), 1.0);
         }
      }
   }

   // save the variable
   infSet[v].push_back(ifs);

   //std::cout << "adding var: " << SCIPvarGetName(var) << std::endl;

   SCIP_CALL(SCIPreleaseVar(scip, &var));
   return SCIP_OKAY;
}

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
SCIP_RETCODE HeurMinInfluence::greedyConstruction(
    SCIP *scip,
    SCIP_SOL *newsol)
{
   set<DNode> actives; // start with a empty solution

   //find the cheapest vertex with Lambda_v_empty in the LP solution
   double minC = 1e+20;
   DNode node = INVALID;
   int idx = -1;
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      for (size_t i = 0; i < infSet[v].size(); i++)
      {
         if (infSet[v][i].getNodes().empty() && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, infSet[v][i].getVar())))
         {
            //choose the minimum cost
            if (infSet[v][i].getCost() < minC)
            {
               minC = infSet[v][i].getCost();
               node = v;
               idx = i;
            }
            //since we fond the variable Lambda_v_empty we don't need to visit the rest
            break; //break the inner loop
         }
      }
   }

   if (idx != -1)
   {
      //fix the variable Lambda_v_empty for the chosen node to start the propagation
      SCIP_CALL(SCIPsetSolVal(scip, newsol, infSet[node][idx].getVar(), 1.0));
      //cout << "\nsetting var: " << SCIPvarGetName(infSet[node][idx].var) << endl;
      SCIP_CALL(SCIPsetSolVal(scip, newsol, x[node], 1.0));
      actives.insert(node);
   }

   while (actives.size() < instance.alpha * instance.n)
   {
      DNode v = INVALID;
      set<DNode> nodes;
      //double cost;
      double minCost = getMinIncentiveNode(scip, actives, v);

      assert(v != INVALID);
      SCIP_CALL(SCIPsetSolVal(scip, newsol, x[v], 1.0));

      double sum = 0;
      // save v's influencing-set and activate it
      for (InArcIt a(instance.g, v); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         if (actives.count(u) && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, z[a])))
         {
            sum += instance.influence[a];

            nodes.insert(u);
            SCIP_CALL(SCIPsetSolVal(scip, newsol, z[a], 1.0));

            // if the minCost was achieved we don't need to increase the influencing-set
            if (minCost >= GLCIPBase::costInfluencingSet(instance, v, nodes))
               break;
         }
      }

      double slack = (sum + minCost) - instance.threshold[v];

      //if influencing-set 'nodes' isn't minimal remove surplus elements
      for (auto u = nodes.begin(); u != nodes.end();)
      {
         Arc a = findArc(instance.g, *u, v);
         assert(a != INVALID);

         if (instance.influence[a] < slack)
         {
            //set isn't minimal, removing vertex
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
         if (nodes == infSet[v][i].getNodes())
         {
            SCIP_CALL(SCIPsetSolVal(scip, newsol, infSet[v][i].getVar(), 1.0));
            //cout << "setting var: " << SCIPvarGetName(infSet[v][i].var) << endl;
            varFound = TRUE;
            break;
         }
      }

      //isn't the variable in the model?
      if (!varFound)
      {
         //var not found";
         //create a new var and add it to the model
         addInfluencingSetVar(scip, instance, v, nodes, infSet, arcCons, vertCons, gpcRows);

         //the new variable was added in the back of the list infSet[v], then its position is infSet[v].size()-1
         int position = infSet[v].size() - 1;
         SCIP_CALL(SCIPsetSolVal(scip, newsol, infSet[v][position].getVar(), 1.0));
      }

      actives.insert(v);
   }

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
   //cout << "SCIP_DECL_HEURINIT\n";

   /* create heuristic data */
   SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurMinInfluence::scip_exit)
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
SCIP_DECL_HEURINITSOL(HeurMinInfluence::scip_initsol)
{
   //cout << "SCIP_DECL_HEURINITSOL\n";
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
   //sortVerticesTest(instance);
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
      /* cout << "NEW SOL\n";
      double totalCost = 0;
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         for (unsigned int i = 0; i < infSet[v].size(); i++)
         {
            if (SCIPgetSolVal(scip, newsol, infSet[v][i].var) > 0)
            {
               cout << SCIPvarGetName(infSet[v][i].var) << " = "
                    << SCIPgetSolVal(scip, newsol, infSet[v][i].var) << endl;
               totalCost += infSet[v][i].cost;
            }
         }
      }

      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         if (SCIPgetSolVal(scip, newsol, z[a]) > 0)
            cout << SCIPvarGetName(z[a]) << " = " << SCIPgetSolVal(scip, newsol, z[a]) << endl;
      } */

      // get the cost of this solution
      //cout << "solution total cost: " << totalCost << endl;
   }

   /* free all local memory */
   SCIP_CALL(SCIPfreeSol(scip, &newsol));
   //cout << "free all local memory \n";

   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
/* SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurMinInfluence::clone)
{
   return new HeurMinInfluence(scip);
} */
