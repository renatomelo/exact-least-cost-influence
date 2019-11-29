#include "pricer_glcip.h"

/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricerGLCIP::ObjPricerGLCIP(
    SCIP *scip,                  /**< SCIP pointer */
    const char *p_name,          /**< name of pricer */
    GLCIPInstance &p_instance,   /**< problem data */
    ArcSCIPVarMap &p_arc_var,    /**< matrix of arc variables */
    DNodeSCIPVarMap &p_vert_var, /**< matrix of arc variables */
    ArcConsMap &p_arc_con,       /**< matrix of arc constraints */
    DNodeConsMap &p_vert_con,    /**< array of partitioning constraints */
    vector<Phi> &p_gpcrows,
    DNodeInfSetsMap &p_inf_set)
    : ObjPricer(
          scip, p_name, "Finds influencing set with negative reduced cost.", 0, TRUE),
      instance(p_instance),
      z(p_arc_var),
      x(p_vert_var),
      arcCons(p_arc_con),
      vertCons(p_vert_con),
      gpcrows(p_gpcrows),
      infSet(p_inf_set)
{
   /* for (ArcIt a(instance.g); a != INVALID; ++a)
      this->isOnSolution[a] = 0; */
}

ObjPricerGLCIP::~ObjPricerGLCIP() {}

/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(ObjPricerGLCIP::scip_init)
{
   //std::cout << "\n---------------------------- PRICER INIT CALLED ----------------------------" << std::endl;

   for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      SCIP_CALL(SCIPgetTransformedVar(scip, z[a], &z[a]));
      SCIP_CALL(SCIPgetTransformedCons(scip, arcCons[a], &arcCons[a]));
   }
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      SCIP_CALL(SCIPgetTransformedVar(scip, x[v], &x[v]));
      SCIP_CALL(SCIPgetTransformedCons(scip, vertCons[v], &vertCons[v]));
   }

   // verify if is necessary to get the tranformed row generated in GPC

   return SCIP_OKAY;
}

/** perform pricing
 *
 *  @todo compute minimum knapsack problem w.r.t. duals
 *  Obs: here we are not considering to perform Farkas pricing
 */
SCIP_RETCODE ObjPricerGLCIP::pricing(SCIP *scip, bool isFarkas) const
{
   /* std::cout << "\n---------------------------- PRICING CALLED ----------------------------"
             << std::endl; */
   DNodeValueMap dualVertValues(instance.g);
   ArcValueMap dualArcValues(instance.g);

   if (isFarkas)
   {
      std::cout << "Dual farkas arc solution at node: " << SCIPgetFocusDepth(scip) << endl;
      /* compute the dual farkas of the variable associated wiht each arc coverage constraints */
      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         DNode v = instance.g.target(a);
         std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v] << ": ";

         dualArcValues[a] = SCIPgetDualfarkasLinear(scip, arcCons[a]);

         std::cout << dualArcValues[a] << std::endl;
      }
      std::cout << "vertex dual farkas soluton: " << endl;
      /* compute the dual farkas of the variable associated wiht each vertex coverage constraints */
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         dualVertValues[v] = SCIPgetDualfarkasLinear(scip, vertCons[v]);
         std::cout << instance.nodeName[v] << ": " << dualVertValues[v] << std::endl;
      }

      for (unsigned int i = 0; i < gpcrows.size(); i++)
      {
         gpcrows[i].dualVal = SCIProwGetDualfarkas(gpcrows[i].row);
      }
   }
   else
   {
      //std::cout << "Dual arc solution: " << endl;
      //get the dual solution of the variable associated wiht each arc coverage constraints
      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         dualArcValues[a] = SCIPgetDualsolLinear(scip, arcCons[a]);
         /* DNode u = instance.g.source(a); // just to print
         DNode v = instance.g.target(a);
         std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v] << ": ";
         std::cout << dualArcValues[a] << std::endl; */
      }
      //std::cout << "vertex dual soluton: " << endl;
      //get the dual solution of the variable associated wiht each vertex coverage constraints
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         dualVertValues[v] = SCIPgetDualsolLinear(scip, vertCons[v]);
         //std::cout << instance.nodeName[v] << ": " << dualVertValues[v] << std::endl;
      }

      for (unsigned int i = 0; i < gpcrows.size(); i++)
      {
         gpcrows[i].dualVal = SCIProwGetDualsol(gpcrows[i].row);
         /* cout << "SCIProwGetDualsol(gpcrows[i].row) = "
              << SCIProwGetDualsol(gpcrows[i].row) << endl; */
      }
   }

   SCIP_RETCODE retcode;
   // iterate over all vertices to compute the pricing
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      /* compute the minimum cost influencing set w.r.t. dual values */
      set<DNode> nodes;
      SCIP_Real reduced_cost = findMinCostInfluencingSet4(scip, v, dualArcValues, dualVertValues[v], nodes);

      /* add influencing set variable */
      if (SCIPisNegative(scip, reduced_cost))
      {
         //std::cout << "negative reduced cost: " << reduced_cost << std::endl;
         //return addInfluencingSetVar(scip, v, nodes);
         retcode = addInfluencingSetVar(scip, v, nodes);
         if (retcode != SCIP_OKAY)
            return retcode;
      }
      //std::cout << " size of infSet: " << nodes.size() << std::endl;
   }

  /*  std::cout << "\n---------------------------- EXITING ----------------------------"
             << std::endl; */

   //SCIP_CALL(SCIPwriteTransProblem(scip, "glcip_transformed.lp", "lp", FALSE));
   return SCIP_OKAY;
}

/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible tour with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP
 *                      solution is optimal
 */
SCIP_DECL_PRICERREDCOST(ObjPricerGLCIP::scip_redcost)
{
   /* std::cout << "\n--------------------- CALL SCIP_PRICER_REDUCED_COST ----------------------"
             << std::endl; */
   /* set result pointer, see above */
   *result = SCIP_SUCCESS;

   /* if (SCIPgetDepth(scip) != 0)
   {
      return SCIP_OKAY;
   } */

   //cout << "SCIPgetDepth(scip) = " << SCIPgetDepth(scip) << endl;
   /* call pricing routine */
   SCIP_CALL(pricing(scip, FALSE));

   return SCIP_OKAY;
}

/** Pricing of additional variables if LP is infeasible.
 *
 *  - get the values of the dual Farks multipliers you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible tour with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 */
SCIP_DECL_PRICERFARKAS(ObjPricerGLCIP::scip_farkas)
{
   std::cout << "---------------------------- CALL SCIP FARKAS ----------------------------\n";

   //SCIP_CALL(SCIPwriteTransProblem(scip, "glcip_transformed.lp", "lp", FALSE));
   /* call pricing routine */
   //SCIP_CALL(pricing(scip, TRUE));

   //SCIP_CALL(SCIPwriteTransProblem(scip, "glcip_transformed2.lp", "lp", FALSE));
   //exit(0);

   return SCIP_OKAY;
}

/** add influencing-set variable to problem */
SCIP_RETCODE ObjPricerGLCIP::addInfluencingSetVar(SCIP *scip, const DNode &v, const set<DNode> &nodes) const
{
   std::string name;
   if (nodes.size() > 0)
   {
      std::stringstream stream;
      const char *separator = "";
      for (DNode u : nodes)
      {
         stream << separator << instance.nodeName[u];
         separator = ",";
      }

      name = "Lambda_" + instance.nodeName[v] + "_{" + stream.str() + "}";
   }
   else
      name = "Lambda_" + instance.nodeName[v] + "_empty";

   double cost = GLCIPBase::costInfluencingSet(instance, v, nodes);
   SCIP_VAR *var;
   SCIP_CALL(SCIPcreateVar(scip, &var,
                           name.c_str(),            // var name
                           0.0,                     // lower bound
                           SCIPinfinity(scip),      // upper bound
                           cost,                    // coeficient in the objective function
                           SCIP_VARTYPE_CONTINUOUS, // continuous variable
                           FALSE,                   // initial variable
                           FALSE,                   // removable variable
                           NULL, NULL, NULL, NULL, NULL));
   // add new variable to the list of variables to price into LP
   SCIP_CALL(SCIPaddPricedVar(scip, var, 1.0));

   SCIP_CALL(SCIPaddCoefLinear(scip, vertCons[v], var, 1.0));

   // data structure to save the variables and associated nodes
   InfluencingSet in;
   in.var = var;
   in.cost = cost;

   /* std::cout << "adding variable for vertex: " << instance.nodeName[v] << " "
             << name << ", cost = " << cost << std::endl; */
   if (nodes.size() != 0)
   {
      for (DNode u : nodes)
      {
         //std::cout << instance.nodeName[u] << " ";
         Arc a = findArc(instance.g, u, v);
         assert(a != INVALID);
         SCIP_CALL(SCIPaddCoefLinear(scip, arcCons[a], var, -1.0));
         in.nodes.insert(u);
         //isAble[a] = TRUE;
      }
   }

   //update the GPC rows by adding the new var on each constraint in which the set X contains v
   for (unsigned int i = 0; i < gpcrows.size(); i++)
   {
      if (gpcrows[i].generalizedSet.count(v))
      {
         //cout << "updating the " << i + 1 << "-th GPC row: ";
         //SCIPprintRow(scip, gpcrows[i].row, NULL);

         if (!GLCIPBase::intersects(gpcrows[i].generalizedSet, in.nodes))
         {
            SCIPaddVarToRow(scip, gpcrows[i].row, in.var, 1.0);
            //cout << " adding " << SCIPvarGetName(in.var) << " to the " << i + 1 << "-th GPC row\n";
         }
         /*  cout << "updated " << i + 1 << "-th GPC row\n";
         SCIPprintRow(scip, gpcrows[i].row, NULL); */
      }
   }

   // save the variable
   infSet[v].push_back(in);

   std::cout << "adding var: " << SCIPvarGetName(var) << std::endl;

   SCIP_CALL(SCIPreleaseVar(scip, &var));
   return SCIP_OKAY;
}

SCIP_Real sumOfPhis(DNode v, vector<Phi> gpcrows)
{
   double sum = 0;

   for (unsigned int i = 0; i < gpcrows.size(); i++)
   {
      if (gpcrows[i].generalizedSet.count(v))
      {
         sum += gpcrows[i].dualVal;
      }
   }

   return sum;
}

SCIP_Real sumOfPhisFromIds(vector<Phi> &gpcrows, vector<int> ids)
{
    double sum = 0;

    for (int i : ids)
    {
        sum += gpcrows[i].dualVal;
    }

    return sum;
}

void getValidPhis(DNode v, vector<Phi> &gpcrows, vector<Phi> &new_gpcrows)
{
   //cout << "gpcrows.size() = " << gpcrows.size() << endl;
   for (unsigned int i = 0; i < gpcrows.size(); i++)
   {
      if (gpcrows[i].generalizedSet.count(v))
      {
         new_gpcrows.push_back(gpcrows[i]);
      }
   }
}

void getValidPhisIds(DNode v, vector<Phi> &gpcrows, vector<int> &new_gpcrows)
{
   //cout << "gpcrows.size() = " << gpcrows.size() << endl;
   for (unsigned int i = 0; i < gpcrows.size(); i++)
   {
      if (gpcrows[i].generalizedSet.count(v))
      {
         new_gpcrows.push_back(i);
      }
   }
}

/** return negative reduced cost influencing set (uses minimum knapsack dynamic 
 *  programming algorithm)
 *  Algorithm:
 *    - Create a matrix min_cost[n+1][W+1], where n is the number of incoming neighbors
 *      with distinct weight of influence and W is the threshold of vertex
 *    - Initialize 0-th row with SCIPinfinity and 0-th column with 0
 *    - Fill the matrix
 *          min_cost[i][j] = min(min_cost[i-1][j], cost[i-1] + min_cost[i-1][j-weight[i-1]])
 */
SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet(
    SCIP *scip,
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &dualArcValue, /**< dual solution of arc constraints */
    const double dualVertValue,      /**< dual solution of vertex constraints */
    set<DNode> &nodes                /**< list of incoming neighbors */
    ) const
{
   SCIP_Real redCost;
   nodes.clear();
   InDegMap<Digraph> inDegrees(instance.g);
   if (inDegrees[v] == 0)
      return (GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue);

   int n = inDegrees[v];               // number of itens to put into the knapsack
   int W = (int)instance.threshold[v]; // capacity to fill

   SCIP_Real **minCost = new SCIP_Real *[n + 1];
   double *wt = new double[n];
   double *costs = new double[n];

   /* std::cout << "looking the influencing-set of node " << instance.nodeName[v]
             << " at BnB node " << SCIPgetFocusDepth(scip) << "\n"; */

   //initialize weight of influence vector and costs vector
   //double extInf = 0; // external influence
   int i = 0;
   vector<DNode> neighbors(n);
   for (InArcIt a(instance.g, v); a != INVALID; ++a)
   {
      //if (isOnSolution[a] == 0)
      {
         //std::cout << "arc " << SCIPvarGetName(z[a]) << " explored at node "
         //          << SCIPgetFocusDepth(scip) << "\n";
         neighbors[i] = instance.g.source(a);
         costs[i] = dualArcValue[a];
         wt[i++] = instance.influence[a];
      }
      /* else if (isOnSolution[a] == -1)
      {
         //std::cout << "arc " << SCIPvarGetName(z[a]) << " marked to NOT be in solution at node "
         //          << SCIPgetFocusDepth(scip) << "\n";
         n--; //reduces the size of the knapsack
         if (n == 0)
            return (GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue);
      }
      else
      {
         //std::cout << "arc marked to be in solution: ";
         // add the source of the arc to the influencing set of v
         nodes.insert(instance.g.source(a));
         n--;
         W = W - instance.influence[a]; // decreases the knapsack capacity
         //std::cout << "knapsack's capacity changed to " << W << "\n";
         //assert(W >= 0);

         extInf += instance.influence[a];

         if (W <= 0)
            return (GLCIPBase::cheapestIncentive(instance, v, extInf) - dualVertValue + dualArcValue[a]);
      } */
   }

   // fill 0th column
   // initial reduced cost no influence from neighbors
   for (i = 0; i <= n; i++)
   {
      minCost[i] = new SCIP_Real[W + 1];
      minCost[i][0] = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue;
   }

   //std::cout << "\nIncentive cost with no influence from neighbors: "
   //<< cheapestIncentive(v, 0) << std::endl;
   // fill 0th row
   for (i = 1; i <= W; i++)
      minCost[0][i] = 1e+9;

   // check the weight of influence from a vertex u over v for each u
   // and fill the matrix according to the condition
   for (i = 1; i <= n; i++)
   {
      for (int j = 0; j <= W; j++)
      {
         //cout << "i = " << i << " j = " << j << ", W = " << W << endl;
         // get minimum cost by including or excluding the i-th node
         int col = max(j - wt[i - 1], 0.0); // to avoid negative index
         minCost[i][j] = min(minCost[i - 1][j],
                             minCost[i - 1][col] + costs[i - 1] +
                                 GLCIPBase::cheapestIncentive(instance, v, j) -
                                 GLCIPBase::cheapestIncentive(instance, v, col));

         // stop whether the reduced cost is negative
         if (SCIPisNegative(scip, minCost[i][j]))
         {
            std::cout << "table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
            for (int p = 0; p <= n; p++)
            {
               for (int q = 0; q <= W; q++)
               {
                  std::cout << minCost[p][q] << "\t";
               }
               std::cout << std::endl;
            }

            //get the solution
            int k = j;
            for (int l = i; l > 0; l--)
            {
               if (minCost[l][k] != minCost[l - 1][k])
               {
                  nodes.insert(neighbors[l - 1]);
                  k = max(k - wt[l - 1], 0.0);
               }
            }

            redCost = minCost[i][j];

            delete[] wt;
            delete[] costs;
            for (i = 0; i <= n; i++)
               delete[] minCost[i];
            delete[] minCost;

            return redCost;
         }
      }
      //std::cout << std::endl;
   }

   std::cout << "Table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
   for (i = 0; i <= n; i++)
   {
      for (int j = 0; j <= W; j++)
      {
         std::cout << minCost[i][j] << "\t";
      }
      std::cout << std::endl;
   }

   redCost = minCost[n][W];

   delete[] wt;
   delete[] costs;
   //delete[] neighbors;
   for (i = 0; i <= n; i++)
      delete[] minCost[i];
   delete[] minCost;

   return redCost;
}

//new dynamic program for solving the pricing subproblem
SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet2(
    SCIP *scip,
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &dualArcValue, /**< dual solution of arc constraints */
    const double dualVertValue,      /**< dual solution of vertex constraints */
    set<DNode> &nodes                /**< list of incoming neighbors */
    ) const
{
   //cout << "pricing vertex: " << instance.nodeName[v] << endl;
   SCIP_Real redCost;
   nodes.clear();
   InDegMap<Digraph> inDegrees(instance.g);

   int n = inDegrees[v]; // number of itens to put into the knapsack
   //int W = (int)instance.threshold[v]; // capacity to fill

   SCIP_Real **minCost = new SCIP_Real *[n + 1];
   double *wt = new double[n];
   double *costs = new double[n];

   //initialize weight of influence vector and costs vector
   int i = 0;
   int totalInfluence = 0;
   vector<DNode> neighbors(n);
   for (InArcIt a(instance.g, v); a != INVALID; ++a)
   {
      totalInfluence += instance.influence[a];
      neighbors[i] = instance.g.source(a);
      costs[i] = dualArcValue[a];
      wt[i++] = instance.influence[a];
   }

   //the number of columns in the matrix will be the total influence from neighbors
   int W = totalInfluence;

   // fill 0th column
   // initial reduced cost no influence from neighbors
   for (i = 0; i <= n; i++)
   {
      minCost[i] = new SCIP_Real[W + 1];
      minCost[i][0] = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue - sumOfPhis(v, gpcrows);
      //minCost[i][0] = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue;
   }

   // fill 0th row
   for (i = 1; i <= W; i++)
      minCost[0][i] = SCIPinfinity(scip);

   //initialization (no influence from neighbors)
   vector<vector<PairOfPhisAndD>> listOfPairs(n + 1);
   PairOfPhisAndD firstPair;
   firstPair.d = 0;
   getValidPhis(v, gpcrows, firstPair.phis);
   listOfPairs[0].push_back(firstPair);

   //cout << "listOfPairs[0][0].phis.size() = " << listOfPairs[0][0].phis.size() << endl;

   /* cout << "rows of gpc in which X contains node " << instance.nodeName[v] << ": \n";
   for (unsigned int l = 0; l < listOfPairs[0][0].phis.size(); l++)
   {
      //SCIPprintRow(scip, L0[j].row, NULL);
      for (DNode w : listOfPairs[0][0].phis[l].generalizedSet)
         cout << instance.nodeName[w] << " ";
      cout << endl;
   } */

   // check the weight of influence from a vertex u over v for each u
   // and fill the matrix according to the condition
   for (i = 1; i <= n; i++)
   {
      //for (unsigned int j = 0; j < listOfPairs[i-1].size(); j++)
      for (int j = 0; j <= W; j++)
      {
         //cout << "i = " << i << ", j = " << j << ", n = " << n << ", W = " << W << endl;
         if (j >= (int)listOfPairs[i - 1].size())
         {
            minCost[i][j] = SCIPinfinity(scip);
            continue;
         }
         PairOfPhisAndD currentPair = listOfPairs[i - 1][j];

         //L_j = L_j \cup {(d, PHI)}
         listOfPairs[i].push_back(currentPair);

         // phi' = {(X,k) \in Phi : i is not in X}
         PairOfPhisAndD prime;
         double sum = 0;
         for (unsigned int l = 0; l < currentPair.phis.size(); l++)
         {
            // remember that neighbors[i-1] is the i-th neigbor
            if (!currentPair.phis[l].generalizedSet.count(neighbors[i - 1]))
               prime.phis.push_back(currentPair.phis[l]);
            else
               //sum all the phis excluding the elements in prime
               sum += currentPair.phis[l].dualVal;
         }

         //cout << "size of prime phis = " << prime.phis.size() << endl;
         //cout << "sum of phis excluding the elements in prime: " << sum << endl;

         // get minimum cost by including or excluding the i-th node
         //int col = max(j - wt[i - 1], 0.0); // to avoid negative index
         int col = j + wt[i - 1];
         minCost[i][j] = minCost[i - 1][j]; //no influence from i
         //minCost[i][j] = min(minCost[i - 1][j], minCost[i][j]);

         minCost[i][col] = min(minCost[i][col],
                               minCost[i - 1][j] + costs[i - 1] + sum +
                                   GLCIPBase::cheapestIncentive(instance, v, col) -
                                   GLCIPBase::cheapestIncentive(instance, v, j));
         prime.d = col;
         listOfPairs[i].push_back(prime);

         // stop whether the reduced cost is negative
         //if (SCIPisNegative(scip, minCost[i][j]))
         if (SCIPisNegative(scip, minCost[i][col]))
         {
            j = col;
            std::cout << "table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
            for (int p = 0; p <= i; p++)
            {
               for (int q = 0; q <= j; q++)
               {
                  std::cout << minCost[p][q] << "\t";
               }
               std::cout << std::endl;
            }

            //get the solution
            int k = j;
            for (int l = i; l > 0; l--)
            {
               if (minCost[l][k] != minCost[l - 1][k])
               {
                  nodes.insert(neighbors[l - 1]);
                  k = max(k - wt[l - 1], 0.0);
               }
            }

            /*  if (sumOfPhis(v, gpcrows) > 0)
            {
               for (DNode y : nodes)
                  cout << instance.nodeName[y] << " ";
               cout << endl;
               exit(0);
            } */

            redCost = minCost[i][j];

            delete[] wt;
            delete[] costs;
            for (i = 0; i <= n; i++)
               delete[] minCost[i];
            delete[] minCost;

            return redCost;
         }
      }
   }

   std::cout << "Table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
   for (i = 0; i <= n; i++)
   {
      for (int j = 0; j <= W; j++)
      {
         std::cout << minCost[i][j] << "\t";
      }
      std::cout << std::endl;
   }
   //cout << "didn't find negative reduced cost\n";
   redCost = minCost[n][W];

   delete[] wt;
   delete[] costs;
   for (i = 0; i <= n; i++)
      delete[] minCost[i];
   delete[] minCost;

   return redCost;
}

SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet3(
    SCIP *scip,
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &dualArcValue, /**< dual solution of arc constraints */
    const double dualVertValue,      /**< dual solution of vertex constraints */
    set<DNode> &nodes                /**< list of incoming neighbors */
    ) const
{
   //cout << "pricing vertex: " << instance.nodeName[v] << endl;
   SCIP_Real redCost;
   nodes.clear();
   int n = 0; // number of itens to put into the knapsack
   for(InArcIt a(instance.g, v); a != INVALID; ++a) n++;
   int W = (int)instance.threshold[v]; // capacity to fill

   SCIP_Real **minCost = new SCIP_Real *[n + 1];
   double *wt = new double[n];
   double *costs = new double[n];

   //initialize weight of influence vector and costs vector
   int i = 0;
   int totalInfluence = 0;
   vector<DNode> neighbors(n);
   for (InArcIt a(instance.g, v); a != INVALID; ++a)
   {
      totalInfluence += instance.influence[a];
      neighbors[i] = instance.g.source(a);
      costs[i] = dualArcValue[a];
      wt[i++] = instance.influence[a];
   }

   // fill 0th column
   // initial reduced cost no influence from neighbors
   for (i = 0; i <= n; i++)
   {
      minCost[i] = new SCIP_Real[W + 1];
      //minCost[i][0] = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue - sumOfPhis(v, gpcrows);
      minCost[i][0] = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue;
   }

   // fill 0th row
   for (i = 1; i <= W; i++)
      minCost[0][i] = SCIPinfinity(scip);

   //initialization (no influence from neighbors)
   vector<vector<PairOfPhisAndD>> listOfPairs(n + 1);
   PairOfPhisAndD firstPair;
   firstPair.d = 0;
   getValidPhis(v, gpcrows, firstPair.phis);
   listOfPairs[0].push_back(firstPair);

   // check the weight of influence from a vertex u over v for each u
   // and fill the matrix according to the condition
   for (i = 1; i <= n; i++)
   {
      //for (unsigned int j = 0; j < listOfPairs[i-1].size(); j++)
      for (int j = 0; j <= W; j++)
      {
         //cout << "i = " << i << ", j = " << j << ", n = " << n << ", W = " << W << endl;
         if (j >= (int) listOfPairs[i - 1].size())
         {
            minCost[i][j] = SCIPinfinity(scip);
            continue;
         }
         PairOfPhisAndD currentPair = listOfPairs[i - 1][j];

         //L_j = L_j \cup {(d, PHI)}
         listOfPairs[i].push_back(currentPair);

         // phi' = {(X,k) \in Phi : i is not in X}
         PairOfPhisAndD prime;
         double sum = 0;
         for (unsigned int l = 0; l < currentPair.phis.size(); l++)
         {
            // remember that neighbors[i-1] is the i-th neigbor
            if (!currentPair.phis[l].generalizedSet.count(neighbors[i - 1]))
               prime.phis.push_back(currentPair.phis[l]);
            else
               //sum all the phis excluding the elements in prime
               sum += currentPair.phis[l].dualVal;
         }

         //cout << "size of prime phis = " << prime.phis.size() << endl;
         //cout << "sum of phis excluding the elements in prime: " << sum << endl;

         // get minimum cost by including or excluding the i-th node
         int col = max(j - wt[i - 1], 0.0); // to avoid negative index
         minCost[i][j] = min(minCost[i - 1][j],
                             minCost[i - 1][col] + costs[i - 1] + sum +
                                 GLCIPBase::cheapestIncentive(instance, v, j) -
                                 GLCIPBase::cheapestIncentive(instance, v, col));

         prime.d = j;
         listOfPairs[i].push_back(prime);

         // stop whether the reduced cost is negative
         if (SCIPisNegative(scip, minCost[i][j]))
         {
            /* std::cout << "table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
            for (int p = 0; p <= i; p++)
            {
               for (int q = 0; q <= j; q++)
               {
                  std::cout << minCost[p][q] << "\t";
               }
               std::cout << std::endl;
            } */

            //get the solution
            int k = j;
            for (int l = i; l > 0; l--)
            {
               if (minCost[l][k] != minCost[l - 1][k])
               {
                  nodes.insert(neighbors[l - 1]);
                  k = max(k - wt[l - 1], 0.0);
               }
            }

            redCost = minCost[i][j];

            delete[] wt;
            delete[] costs;
            for (i = 0; i <= n; i++)
               delete[] minCost[i];
            delete[] minCost;

            return redCost;
         }
      }
   }

   redCost = minCost[n][W];

   delete[] wt;
   delete[] costs;
   for (i = 0; i <= n; i++)
      delete[] minCost[i];
   delete[] minCost;

   return redCost;
}

vector<set<DNode>> ObjPricerGLCIP::powerSet(
    vector<DNode> neighbors,
    DNode v) const
{
   unsigned int pSize = pow(2, neighbors.size());

   vector<set<DNode>> pSet;

   // run from 000...0 to 111...1
   for (unsigned int i = 0; i < pSize; i++)
   {
      set<DNode> subset;

      for (unsigned int j = 0; j < neighbors.size(); j++)
      {
         if (i & (1 << j))
         {
            subset.insert(neighbors[j]);
         }
      }

      double minCost = GLCIPBase::costInfluencingSet(instance, v, subset);

      int sum = 0;
      for (DNode u : subset)
      {
         Arc a = findArc(instance.g, u, v);
         assert(a != INVALID);

         sum += instance.influence[a];
      }
      double slack = (sum + minCost) - instance.threshold[v];
      //cout << "slack = " << slack << endl;

      bool isMinimal = TRUE;
      for (DNode u : subset)
      {
         Arc a = findArc(instance.g, u, v);
         assert(a != INVALID);

         if (instance.influence[a] < slack)
            isMinimal = FALSE;
      }

      //add to the list only if is minimal
      if (isMinimal)
         pSet.push_back(subset);
   }

   return pSet;
}

/**
 * brute force algorithm for the pricing subproblem
 */
SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet4(
    SCIP *scip,
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &pi, /**< dual solution of arc constraints */
    const double mi,      /**< dual solution of vertex constraints */
    set<DNode> &nodes                /**< list of incoming neighbors */
    ) const
{
   //cout << "(brute force) pricing vertex: " << instance.nodeName[v] << endl;
   SCIP_Real redCost = 0;
   nodes.clear();

   //generate all minimal influencing-set for v
   vector<DNode> neighbors;
   for (InArcIt a(instance.g, v); a != INVALID; ++a)
   {
      neighbors.push_back(instance.g.source(a));
   }
   vector<set<DNode>> subsets = powerSet(neighbors, v);

   vector<Phi> listOfPhis;
   getValidPhis(v, gpcrows, listOfPhis);
   
   for (size_t i = 0; i < subsets.size(); i++)
   {
      double sum = 0, sumOfArcDuals = 0, sumOfGPCduals = 0;
      for (DNode u : subsets[i])
      {
         //cout << instance.nodeName[u] << " ";
         
         Arc a = findArc(instance.g, u, v);
         assert(a != INVALID);
         
         //total influence exerted
         sum += instance.influence[a];
         sumOfArcDuals += pi[a];
      }
      //cout << endl;

      //iterate over all GPCs already in the model to get the duals
      for (size_t j = 0; j < listOfPhis.size(); j++)
      {
         //consider only the GPCs in which v belongs to X
         if(!GLCIPBase::intersects(listOfPhis[j].generalizedSet, subsets[i]))
         {
            sumOfGPCduals += listOfPhis[j].dualVal;
         }
      }
      
      redCost = GLCIPBase::cheapestIncentive(instance, v, sum) - mi + sumOfArcDuals - sumOfGPCduals;

      if (SCIPisNegative(scip, redCost))
      {
         for (DNode u : subsets[i])
         {
            nodes.insert(u);
         }
         return redCost;
      }
   }

   //if no reduced cost was found return 0
   return redCost;
}

SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet5(
    SCIP *scip,
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &dualArcValue, /**< dual solution of arc constraints */
    const double dualVertValue,      /**< dual solution of vertex constraints */
    set<DNode> &nodes                /**< list of incoming neighbors */
    ) const
{
    nodes.clear();

    int n = 0; // number of itens to put into the knapsack
    for(InArcIt a(instance.g, v); a != INVALID; ++a) n++;
    int W = (int)instance.threshold[v]; // capacity to fill

    //initialize weight of influence vector and costs vector
    int i = 0;
    double *wt = new double[n];
    double *costs = new double[n];
    vector<DNode> neighbors(n);
    for (InArcIt a(instance.g, v); a != INVALID; ++a)
    {
        neighbors[i] = instance.g.source(a);
        costs[i] = dualArcValue[a];
        wt[i++] = instance.influence[a];
    }
    unordered_map<pair<double, vector<int>>, double> L1;
    unordered_map<pair<double, vector<int>>, double> L2;

    // initialize L1
    vector<int> allPhisIds;
    getValidPhisIds(v, gpcrows, allPhisIds);
    double initialValue = GLCIPBase::cheapestIncentive(instance, v, 0) - dualVertValue
           - sumOfPhisFromIds(gpcrows, allPhisIds);

    // found negative reduced cost
    if(SCIPisNegative(scip, initialValue)){
        cout << "found reduced cost" << endl;
    }

    // create first pair and put into L1
    pair<double, vector<int>> p;
    p.first = 0.0;
    p.second = allPhisIds;
    L1[p] = initialValue;

    for(int j = 0; j < n; j++)
    {
        for(const auto& p: L1){
            double curr, aux;
            if(L2.count(p)){
                aux = L2[p];
            }
            else{
                aux = SCIPinfinity(scip);
            }

            L2[p] = min(aux, L1[p]);

            // construct difference set
            vector<int> phiPrimeIds;
            double diffSum = 0.0;
            for(auto i: p.second){
                if(gpcrows[i].generalizedSet.count(neighbors[j])){
                    diffSum += gpcrows[i].dualVal;
                }
                else{
                    phiPrimeIds.push_back(i);
                }
            }

            pair<double, vector<int>> p2;
            p2.first = p.first + wt[j];
            p2.second = phiPrimeIds;

            if(L2.count(p2)){
                aux = L2[p2];
            }
            else{
                aux = SCIPinfinity(scip);
            }

            L2[p2] = min(aux, L1[p] + diffSum + cost[j] + GLCIPBase::cheapestIncentive(instance, v, p2.first)
                         - GLCIPBase::cheapestIncentive(instance, v, p.first));

            if(SCIPisNegative(scip, L2[p2])){
                cout << "found negative reduced cost" << endl;
            }
        }
        L1 = L2;
        L1.clear();
    }

    delete[] wt;
    delete[] costs;

    return redCost;
}
