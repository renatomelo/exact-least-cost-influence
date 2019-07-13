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
    //DNodeSCIPVarsMap &p_inf_set
    DNodeInfSetsMap &p_inf_set)
    : ObjPricer(scip, p_name, "Finds influencing set with negative reduced cost.", 0, TRUE),
      instance(p_instance),
      z(p_arc_var),
      x(p_vert_var),
      arcCons(p_arc_con),
      vertCons(p_vert_con),
      infSet(p_inf_set)
{
   // arcMarker isOnSolution(instance.g);
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

   return SCIP_OKAY;
}

/** perform pricing
 *
 *  @todo compute minimum knapsack problem w.r.t. duals
 *  Obs: here we are not considering to perform Farkas pricing
 */
SCIP_RETCODE ObjPricerGLCIP::pricing(SCIP *scip, bool isFarkas) const
{
   /*  std::cout << "\n---------------------------- PRICING CALLED ----------------------------"
             << std::endl; */
   DNodeValueMap dualVertValues(instance.g);
   ArcValueMap dualArcValues(instance.g);

   if (isFarkas)
   {
      std::cout << "Dual farkas arc solution: " << endl;
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
   }
   else
   {
      //std::cout << "Dual arc solution: " << endl;
      /* compute the dual solution of the variable associated wiht each arc coverage constraints */
      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         dualArcValues[a] = SCIPgetDualsolLinear(scip, arcCons[a]);
         /* DNode u = instance.g.source(a); // just to print
         DNode v = instance.g.target(a);
         std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v] << ": ";
         std::cout << dualArcValues[a] << std::endl; */
      }
      //std::cout << "vertex dual soluton: " << endl;
      /* compute the dual solution of the variable associated wiht each vertex coverage constraints */
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         dualVertValues[v] = SCIPgetDualsolLinear(scip, vertCons[v]);
         //std::cout << instance.nodeName[v] << ": " << dualVertValues[v] << std::endl;
      }
   }

   SCIP_RETCODE retcode;
   // iterate over all vertices to compute the pricing
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      /* compute the minimum cost influencing set w.r.t. dual values */
      set<DNode> nodes;
      SCIP_Real reduced_cost = findMinCostInfluencingSet(scip, v, dualArcValues, dualVertValues[v], nodes);

      //std::cout << "Pricing the vertex: " << instance.nodeName[v] << "  ";
      /* add influencing set variable */
      if (SCIPisNegative(scip, reduced_cost))
      {
         //std::cout << "Negative reduced cost: " << reduced_cost << std::endl;
         //return addInfluencingSetVar(scip, v, nodes);
         retcode = addInfluencingSetVar(scip, v, nodes);
         if (retcode != SCIP_OKAY)
            return retcode;
      }
      //std::cout << " size of infSet: " << nodes.size() << std::endl;
   }

   /* std::cout << "\n---------------------------- EXITING PRICING ----------------------------"
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
   SCIPdebugMsg(scip, "call scip_farkas ...\n");

   std::cout << "\n---------------------------- CALL SCIP FARKAS ----------------------------"
             << std::endl;

   /* call pricing routine */
   SCIP_CALL(pricing(scip, TRUE));

   return SCIP_OKAY;
}

/** add influencing-set variable to problem */
SCIP_RETCODE ObjPricerGLCIP::addInfluencingSetVar(SCIP *scip, const DNode &v, const set<DNode> &nodes) const
{
   std::string name;
   if (nodes.size() > 0)
   {
      std::stringstream stream;
      for (DNode u : nodes)
         stream << instance.nodeName[u];
      name = "infSetVar_" + instance.nodeName[v] + "_" + stream.str();
   }
   else
      name = "infSetVar" + instance.nodeName[v] + "_empty";

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
   /* std::cout << "Adding variable for vertex: " << instance.nodeName[v] << " "
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
      }
   }

   // save the variable
   infSet[v].push_back(in);

   SCIP_CALL(SCIPreleaseVar(scip, &var));
   return SCIP_OKAY;
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

   //initialize weight of influence vector and costs vector
   int i = 0;
   vector<DNode> neighbors(n);
   for (InArcIt a(instance.g, v); a != INVALID; ++a)
   {
      neighbors[i] = instance.g.source(a);
      costs[i] = dualArcValue[a];
      wt[i++] = instance.influence[a];
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
         // get minimum cost by including or excluding the i-th node
         int col = max(j - wt[i - 1], 0.0); // to avoid negative index
         minCost[i][j] = min(minCost[i - 1][j],
                             minCost[i - 1][col] + costs[i - 1] +
                                 GLCIPBase::cheapestIncentive(instance, v, j) -
                                 GLCIPBase::cheapestIncentive(instance, v, col));
         /* std::cout << "Cheapest incentive of " << instance.nodeName[v] << " with j = " << j << ": " << cheapestIncentive(v, j)
                   << ", and col = " << col << ": " << cheapestIncentive(v, col) << std::endl; */
         // stop whether the reduced cost is negative
         if (SCIPisNegative(scip, minCost[i][j]))
         {
            /* std::cout << "Table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
            for (int p = 0; p <= n; p++)
            {
               for (int q = 0; q <= W; q++)
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
                  //j -= wt[i - 1];
                  k = max(k - wt[l - 1], 0.0);
               }
            }

            redCost = minCost[i][j];

            delete[] wt;
            delete[] costs;
            //delete[] neighbors;
            for (i = 0; i <= n; i++)
               delete[] minCost[i];
            delete[] minCost;

            return redCost;
         }
      }
      //std::cout << std::endl;
   }

   /* std::cout << "Table of dynamic program for vertex " << instance.nodeName[v] << std::endl;
   for (i = 0; i <= n; i++)
   {
      for (int j = 0; j <= W; j++)
      {
         std::cout << minCost[i][j] << "\t";
      }
      std::cout << std::endl;
   } */

   redCost = minCost[n][W];

   delete[] wt;
   delete[] costs;
   //delete[] neighbors;
   for (i = 0; i <= n; i++)
      delete[] minCost[i];
   delete[] minCost;

   return redCost;
}