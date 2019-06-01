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
    DNodeConsMap &p_vert_con     /**< array of partitioning constraints */
    ) : ObjPricer(scip, p_name, "Finds inf-set with negative redcost.", 0, TRUE),
        instance(p_instance),
        z(p_arc_var),
        x(p_vert_var),
        arcCons(p_arc_con),
        vertCons(p_vert_con)
{
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
   std::cout << "\n PRICER INIT CALLED DHFLKAHFLAKSDHÇFLHSAÇFLÇALSDFLÇASFÇLSAFÇSDFSF" << std::endl;
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
   std::cout << "\n PRICING CALLED DHFLKAHFLAKSDHÇFLHSAÇFLÇALSDFLÇASFÇLSAFÇSDFSF\n"
             << std::endl;

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
      std::cout << "Dual farkas arc solution: " << endl;
      /* compute the dual solution of the variable associated wiht each arc coverage constraints */
      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a); // just to print
         DNode v = instance.g.target(a);
         std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v] << ": ";

         dualArcValues[a] = SCIPgetDualsolLinear(scip, arcCons[a]);

         std::cout << dualArcValues[a] << std::endl;
      }
      std::cout << "vertex dual soluton: " << endl;
      /* compute the dual solution of the variable associated wiht each vertex coverage constraints */
      for (DNodeIt v(instance.g); v != INVALID; ++v)
      {
         dualVertValues[v] = SCIPgetDualsolLinear(scip, vertCons[v]);
         std::cout << instance.nodeName[v] << ": " << dualVertValues[v] << std::endl;
      }
   }

   // iterate over all vertices to compute the pricing
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      /* compute the minimum cost influencing set w.r.t. dual values */
      list<DNode> infSet;
      SCIP_Real reduced_cost = findMinCostInfluencingSet(v, dualArcValues, dualVertValues[v], infSet);

      /* add influencing set variable */
      if (SCIPisNegative(scip, reduced_cost))
      {
         std::cout << "Negative reduced cost found" << endl;
         return addInfluencingSetVar(scip, v, infSet);
      }
   }
   std::cout << "Negative reduced cost NOT found" << endl;

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

   std::cout << "CALL SCIP FARKAS ÇDALH ONLDKASHFLASDKHFÇLKSAD ...\n"
             << std::endl;

   /* call pricing routine */
   SCIP_CALL(pricing(scip, TRUE));

   return SCIP_OKAY;
}

/** add influencing-set variable to problem */
SCIP_RETCODE ObjPricerGLCIP::addInfluencingSetVar(SCIP *scip, const DNode &v, const list<DNode> &infSet) const
{
   double cost = costInfluencingSet(v, infSet);
   ScipVar *var = new ScipPriceVar(scip, cost);

   SCIP_CALL(SCIPaddCoefLinear(scip, vertCons[v], var->var, 1.0));

   for (DNode u : infSet)
   {
      Arc a = findArc(instance.g, u, v);
      SCIP_CALL(SCIPaddCoefLinear(scip, arcCons[a], var->var, -1.0));
   }

   return SCIP_OKAY;
}

/**
 * Computes the cost paid to activate a vertex v with a given weight of influence
 */
double ObjPricerGLCIP::cheapestIncentive(const DNode &v, double exertedInfluence) const
{
   double cost = 0;
   // assuming that the incentives are sorted in an increasing order
   // uses the first incentive that overcomes the threshold of v
   for (unsigned int i = 0; i < instance.incentives[v].size(); i++)
   {
      if (exertedInfluence + instance.incentives[v][i] >= instance.threshold[v])
      {
         cost = instance.incentives[v][i];
         break;
      }
   }
   return cost;
}

/**
 * Computes the cost paid to activate a vertex v with a given set of incoming neigobors
 */
double ObjPricerGLCIP::costInfluencingSet(const DNode &v, const list<DNode> &nodes) const
{
   int thr = instance.threshold[v];
   double cost = 0;

   // activation function
   double exertedInfluence = 0;
   for (DNode u : nodes)
   {
      Arc a = findArc(instance.g, u, v);
      exertedInfluence += instance.influence[a];
   }

   for (unsigned int i = 0; i < instance.incentives[v].size(); i++)
   {
      if (exertedInfluence + instance.incentives[v][i] >= thr)
      {
         cost = instance.incentives[v][i];
         break;
      }
   }
   return cost;
}

/** return negative reduced cost tour (uses restricted shortest path dynamic programming algorithm) 
 *
 *  Algorithm:
 *    - Create a matrix min_cost[n+1][W+1], where n is the number of incoming neighbors
 *      with distinct weight of influence and W is the threshold of vertex
 *    - Initialize 0-th row with SCIPinfinity and 0-th column with 0
 *    - Fill the matrix
 *          min_cost[i][j] = min(min_cost[i-1][j], cost[i-1] + min_cost[i-1][j-weight[i-1]])
 */
SCIP_Real ObjPricerGLCIP::findMinCostInfluencingSet(
    const DNode &v,                  /**< vertex to be influenced */
    const ArcValueMap &dualArcValue, /**< dual solution of arc constraints */
    const double dualVertValue,      /**< dual solution of vertex constraints */
    const list<DNode> &infSet        /**< list of incoming neighbors */
    ) const
{
   //infSet.clear();
   InDegMap<Digraph> inDegrees(instance.g);
   int n = inDegrees[v];               // number of itens to put into the knapsack
   int W = (int)instance.threshold[v]; // capacity to fill

   SCIP_Real **minCost = new SCIP_Real *[n + 1];
   double *wt = new double[n];
   double *costs = new double[n];

   //initialize weight of influence vector and costs vector
   int i = 0;
   DNode *neighbors = new DNode[n];
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
      minCost[i][0] = cheapestIncentive(v, 0) - dualVertValue;
   }

   // fill 0th row
   for (i = 1; i <= W; i++)
      minCost[0][i] = 9999999;

   // check the weight of influence from a vertex u over v for each u
   // and fill the matrix according to the condition
   for (i = 1; i <= n; i++)
   {
      for (int j = 1; j <= W; j++)
      {
         // get minimum cost by including or excluding the i-th node
         int col = max(j - wt[i - 1], 0.0); // to avoid negative index
         minCost[i][j] = min(minCost[i - 1][j], minCost[i - 1][col] + costs[i - 1] + cheapestIncentive(v, j) - cheapestIncentive(v, col));
      }
   }

   //get the solution
   int j = W;
   for (i = n; i > 0; i--)
   {
      if (minCost[i][j] != minCost[i - 1][j])
      {
         // infSet.push_front(neighbors[i - 1]);
         j -= wt[i];
      }
   }

   SCIP_Real redCost = minCost[n][W];

   delete[] wt;
   delete[] costs;
   delete[] neighbors;
   for (i = 0; i <= n; i++)
      delete[] minCost[i];
   delete[] minCost;
   return redCost;
}