#ifndef __SCIP_PRICER_GLCIP_H__
#define __SCIP_PRICER_GLCIP_H__

#include "GLCIPBase.h"

typedef Digraph::ArcMap<SCIP_CONS *> ArcConsMap;
typedef Digraph::NodeMap<SCIP_CONS *> DNodeConsMap;

class ObjPricerGLCIP : public ObjPricer
{
private:
   GLCIPInstance instance; /**< influence matrix */

   ArcSCIPVarMap &z;       /**< map of arc variables */
   DNodeSCIPVarMap &x;     /**< map of vertex variables */
   ArcConsMap &arcCons;    /**< map of arc constraints */
   DNodeConsMap &vertCons; /**< map of partitioning constraints */
   //DNodeSCIPVarsMap&             infSetVar;
   DNodeInfSetsMap &infSet;

public:
   /** Constructs the pricer object with the data needed */
   ObjPricerGLCIP(
       SCIP *scip,                  /**< SCIP pointer */
       const char *p_name,          /**< name of pricer */
       GLCIPInstance &p_instance,   /**< problem data */
       ArcSCIPVarMap &p_arc_var,    /**< map of arc variables */
       DNodeSCIPVarMap &p_vert_var, /**< map of arc variables */
       ArcConsMap &p_arc_con,       /**< map of arc constraints */
       DNodeConsMap &p_vert_con,    /**< map of vertex constraints */
       //DNodeSCIPVarsMap&             p_inf_set_var
       DNodeInfSetsMap &p_inf_set);

   /** Destructs the pricer object. */
   virtual ~ObjPricerGLCIP();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** farkas pricing method of variable pricer for infeasible LPs */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas); /** farkas pricing method of variable pricer for infeasible LPs */

   /** perform pricing */
   SCIP_RETCODE pricing(SCIP *scip, bool isFarkas) const;

   /** add influencing-set variable to problem */
   SCIP_RETCODE addInfluencingSetVar(SCIP *scip, const DNode &v, const list<DNode> &nodes) const;

   /** return negative reduced cost influencing set (uses minimum knapsack dynamic programming algorithm)*/
   double findMinCostInfluencingSet(
       SCIP*                  scip,
       const DNode&           v,                /**< vertex to be influenced */
       const ArcValueMap&     dualValues,       /**< map of dual values associated to arc-constraints */
       const double           dualVertValue,    /**< dual solution of vertex constraints */
       list<DNode>&           nodes             /**< list of influencing neighbors */
       ) const;

   double cheapestIncentive(const DNode &v, double exertedInfluence) const;
   double costInfluencingSet(const DNode &v, const list<DNode> &nodes) const;
};

#endif
