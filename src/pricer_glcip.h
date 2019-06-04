#ifndef __SCIP_PRICER_GLCIP_H__
#define __SCIP_PRICER_GLCIP_H__

#include "objscip/objscip.h"
#include "scip/pub_var.h"

#include <float.h>
#include <math.h>
#include <set>
#include <lemon/list_graph.h>
#include <lemon/unionfind.h>
#include <lemon/gomory_hu.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include "mygraphlib.h"
#include <lemon/preflow.h>
#include "easyscip.h"
#include <scip/scip.h>
#include <scip/cons_linear.h>
#include <scip/scipdefplugins.h>
#include <scip/scip_solvingstats.h>
#include <scip/pub_misc.h>
#include "glcipinstance.h"
#include "glcipsolution.h"
#include "cyclecutsgenerator.h"
#include <deque>
#include <vector>
#include <list>

using namespace std;
using namespace scip;

typedef Digraph::ArcMap<SCIP_CONS*> ArcConsMap;
typedef Digraph::NodeMap<SCIP_CONS*> DNodeConsMap;

class ObjPricerGLCIP : public ObjPricer
{
private:
   GLCIPInstance                 instance;   /**< influence matrix */

   ArcSCIPVarMap&                z;           /**< map of arc variables */
   DNodeSCIPVarMap&              x;           /**< map of vertex variables */
   ArcConsMap&                   arcCons;     /**< map of arc constraints */
   DNodeConsMap&                 vertCons;    /**< map of partitioning constraints */

public:
   /** Constructs the pricer object with the data needed */
   ObjPricerGLCIP(  
      SCIP*                               scip,         /**< SCIP pointer */
      const char*                         p_name,       /**< name of pricer */
      GLCIPInstance&                p_instance,   /**< problem data */
      ArcSCIPVarMap&                p_arc_var,    /**< map of arc variables */
      DNodeSCIPVarMap&              p_vert_var,   /**< map of arc variables */
      ArcConsMap&                   p_arc_con,    /**< map of arc constraints */
      DNodeConsMap&                 p_vert_con    /**< map of vertex constraints */
      );

   /** Destructs the pricer object. */
   virtual ~ObjPricerGLCIP();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

    /** farkas pricing method of variable pricer for infeasible LPs */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas); /** farkas pricing method of variable pricer for infeasible LPs */

   /** perform pricing */
   SCIP_RETCODE pricing(SCIP* scip, bool isFarkas) const;

   /** add influencing-set variable to problem */
   SCIP_RETCODE addInfluencingSetVar(SCIP *scip, const DNode &v, const list<DNode> &infSet) const;

   /** return negative reduced cost influencing set (uses minimum knapsack dynamic programming algorithm)*/
  double findMinCostInfluencingSet(
      const DNode&                      v,            /**< vertex to be influenced */
      const ArcValueMap&                dualValues,   /**< map of dual values associated to arc-constraints */
      const double                      dualVertValue,  /**< dual solution of vertex constraints */
      list<DNode>&                infSet        /**< list of influencing neighbors */
      ) const;

  // SCIP_RETCODE incentivesForAll(SCIP* scip) const;
   double cheapestIncentive(const DNode& v, double exertedInfluence) const;
   double costInfluencingSet(const DNode& v, const list<DNode>& nodes) const;
};

#endif
