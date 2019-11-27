#ifndef __SCIP_PRICER_GLCIP_H__
#define __SCIP_PRICER_GLCIP_H__

#include "GLCIPBase.h"

typedef Digraph::ArcMap<SCIP_CONS *> ArcConsMap;
typedef Digraph::NodeMap<SCIP_CONS *> DNodeConsMap;

class ObjPricerGLCIP : public ObjPricer
{
private:
    GLCIPInstance instance; /**< problem data*/

    ArcSCIPVarMap &z;       /**< map of arc variables */
    DNodeSCIPVarMap &x;     /**< map of vertex variables */
    ArcConsMap &arcCons;    /**< map of arc constraints */
    DNodeConsMap &vertCons; /**< map of partitioning constraints */
    vector<Phi> &gpcrows;
    DNodeInfSetsMap &infSet;
    //ArcBoolMap &isAble;

    //ArcBoolMap arcMarker;   // map to signals the decision about the arcs in the branching rule

public:
    // indicates whether a arc variable is fixed to be in soluton (+1),
    // to not be in solution (-1) or free (0). To be used by the branching rule
    //ArcIntMap &isOnSolution;

    /** Constructs the pricer object with the data needed */
    ObjPricerGLCIP(
        SCIP *scip,                  /**< SCIP pointer */
        const char *p_name,          /**< name of pricer */
        GLCIPInstance &p_instance,   /**< problem data */
        ArcSCIPVarMap &p_arc_var,    /**< map of arc variables */
        DNodeSCIPVarMap &p_vert_var, /**< map of arc variables */
        ArcConsMap &p_arc_con,       /**< map of arc constraints */
        DNodeConsMap &p_vert_con,    /**< map of vertex constraints */
        vector<Phi> &p_gpcrows,
        DNodeInfSetsMap &p_inf_set);
    //ArcBoolMap&      p_isAble,
    //ArcIntMap&      p_isOnSolution);

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
    SCIP_RETCODE addInfluencingSetVar(SCIP *scip, const DNode &v, const set<DNode> &nodes) const;

    /** return negative reduced cost influencing set (uses minimum knapsack dynamic programming algorithm)*/
    double findMinCostInfluencingSet(
        SCIP *scip,
        const DNode &v,                /**< vertex to be influenced */
        const ArcValueMap &dualValues, /**< map of dual values associated to arc-constraints */
        const double dualVertValue,    /**< dual solution of vertex constraints */
        set<DNode> &nodes              /**< list of influencing neighbors */
        ) const;

    double findMinCostInfluencingSet2(
        SCIP *scip,
        const DNode &v,                /**< vertex to be influenced */
        const ArcValueMap &dualValues, /**< map of dual values associated to arc-constraints */
        const double dualVertValue,    /**< dual solution of vertex constraints */
        set<DNode> &nodes              /**< list of influencing neighbors */
        ) const;

    double findMinCostInfluencingSet3(
        SCIP *scip,
        const DNode &v,                /**< vertex to be influenced */
        const ArcValueMap &dualValues, /**< map of dual values associated to arc-constraints */
        const double dualVertValue,    /**< dual solution of vertex constraints */
        set<DNode> &nodes              /**< list of influencing neighbors */
        ) const;

    double findMinCostInfluencingSet4(
        SCIP *scip,
        const DNode &v,                /**< vertex to be influenced */
        const ArcValueMap &pi, /**< map of dual values associated to arc-constraints */
        const double mi,    /**< dual solution of vertex constraints */
        set<DNode> &nodes              /**< list of influencing neighbors */
        ) const;

    vector< set<DNode> > powerSet(
    vector<DNode> neighbors,
    DNode v) const;

    /* double cheapestIncentive(const DNode &v, double exertedInfluence) const;
   double costInfluencingSet(const DNode &v, const list<DNode> &nodes) const; */
};

#endif
