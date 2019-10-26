#ifndef __GLCIPBASE_H__
#define __GLCIPBASE_H__

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


using namespace easyscip;
using namespace scip;

typedef Digraph::NodeMap<SCIP_VAR*> DNodeSCIPVarMap;
typedef Digraph::ArcMap<SCIP_VAR*> ArcSCIPVarMap;
typedef Digraph::NodeMap<vector<SCIP_VAR*>> DNodeSCIPVarsMap;

/**
 * Class base containing the declaration methods that are common for all the implemented models
 * for the problems and algorithms related to the GLCIP problem 
 */
class GLCIPBase
{
public:
    
    static void addLinkingConstraints(
        SCIP *scip, 
        GLCIPInstance &instance, 
        DNodeSCIPVarMap &x, 
        ArcSCIPVarMap &z);

    static void addCoverageConstraints (
        SCIP *scip, 
        GLCIPInstance &instance, 
        DNodeSCIPVarMap &x);

    static void addCycleConstraints(
        SCIP *scip,
        GLCIPInstance &instance, 
        DNodeSCIPVarMap &x, 
        ArcSCIPVarMap &z, 
        DNodeIntMap &predMap, 
        Arc &backArc);

    static void dfsSmallCycles(
        SCIP *scip, 
        GLCIPInstance &instance, 
        DNodeSCIPVarMap &x, 
        ArcSCIPVarMap &z, 
        DNodeIntMap &colors,
        DNodeIntMap &predMap, 
        DNode curr, 
        int level, 
        int rootId);
                                
    static void addSmallCycleConstraints(
        SCIP *scip, GLCIPInstance &instance,
        DNodeSCIPVarMap &x,
        ArcSCIPVarMap &z);
    /**
     * method created to substitute the methods above
     */ 
    static void addAllSmallDirectedCycles(
        SCIP *scip, 
        GLCIPInstance &instance, 
        DNodeSCIPVarMap &x, 
        ArcSCIPVarMap &z);

    static double cheapestIncentive(
        const GLCIPInstance &instance,
        const DNode &v,
        double exertedInfluence);
    
    static double costInfluencingSet(
        const GLCIPInstance &instance,
        const DNode &v, 
        const set<DNode> &nodes);

    static bool intersects(set<DNode> &set1, set<DNode> &set2);
//     static void addCuttingPlanes(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z);
};

/**
 * Models the integer linear program to solve the special case of GLCIP with additivelly separated activation function
 * and branch and cut with the generalized cycle elimination constraints
 */
class ArcModel: public GLCIPBase
{
public:
 //   static SCIP_RETCODE addCuttingPlanes(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z);
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

// structure to storage the set of nodes belonging to a influencing set
typedef struct influencing_set{
    set<DNode> nodes;
    SCIP_VAR* var;
    double cost;
}InfluencingSet;

// data structure for each generalized set X and a vertex i of GPCs
typedef struct phi{
    set<DNode> generalizedSet;
    SCIP_ROW* row;
    DNode k;
    SCIP_Real dualVal;
}Phi;

typedef Digraph::NodeMap< vector<InfluencingSet> > DNodeInfSetsMap;

/**
 * Exact method based on COV model which add all influencing-set variables to the model and performs a branch-and-cut. 
 * Here we use the generalized cycle elimination constraints separated only in the integer case.
 */
class CovModelAllVariables: public GLCIPBase
{
public:
    CovModelAllVariables();
    ~CovModelAllVariables();
    static vector<InfluencingSet> powerSet(
        GLCIPInstance &instance,
        vector<DNode> neighbors, 
        DNode &v);
    
    static double costInfluencingSet(GLCIPInstance &instance, DNode v, set<DNode> nodes);
    
    static void addPropagationConstraints(SCIP *scip,
                                          GLCIPInstance &instance, 
                                          DNodeSCIPVarMap &x, 
                                          DNodeInfSetsMap &infSet);
    
    static void addChosenArcsConstraints(SCIP *scip, 
                                         GLCIPInstance &instance, 
                                         ArcSCIPVarMap &z, 
                                         DNodeInfSetsMap &infSet);
    
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

class CovModel: public GLCIPBase
{
public:
    static bool isFeasible(GLCIPInstance &instance, GLCIPSolution &solution);
   
    static void constructSoltion(SCIP *scip, 
                                  GLCIPInstance &instance, 
                                  GLCIPSolution &solution, 
                                  ArcSCIPVarMap& z, 
                                  DNodeInfSetsMap& infSet);
    
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

#endif