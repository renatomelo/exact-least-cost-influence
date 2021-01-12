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

typedef Digraph::NodeMap<SCIP_VAR *> DNodeSCIPVarMap;
typedef Digraph::ArcMap<SCIP_VAR *> ArcSCIPVarMap;
typedef Digraph::NodeMap<vector<SCIP_VAR *>> DNodeSCIPVarsMap;

typedef Digraph::ArcMap<SCIP_CONS *> ArcConsMap;
typedef Digraph::NodeMap<SCIP_CONS *> DNodeConsMap;

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

    static void addCoverageConstraints(
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

    static void addAllSmallDirectedCycles2(
        SCIP *scip,
        GLCIPInstance &instance,
        DNodeSCIPVarMap &x,
        ArcSCIPVarMap &z);

    static double cheapestIncentive(
        const GLCIPInstance &instance,
        const DNode &v,
        double exertedInfluence);

    static int getIndexOfChepeastIncentive(GLCIPInstance &instance, DNode &node);
    static int getIndexOfChepeastIncentive(GLCIPInstance &instance, DNode &node, double exertedInf);

    static double costInfluencingSet(
        const GLCIPInstance &instance,
        const DNode &v,
        const set<DNode> &nodes);

    static bool intersects(set<DNode> set1, set<DNode> set2);

    static void getSuportGraph(
        SCIP *scip,
        GLCIPInstance &instance,
        SCIP_SOL *sol,
        ArcSCIPVarMap &z,
        Digraph &new_graph);
};

/**
 * Models the integer linear program to solve the special case of GLCIP with additivelly separated activation function
 * and branch and cut with the generalized cycle elimination constraints
 */
class ArcModel : public GLCIPBase
{
public:
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

class ArcModelWithBounds : public GLCIPBase
{
public:
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

// data structure for each generalized set X and a vertex i of GPCs
typedef struct phi
{
    set<DNode> generalizedSet;
    SCIP_ROW *row;
    DNode k;
    SCIP_Real dualVal;
} Phi;

// structure to storage the set of nodes belonging to a influencing set
class InfluencingSet
{
private:
    GLCIPInstance &instance;
    DNode v;
    set<DNode> nodes;
    SCIP_VAR *var;
    double cost;
    string name;

public:
    InfluencingSet(GLCIPInstance &_instance, DNode _v) : instance(_instance), v(_v){};
    InfluencingSet(
        GLCIPInstance &_instance,
        DNode _v,
        set<DNode> _nodes) : instance(_instance), v(_v), nodes(_nodes){};
    ~InfluencingSet(){};

    void addNode(DNode u)
    {
        nodes.insert(u);
    }

    const set<DNode> &getNodes()
    {
        return nodes;
    }

    /**
    * give a significative name to the associated var using the (previously) given nodes
    */
    void giveName()
    {
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
    }

    void setName(string newname)
    {
        name = newname;
    }
    string getName()
    {
        if (name.empty())
        {
            giveName();
        }
        return name;
    }

    void setCost(double c)
    {
        cost = c;
    }

    double getCost()
    {
        return cost;
    }

    void setVar(SCIP_VAR *variable)
    {
        var = variable;
    }

    SCIP_VAR *getVar()
    {
        return var;
    }
};

// structure to storage the set of nodes belonging to a influencing set
/* typedef struct influencing_set
{
    set<DNode> nodes;
    SCIP_VAR *var;
    double cost;
} InfluencingSet; */

// structure to represent each pair (d, phis) of the dynamic program for the pricing
typedef struct pair_of_phis_and_d
{
    int d;
    vector<Phi> phis; // remove later
    vector<int> ids;  //index for the phis
} PairOfPhisAndD;

typedef Digraph::NodeMap<vector<InfluencingSet>> DNodeInfSetsMap;

/**
 * Exact method based on COV model which add all influencing-set variables to the model and performs a branch-and-cut. 
 * Here we use the generalized cycle elimination constraints separated only in the integer case.
 */
class CovModelAllVariables : public GLCIPBase
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

class CovModel : public GLCIPBase
{
public:
    static bool isFeasible(GLCIPInstance &instance, GLCIPSolution &solution);

    static void constructSoltion(SCIP *scip,
                                 GLCIPInstance &instance,
                                 GLCIPSolution &solution,
                                 ArcSCIPVarMap &z,
                                 DNodeInfSetsMap &infSet);

    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

#endif
