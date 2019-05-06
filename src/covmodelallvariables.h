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
#include <vector>

using namespace easyscip;
using namespace scip;

struct influencingSet{
    set<DNode> elements;
};

//typedef Digraph::ArcMap<SCIP_VAR*> ArcSCIPVarMap;
//typedef Digraph::NodeMap<SCIP_VAR*> DNodeSCIPVarMap;
//typedef Digraph::NodeMap< set<SCIP_VAR*> > DNodeSCIPVarsMap;
typedef Digraph::NodeMap< vector<influencingSet> > DNodeInfSetsMap;

class CovModelAllVariables
{
public:
    CovModelAllVariables();
    ~CovModelAllVariables();
    static vector<influencingSet> powerSet(vector<DNode> neighbors);
    static double costInfluencingSet(GLCIPInstance &instance, DNode v, set<DNode> elements);
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};
