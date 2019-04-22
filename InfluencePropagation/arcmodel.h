#ifndef ARCMODEL_H
#define ARCMODEL_H

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

using namespace easyscip;
using namespace scip;

typedef Digraph::ArcMap<SCIP_VAR*> ArcSCIPVarMap;
typedef Digraph::NodeMap<vector<SCIP_VAR*>> DNodeSCIPVarsMap;

class ArcModel
{
    public:
        static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

#endif // ARCMODEL_H
