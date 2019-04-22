#ifndef CUTGENERATOR_H
#define CUTGENERATOR_H
//#define SCIP_DEBUG
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include "objscip/objscip.h"
#include <scip/cons_linear.h>
#include "mygraphlib.h"
#include <list>
#include <lemon/list_graph.h>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <set>
#include <stack>
#include <scip/scip_sol.h>
#include <scip/lp.h>
#include <lemon/dijkstra.h>
#include "glcipinstance.h"

using namespace easyscip;
using namespace scip;

typedef Digraph::NodeMap<vector<SCIP_VAR*>> DNodeSCIPVarsMap;
typedef Digraph::ArcMap<SCIP_VAR*> ArcSCIPVarMap;
typedef lemon::Dijkstra<Digraph, ArcValueMap> SptSolver;

class CycleCutsGenerator: public scip::ObjConshdlr{
    public:
        GLCIPInstance &instance;
        DNodeSCIPVarsMap &x;
        ArcSCIPVarMap &z;

        virtual SCIP_DECL_CONSTRANS(scip_trans);
        virtual SCIP_DECL_CONSSEPALP(scip_sepalp);
        virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);
        virtual SCIP_DECL_CONSENFOLP(scip_enfolp);
        virtual SCIP_DECL_CONSENFOPS(scip_enfops);
        virtual SCIP_DECL_CONSCHECK(scip_check);
        virtual SCIP_DECL_CONSLOCK(scip_lock);

        CycleCutsGenerator(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarsMap &x, ArcSCIPVarMap &z);
        ~CycleCutsGenerator();

        SCIP_RETCODE createCycleCuts(
           SCIP*                 scip,               /**< SCIP data structure */
           SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
           const char*           name,               /**< name of constraint */
           SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
           SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
           SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
           SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
           SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
           SCIP_Bool             local,              /**< is constraint only valid locally? */
           SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
           SCIP_Bool             dynamic,            /**< is constraint dynamic? */
           SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
           );

    private:
        double getXValue(SCIP* scip, SCIP_SOL* sol, DNode v);
        bool isValid(SCIP* scip, SCIP_SOL* sol);
        SCIP_RETCODE addCycleInequality(SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLR* conshdlr, vector<DNode> &cycle);
        SCIP_RETCODE findCycleCuts(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_SOL* sol, SCIP_RESULT* result, bool feasible);
};
#endif // CUTGENERATOR_H
