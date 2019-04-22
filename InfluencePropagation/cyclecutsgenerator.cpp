#include "cyclecutsgenerator.h"

#define SCIP_DEBUG

CycleCutsGenerator::CycleCutsGenerator(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z) :
    instance(instance), x(x), z(z),
    ObjConshdlr(scip, "GLCIP Cycle Cuts", "GLCIP Cycle callback constraints", 10000, 10000, 10000, 1, -1, 1, 0,
        FALSE, FALSE, TRUE, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST) {}

// may be used to free data structures
CycleCutsGenerator::~CycleCutsGenerator(){}

/** creates and captures a GLCIP Cycle constraints */
SCIP_RETCODE CycleCutsGenerator::createCycleCuts(
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
){
    SCIP_CONSHDLR* conshdlr;
    SCIP_CONSDATA* consdata = NULL;

    /* find the subtour constraint handler */
    conshdlr = SCIPfindConshdlr(scip, "GLCIP Cycle Cuts");
    if( conshdlr == NULL ){
      SCIPerrorMessage("GLCIP Cycle constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
    }

    /* create constraint */
    SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

    return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(CycleCutsGenerator::scip_trans)
{
   SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, NULL,
       SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
       SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
       SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
       SCIPconsIsStickingAtNode(sourcecons)));

   return SCIP_OKAY;
}

// separation method of constraint handler for LP solution
SCIP_DECL_CONSSEPALP(CycleCutsGenerator::scip_sepalp)
{
    //cout << "consepalp\n";
    bool feasible = true;
    SCIP_CALL(findCycleCuts(scip, conshdlr, NULL, result, &feasible));
    return SCIP_OKAY;
}

// separation method of constraint handler for arbitrary primal solution
SCIP_DECL_CONSSEPASOL(CycleCutsGenerator::scip_sepasol)
{
    //cout << "consepasol\n";
    bool feasible = true;
    *result = SCIP_DIDNOTFIND;
    SCIP_CALL(findCycleCuts(scip, conshdlr, sol, result, &feasible));
    return SCIP_OKAY;
}


// constraint enforcing method of constraint handler for LP solutions
SCIP_DECL_CONSENFOLP(CycleCutsGenerator::scip_enfolp)
{
    //cout << "consenfolp\n";
    bool check = isValid(scip, NULL);
    if(check)
        *result = SCIP_FEASIBLE;
    else{
        *result = SCIP_INFEASIBLE;
        bool feasible;
        SCIP_CALL(findCycleCuts(scip, conshdlr, NULL, result, &feasible));

        if(*result == SCIP_DIDNOTFIND)
            *result = SCIP_INFEASIBLE;
    }

    return SCIP_OKAY;
}/*lint !e715*/

// constraint enforcing method of constraint handler for pseudo solutions
SCIP_DECL_CONSENFOPS(CycleCutsGenerator::scip_enfops)
{
    //cout << "consenfops\n";
    bool check = isValid(scip, NULL);
    if(check)
        *result = SCIP_FEASIBLE;
    else
        *result = SCIP_SOLVELP;

    return SCIP_OKAY;
} /*lint !e715*/

// feasibility check method of constraint handler for primal solutions
SCIP_DECL_CONSCHECK(CycleCutsGenerator::scip_check)
{
    //cout << "conscheck\n";
    bool check = isValid(scip, sol);
    if(check)
        *result = SCIP_FEASIBLE;
    else
        *result = SCIP_INFEASIBLE;

    return SCIP_OKAY;
} /*lint !e715*/

// variable rounding lock method of constraint handler
/*
 * The CONSLOCK callback provides dual information for a single constraint. It has to tell SCIP, which variables are existing in the given constraint, and in which way modifications of these variables may affect the feasibility of the constraint.

    For each variable that is affected by the constraint, the callback should call SCIPaddVarLocks():

        If the constraint may become violated by decreasing the value of a variable, it should call
            SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the (positive)
            constraint infeasible and rounding up is potentially rendering the negation of the constraint infeasible.
        If the constraint may become violated by increasing the value of a variable, it should call
            SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the constraint's
            negation infeasible and rounding down is potentially rendering the constraint itself infeasible.
        If the constraint may become violated by changing the variable in any direction, it should call
            SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
*/
SCIP_DECL_CONSLOCK(CycleCutsGenerator::scip_lock){
    // round z up may cause cycle inequalities to be infeasible
    for(ArcIt a(instance.g); a != INVALID; ++a){
        SCIP_CALL(SCIPaddVarLocksType(scip, z[a], locktype, nlocksneg, nlockspos));
    }

    // round x down may cause cycle inequalities to be infeasible
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        SCIP_CALL(SCIPaddVarLocksType(scip, x[v], locktype, nlockspos, nlocksneg));
    }

    return SCIP_OKAY;
} /*lint !e715*/

// check if current solution respect cycle constraints
bool CycleCutsGenerator::isValid(SCIP* scip, SCIP_SOL* sol){
    assert(result != NULL);
    *result = SCIP_DIDNOTFIND;

    ArcValueMap weight(instance.g);

    // construct edges weights according to w[(i, j)] = x_i - z_{i,j}
    for(ArcIt a(instance.g); a != INVALID; ++a){
        weight[a] = SCIPgetSolVal(scip, sol, x[instance.g.source(a)]) - SCIPgetSolVal(scip, sol, z[a]);
    }

    // for each node u in the graph, we will run dijkstra to get the shortest path from u to its neighbors
    for(DNodeIt u(instance.g); u != INVALID; ++u){
        // run dijkstra starting from u
        SptSolver spt(instance.g, weight);
        spt.run(u);

        // for each neighbor of u, we check if it violates the cycle inequality
        for(InArcIt a(instance.g, u); a != INVALID; ++a){
            DNode v = instance.g.source(a);

            // in order to violate cycle inequalities we need that L_{u, v} + w_{v, u} < x_u
            if(spt.dist(v) + weight[a] < SCIPgetSolVal(scip, sol, x[u])){
                return false;
            }
        }
    }

    return true;
}

// add capacity cuts
// \sum_{(i, j) \in C} z_{i, j} - \sum_{i \in V(C) \backslash \{k\}} x_i <= 0 \forall k \in V(C), \text{cycle } C \subseteq A,
SCIP_RETCODE CycleCutsGenerator::addCycleInequality(SCIP* scip, SCIP_CONSHDLR* conshdlr, list<int> &cycle){
    SCIP_ROW* row;

    // (the bounds in SCIP are in the following format LHS <= expr <= RHS)
    SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "greater1Cut", -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE));
    SCIP_CALL(SCIPcacheRowExtensions(scip, row));

    // traverse the cycle adding the x and z variables
    for(int i = cycle.size() - 1; i >= 0; i--){
        // pick adjacent nodes in the cycle arc
        DNode u = cycle[i];
        DNode v;

        if(i == 0){
            v = cycle[cycle.size() - 1];
        }
        else{
            v = cycle[i - 1];
        }

        Arc a = findArc(instance.g, u, v);

        // add x and z to the constraint
        if(i > 0){
            SCIPaddVarToRow(scip, row, x[u], -1.0);
        }
        SCIPaddVarToRow(scip, row, z[a], 1.0);
    }

    //Add the cut to the LP
    SCIP_CALL(SCIPflushRowExtensions(scip, row));
    SCIP_Bool infeasible;
    SCIP_CALL(SCIPaddCut(scip, sol, row, TRUE, &infeasible));

    return SCIP_OKAY;
}

// this is the main separation routine
SCIP_RETCODE CycleCutsGenerator::findCycleCuts(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_SOL* sol, SCIP_RESULT* result, bool feasible){
    assert(result != NULL);
    *result = SCIP_DIDNOTFIND;

    ArcValueMap weight(instance.g);

    // construct edges weights according to w[(i, j)] = x_i - z_{i,j}
    for(ArcIt a(instance.g); a != INVALID; ++a){
        weight[a] = SCIPgetSolVal(scip, sol, x[instance.g.source(a)]) - SCIPgetSolVal(scip, sol, z[a]);
    }

    // for each node u in the graph, we will run dijkstra to get the shortest path from u to its neighbors
    for(DNodeIt u(instance.g); u != INVALID; ++u){
        // run dijkstra starting from u
        SptSolver spt(instance.g, weight);
        spt.run(u);

        // for each neighbor of u, we check if it violates the cycle inequality
        for(InArcIt a(instance.g, u); a != INVALID; ++a){
            DNode v = instance.g.source(a);

            // in order to violate cycle inequalities we need that L_{u, v} + w_{v, u} < x_u
            if(spt.dist(v) + weight[a] < SCIPgetSolVal(scip, sol, x[u])){
                // pick the nodes in the violated cycle
                list<DNode> cycle;

                for (DNode curr = v; curr != u; curr = spt.predNode(curr)){
                    if (curr != INVALID && spt.reached(curr)){
                        cycle.push_back(curr);
                        //cout << g.id(curr) << ", ";
                    }
                }

                cycle.push_back(u);

                // now that we have a cycle, add the corresponding inequality
                addCycleInequality(scip, conshdlr, cycle);
            }
        }
    }

    return SCIP_OKAY;
}
