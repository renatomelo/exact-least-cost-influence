/**@file   GeneralizedPropagationCons.cpp
 * @brief  Generalized propagation constraint handler for GLCIP problems
 */
#include "generalizedpropagationcons.h"

struct SCIP_ConsData
{
   Digraph *graph;
};

/** 
 * separates generalized propagation constraints
 */
SCIP_RETCODE GeneralizedPropagation::sepaGeneralizedPropCons(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr,  //the constraint handler itself
    SCIP_CONS **conss,        //array of constraint to process
    int nConss,               //number of constraints to process
    int nUsefulConss,         //number of useful (non-obsolete) constraints to process
    SCIP_SOL *sol,            //primal solution that should be separated
    SCIP_RESULT *result,      //pointer to store the result of the separation call
    set<DNode> generalizedSet // set of vertices to be separated
)
{
   *result = SCIP_DIDNOTFIND;

   // add a constraint for each vertex k in the generilized set (k \in X)
   for (DNode k : generalizedSet)
   {
      SCIP_ROW *row;

      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "GPC_separation", 0, SCIPinfinity, FALSE, FALSE, TRUE));
      SCIP_CALL(SCIPcacheRowExtensions(scip, row));
      SCIPaddVarToRow(scip, row, x[k], -1.0);
      for (DNode v : generalizedSet)
      {
         for (int i = 0; i < infSet[v].size(); i++)
         {
            set<DNode> intersection;
            set_intersection(generalizedSet.begin(), generalizedSet.end(),
                             infSet[v][i].nodes.begin(), infSet[v][i].nodes.end(),
                             std::inserter(intersection, intersection.begin()));

            if (intersection.size() == 0)
               SCIPaddVarToRow(scip, row, infSet[v][i].var, 1.0);
         }
      }
      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      SCIP_Bool infeasible;
      SCIP_CALL(SCIPaddCut(scip, sol, row, TRUE, &infeasible));
   }

   if (infeasible)
      *result = SCIP_CUTOFF;
   else
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** 
 * frees specific constraint data
 */
SCIP_DECL_CONSDELETE(GeneralizedPropagation::scip_delete)
{
   assert(consdata != NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(GeneralizedPropagation::scip_trans)
{
   SCIP_CONSDATA *sourcedata;
   SCIP_CONSDATA *targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   targetdata = NULL;

   SCIP_CALL(SCIPallocBlockMemory(scip, &targetdata));
   targetdata->graph = sourcedata->graph;
   capture_graph(targetdata->graph);

   /* create target constraint */
   SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
                            SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
                            SCIPconsIsStickingAtNode(sourcecons)));

   return SCIP_OKAY;
}

/**TODO call here the separation heuristics for the GPC and/or the exact separation
 *  
 * separation method of constraint handler for LP solution
 *
 *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
 *  which means that a valid LP solution exists.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
SCIP_DECL_CONSSEPALP(GeneralizedPropagation::scip_sepalp)
{
   //implement the fischetti's model for separation
   SCIP_CALL(exactSeparation(scip, conshdlr, conss, nconss, nusefulconss, NULL, result));

   return SCIP_OKAY;
}

/** TODO maybe not necessary but verify 
 * separation method of constraint handler for arbitrary primal solution
 *
 *  Separates all constraints of the constraint handler. The method is called outside the LP solution loop (e.g., by
 *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *  Instead, the method should produce cuts that separate the given solution.
 *
 */
SCIP_DECL_CONSSEPASOL(GeneralizedPropagation::scip_sepasol)
{
   // put here some separation method for an primal solution
   SCIP_CALL(heuristicSeparation(scip, conshdlr, conss, nconss, nusefulconss, sol, result));

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was solved.
 *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
 *  cutting plane.
 *
 *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
 *  priorities until the first constraint handler returned with the value SCIP_CUTOFF, SCIP_SEPARATED,
 *  SCIP_REDUCEDDOM, SCIP_CONSADDED, or SCIP_BRANCHED.
 *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
 *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
 *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
 *  SOS-branching on non-integral solutions).
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
SCIP_DECL_CONSENFOLP(GeneralizedPropagation::scip_enfolp)
{
   *result = SCIP_FEASIBLE;

   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CONSDATA *consdata;
      Digraph *graph;
      SCIP_Bool found;
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);

      //implemente some funtion here in this sense
      found = findTour(scip, graph, NULL);

      // if a tour was found, we generate a cut constraint saying that there must be at least two outgoing edges
      if (found)
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was not solved.
 *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or adding an additional constraint.
 *  Separation is not possible, since the LP is not processed at the current node. All LP informations like
 *  LP solution, slack values, or reduced costs are invalid and must not be accessed.
 *
 *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
 *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
 *  the value SCIP_CUTOFF, SCIP_REDUCEDDOM, SCIP_CONSADDED, SCIP_BRANCHED, or SCIP_SOLVELP.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
 *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
 *  its constraints and return any other possible result code.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the SCIP_LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 */
SCIP_DECL_CONSENFOPS(GeneralizedPropagation::scip_enfops)
{
   *result = SCIP_FEASIBLE;

   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CONSDATA *consdata;
      Digraph *graph;
      SCIP_Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);

      // if a tour is found, the solution must be infeasible
      found = findTour(scip, graph, NULL);
      if (found)
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions
 *
 *  The given solution has to be checked for feasibility.
 *  
 *  The check methods of the active constraint handlers are called in decreasing order of their check
 *  priorities until the first constraint handler returned with the result SCIP_INFEASIBLE.
 *  The integrality constraint handler has a check priority of zero. A constraint handler which can
 *  (or wants) to check its constraints only for integral solutions should have a negative check priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to check feasibility even on non-integral solutions must have a
 *  check priority greater than zero (e.g. if the check is much faster than testing all variables for
 *  integrality).
 *
 *  In some cases, integrality conditions or rows of the current LP don't have to be checked, because their
 *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
 *  'checklprows' is FALSE.
 *
 *  possible return values for *result:
 *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
SCIP_DECL_CONSCHECK(GeneralizedPropagation::scip_check)
{
   *result = SCIP_FEASIBLE;

   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CONSDATA *consdata;
      Digraph *graph;
      SCIP_Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);

      // if a subtour is found, the solution must be infeasible
      found = findTour(scip, graph, sol);
      if (found)
      {
         *result = SCIP_INFEASIBLE;
         if (printreason)
         {
            SCIP_CALL(SCIPprintCons(scip, conss[i], NULL));
            SCIPinfoMessage(scip, NULL, "violation: graph has a subtour\n");
         }
      }
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
SCIP_DECL_CONSPROP(GeneralizedPropagation::scip_prop)
{
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler
 *
 *  This method is called, after a constraint is added or removed from the transformed problem.
 *  It should update the rounding locks of all associated variables with calls to SCIPaddVarLocksType(),
 *  depending on the way, the variable is involved in the constraint:
 *  - If the constraint may get violated by decreasing the value of a variable, it should call
 *    SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg), saying that rounding down is
 *    potentially rendering the (positive) constraint infeasible and rounding up is potentially rendering the
 *    negation of the constraint infeasible.
 *  - If the constraint may get violated by increasing the value of a variable, it should call
 *    SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlocksneg, nlockspos), saying that rounding up is
 *    potentially rendering the constraint's negation infeasible and rounding up is potentially rendering the
 *    constraint itself infeasible.
 *  - If the constraint may get violated by changing the variable in any direction, it should call
 *    SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg).
 */
SCIP_DECL_CONSLOCK(GeneralizedPropagation::scip_lock)
{
   for (DNodeIt v(instance.g); v != INVALID; ++v)
      SCIP_CALL(SCIPaddVarLocksType(scip, x[v], locktype, nlockspos, nlocksneg));

   return SCIP_OKAY;
}

/** variable deletion method of constraint handler
 *
 *  This method should iterate over all constraints of the constraint handler and delete all variables
 *  that were marked for deletion by SCIPdelVar().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
SCIP_DECL_CONSDELVARS(GeneralizedPropagation::scip_delvars)
{
   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable *GeneralizedPropagation::clone)
{
   *valid = true;
   return new GeneralizedPropagation(scip);
}

/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method, which copies a constraint from one SCIP data structure into a other
 *  SCIP data structure.
 */
SCIP_DECL_CONSCOPY(GeneralizedPropagation::scip_copy)
{
   SCIP_CONSHDLR *conshdlr;
   SCIP_CONSDATA *consdata;

   /* find the subtour constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "GPC");
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("CPC constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   SCIP_CALL(SCIPallocBlockMemory(scip, &consdata));

   Digraph *graph = instance.g;
   consdata->graph = graph;
   capture_graph(consdata->graph);

   /* create constraint */
   SCIP_CALL(SCIPcreateCons(scip, cons, (name == NULL) ? SCIPconsGetName(sourcecons) : name,
                            conshdlr, consdata, initial, separate, enforce, check,
                            propagate, local, modifiable, dynamic, removable, FALSE));

   *valid = true;
   return SCIP_OKAY;
}

SCIP_RETCODE createGenPropagationCons(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name)
{

   //std::cout << "Entenring in createConsArcMarker()\n";
   SCIP_CONSHDLR *conshdlr;
   SCIP_CONSDATA *consdata;

   // find the arc marker constraint handler
   conshdlr = SCIPfindConshdlr(scip, "GPC");
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("GPC constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   // create constraint data
   //SCIP_CALL(createConsData(scip, &consdata, arc, arcVar, type, node));

   //create constraint
   SCIP_CALL(
       SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, TRUE, TRUE, TRUE,
                      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE));

   /* std::cout << "created constraint: ";
    printConsData(consdata);
 */
   return SCIP_OKAY;
}