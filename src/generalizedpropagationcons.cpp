/**@file   GeneralizedPropagationCons.cpp
 * @brief  Generalized propagation constraint handler for GLCIP problems
 */
#include "generalizedpropagationcons.h"

typedef lemon::Dijkstra<Digraph, ArcValueMap> SptSolver;

struct SCIP_ConsData
{
   Digraph *graph;
};

//TODO create the suport graph and test if it is a DAG
//peforms a DFS to verify if there is cycle in the graph suport of the feasible solution
SCIP_Bool findTour(
    SCIP *scip,
    SCIP_SOL *sol,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z)
{
   DNodeIntMap color(instance.g); //0=white, 1=gray, 2=black

   //initialization
   for (DNodeIt v(instance.g); v != INVALID; ++v)
      color[v] = 0;

   vector<DNode> stack;
   for (DNodeIt u(instance.g); u != INVALID; ++u)
   {
      double value = SCIPgetVarSol(scip, x[u]);
      //if (SCIPisIntegral(scip, value) && SCIPisEQ(scip, value, 1))
      if ((color[u] == 0) && SCIPisEQ(scip, value, 1))
      {
         stack.push_back(u);
         color[u] = 1;

         while (stack.size() > 0)
         {
            DNode v = stack.back();
            stack.pop_back();

            for (OutArcIt a(instance.g, v); a != INVALID; ++a)
            {
               if (!SCIPisEQ(scip, SCIPgetVarSol(scip, z[a]), 1))
                  continue;

               DNode w = instance.g.target(a);
               if (SCIPisEQ(scip, SCIPgetVarSol(scip, x[w]), 1) && (color[w] == 0))
               {
                  stack.push_back(w);
                  color[w] = 1;
               }
               else if (color[w] == 1)
               {
                  cout << "tour founded!" << endl;
                  return true;
               }
            }

            color[v] = 2;
         }
      }
   }

   return false;
}

set<DNode> greedSetExtensionHeur(
    SCIP *scip,
    SCIP_SOL *sol,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    DNodeInfSetsMap &infSet)
{
   set<DNode> generalizedSet;

   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      double infSetValue = 0;
      //get the empty influensin-set var
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (infSet[v][i].nodes.empty())
         {
            infSetValue = SCIPgetSolVal(scip, sol, infSet[v][i].var);
            break;
         }
      }

      // x[v] > 0
      if (SCIPgetSolVal(scip, sol, x[v]) > SCIPepsilon(scip) && (infSetValue < SCIPgetSolVal(scip, sol, x[v])))
      {
      }
   }

   return generalizedSet;
}

/**
 * Create a copy of a given graph
 */
void digraph_copy(Digraph &old_grpah, Digraph &new_graph)
{
   Digraph::NodeMap<DNode> nodeReference(old_grpah);
   Digraph::ArcMap<Arc> arcReference(new_graph);
   digraphCopy(old_grpah, new_graph).nodeRef(nodeReference).arcCrossRef(arcReference).run();
}
/** 
 * separates generalized propagation constraints
 */
/* SCIP_RETCODE GeneralizedPropagation::sepaGeneralizedPropCons(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr,  //the constraint handler itself
    SCIP_CONS **conss,        //array of constraint to process
    int nConss,               //number of constraints to process
    int nUsefulConss,         //number of useful (non-obsolete) constraints to process
    SCIP_SOL *sol,            //primal solution that should be separated
    SCIP_RESULT *result,      //pointer to store the result of the separation call
    set<DNode> generalizedSet // set of vertices to be separated
) */
SCIP_RETCODE GeneralizedPropagation::addGeneralizedPropCons(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr, //the constraint handler itself
    SCIP_SOL *sol,
    SCIP_RESULT *result,
    set<DNode> generalizedSet)
{
   cout << "ADDING CONSTRAINT" << endl;

   //choose the maximum x[k]
   double max_k = -1;
   int idOfNodek = -1;
   for (DNode v : generalizedSet)
   {
      //cout << "value of x_v: " << SCIPgetSolVal(scip, sol, x[v]) << endl;
      if (SCIPgetSolVal(scip, sol, x[v]) > max_k)
      {
         idOfNodek = instance.g.id(v);
         max_k = SCIPgetSolVal(scip, sol, x[v]);
         if (SCIPisEQ(scip, max_k, 1))
         {
            break;
         }
      }
   }
   assert(idOfNodek != -1);
   DNode k = instance.g.nodeFromId(idOfNodek);

   // add a constraint for each vertex k in the generilized set (k \in X)
   SCIP_ROW *row;

   SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "GPC_separation", 0, SCIPinfinity(scip), FALSE, FALSE, TRUE));
   SCIP_CALL(SCIPcacheRowExtensions(scip, row));

   SCIPaddVarToRow(scip, row, x[k], -1.0);
   for (DNode v : generalizedSet)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
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

   // add cut
   if (SCIPisCutEfficacious(scip, sol, row))
   {
      SCIP_Bool infeasible;
      SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));
      if (infeasible)
         *result = SCIP_CUTOFF;
   }
   SCIP_CALL(SCIPreleaseRow(scip, &row));

   return SCIP_OKAY;
}

SCIP_RETCODE GeneralizedPropagation::sepaGeneralizedPropCons(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr, //the constraint handler itself
                             //    SCIP_CONS **conss,        //array of constraint to process
                             //    int nConss,               //number of constraints to process
                             //    int nUsefulConss,         //number of useful (non-obsolete) constraints to process
    SCIP_SOL *sol,           //primal solution that should be separated
    SCIP_RESULT *result      //pointer to store the result of the separation call
                             //    set<DNode> generalizedSet // set of vertices to be separated
)
{
   cout << "HEURISITC SEPARATION()" << endl;

   assert(result != NULL);
   *result = SCIP_DIDNOTFIND;

   //create a copy of the graph in order to avoid adding repeated cycles in the GPC
   /* Digraph g;
   digraph_copy(instance.g, g); */

   ArcValueMap weight(instance.g);

   // construct edges weights according to w[(i, j)] = x_i - z_{i,j}
   for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      weight[a] = SCIPgetSolVal(scip, sol, x[instance.g.source(a)]) - SCIPgetSolVal(scip, sol, z[a]);
   }

   for (DNodeIt u(instance.g); u != INVALID; ++u)
   {
      //if (SCIPisEQ(scip, SCIPgetSolVal(scip, sol, x[u]), 1.0))
      //if (SCIPgetSolVal(scip, sol, x[u]) > 0.5)
      {
         //cout << "ENTERED IN FIRST IF CONDITION" << endl;
         SptSolver spt(instance.g, weight);
         spt.run(u);

         for (InArcIt a(instance.g, u); a != INVALID; ++a)
         {
            DNode v = instance.g.source(a);

            //TODO add the constraint only the cycle that violates the GPC
            // in order to violate cycle inequalities we need that L_{u, v} + w_{v, u} < x_u
            if (spt.reached(v) && (spt.dist(v) + weight[a] < SCIPgetSolVal(scip, sol, x[u])))
            //if (spt.reached(v) && (spt.dist(v) + weight[a] < SCIPepsilon(scip)))
            {
               *result = SCIP_SEPARATED;

               cout << "Tour found: ";
               set<DNode> tour;
               for (DNode w = v; w != u; w = spt.predNode(w))
               {
                  if (w != INVALID && spt.reached(w))
                  {
                     //tour.insert(instance.g.nodeFromId(g.id(w)));
                     tour.insert(w);
                     cout << instance.nodeName[w] << ", ";
                  }
               }
               tour.insert(u);
               //tour.insert(instance.g.nodeFromId(g.id(u)));
               cout << instance.nodeName[u] << endl;

               addGeneralizedPropCons(scip, conshdlr, sol, result, tour);

               if (*result == SCIP_CUTOFF)
                  return SCIP_OKAY;
            }
         }

         //remove the vertex to avoid add repeated cycles
         //g.erase(u);
      }
   }

   return SCIP_OKAY;
}
/**
 * Integer linear program to find the maximally violated inequality 
 */
SCIP_RETCODE GeneralizedPropagation::exactSeparation(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr, //the constraint handler itself
    SCIP_SOL *sol,           //primal solution that should be separated
    SCIP_RESULT *result      //pointer to store the result of the separation call
)
{
   SCIP *new_scip = NULL;

   // initialize SCIP enviroment
   SCIP_CALL(SCIPcreate(&new_scip));

   SCIP_CALL(SCIPincludeDefaultPlugins(new_scip));
   SCIP_CALL(SCIPsetSeparating(new_scip, SCIP_PARAMSETTING_OFF, TRUE));

   // create empty problem
   SCIP_CALL(SCIPcreateProb(new_scip, "GPC_separation", 0, 0, 0, 0, 0, 0, 0));
   SCIP_CALL(SCIPsetIntParam(new_scip, "display/verblevel", 3));

   // declaring variables
   DNodeSCIPVarMap belongsToX(instance.g);  //indicates membership of nodes to set X
   DNodeSCIPVarMap isOnRHS(instance.g);     //equal to one node i is on the right-hand side of GPC
   DNodeInfSetsMap validInfSet(instance.g); //equal to one iff node i \in X and U has no
                                            // intersection with X
   //decides whether the lifting to one on the right-hand side is applied or not
   ScipVar *var_lifting = new ScipBinVar(scip, "lifting_rhs", -1);
   SCIP_VAR *liftingRHS = var_lifting->var;

   //create "belong to X" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      ScipVar *var = new ScipBinVar(new_scip, "v_" + instance.nodeName[v], 0.0);
      belongsToX[v] = var->var;

      //fix the variable in zero in order to consider only variables which x[v] > 0
      if (SCIPisEQ(new_scip, SCIPgetSolVal(new_scip, sol, x[v]), 0.0))
      {
         SCIPfixVar(new_scip, belongsToX[v], 0.0, NULL, NULL);
      }
   }

   //create "is on the right-hand side" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //get solution value of var x[v]
      double value = SCIPgetSolVal(scip, sol, x[v]);
      ScipVar *var = new ScipBinVar(new_scip, "y_" + instance.nodeName[v], -value);
      isOnRHS[v] = var->var;
   }

   //create "valid influencing-set" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //go through all the influencing-sets of v
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         //get the value of infSet[v][i] variable
         double value = SCIPgetSolVal(scip, sol, infSet[v][i].var);

         //give a significative name
         std::stringstream stream;
         for (DNode u : infSet[v][i].nodes)
            stream << instance.nodeName[u];
         string name;
         if (infSet[v][i].nodes.empty())
            name = "u_" + instance.nodeName[v] + "_empty";
         else
            name = "u_" + instance.nodeName[v] + "_" + stream.str();

         ScipVar *var = new ScipBinVar(new_scip, name, value);

         InfluencingSet ini;
         ini.cost = value;
         ini.nodes = infSet[v][i].nodes;
         ini.var = var->var;

         //fix the variable in zero if the value of infSet is zero
         if (SCIPisEQ(new_scip, value, 0.0))
         {
            SCIPfixVar(new_scip, ini.var, 0.0, NULL, NULL);
         }

         validInfSet[v].push_back(ini);
      }
   }

   //create constraint to force a minimum size of two for set X
   //or size (1 - alpha) * n if the lifting is applied.
   //double lb = 2 + (floor((1 - instance.alpha) * instance.n) - 2) * liftingRHS;
   ScipCons *cons= new ScipCons(new_scip, 2.0, SCIPinfinity(new_scip));
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cons->addVar(belongsToX[v], 1.0);
   }
   cons->addVar(liftingRHS, floor((1 - instance.alpha) * instance.n) - 2);
   cons->commit();

   //constraint to decide for the node or the lifting on the right-hand side
   ScipCons *cons2 = new ScipCons(new_scip, 1.0, 1.0);
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cons2->addVar(isOnRHS[v], 1.0);
   }
   cons2->addVar(liftingRHS, 1.0);
   cons2->commit();

   for(DNodeIt v(instance.g); v != INVALID; ++v)
   {
      ScipCons *cons = new ScipCons(new_scip, -SCIPinfinity(new_scip), 0.0);
      cons->addVar(isOnRHS[v], 1.0);
      cons->addVar(belongsToX[v], -1.0);
      cons->commit();
   }

   //constraint to ensure that all variables corresponding to relevant 
   //minimal influencing sets are set to one
   for(DNodeIt v(instance.g); v != INVALID; ++v)
   {
   //go through all the influencing-sets of v
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {

      }
   }

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
   //SCIP_CONSDATA *sourcedata;
   SCIP_CONSDATA *targetdata;

   /* sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL); */
   targetdata = NULL;

   SCIP_CALL(SCIPallocBlockMemory(scip, &targetdata));
   /* targetdata->graph = sourcedata->graph;
   capture_graph(targetdata->graph); */

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
   cout << "CONSSEPALP()" << endl;
   //implement the fischetti's model for separation
   //SCIP_CALL(exactSeparation(scip, conshdlr, conss, nconss, nusefulconss, NULL, result));
   SCIP_CALL(sepaGeneralizedPropCons(scip, conshdlr, NULL, result));

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

   cout << "CONSSEPASOL()" << endl;
   // heuristic separation method for an primal solution
   SCIP_CALL(sepaGeneralizedPropCons(scip, conshdlr, sol, result));

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
   cout << "CONSENFOLP()" << endl;

   //cout << "Construct the suport graph" << endl;
   Digraph new_graph;

   digraph_copy(instance.g, new_graph);

   for (ArcIt a(new_graph); a != INVALID; ++a)
   {
      //if (SCIPgetVarSol(scip, z[a]) < 1 + SCIPepsilon(scip))
      if (SCIPisLT(scip, SCIPgetVarSol(scip, z[a]), 1.0))
      {
         new_graph.erase(a);
      }
   }

   //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph");

   // if a tour was found, we generate a cut constraint
   if (!dag(new_graph))
   {
      *result = SCIP_INFEASIBLE;
      //cout << "isn't acyclic\n";
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
   cout << "CONSELFOPS" << endl;
   *result = SCIP_FEASIBLE;

   //cout << "Construct the suport graph" << endl;
   Digraph new_graph;

   digraph_copy(instance.g, new_graph);

   for (ArcIt a(new_graph); a != INVALID; ++a)
   {
      //if (SCIPgetVarSol(scip, z[a]) < 1 + SCIPepsilon(scip))
      if (SCIPisLT(scip, SCIPgetVarSol(scip, z[a]), 1.0))
      {
         new_graph.erase(a);
      }
   }

   //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph");

   // if a tour is found, the solution must be infeasible
   if (!dag(new_graph))
      *result = SCIP_INFEASIBLE;

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

   cout << "CONSCHECK()\n";
   //cout << "Construct the suport graph" << endl;
   Digraph new_graph;

   digraph_copy(instance.g, new_graph);

   for (ArcIt a(new_graph); a != INVALID; ++a)
   {
      //if (SCIPgetSolVal(scip, sol, z[a]) < 1 + SCIPepsilon(scip))
      if (SCIPisLT(scip, SCIPgetSolVal(scip, sol, z[a]), 1.0))
      {
         new_graph.erase(a);
      }
   }

   //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph");

   // if a subtour is found, the solution must be infeasible
   if (!dag(new_graph))
   {
      *result = SCIP_INFEASIBLE;
      cout << "violation: solution has a subtour\n";
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
   /* for (DNodeIt v(instance.g); v != INVALID; ++v)
      SCIP_CALL(SCIPaddVarLocksType(scip, x[v], locktype, nlockspos, nlocksneg)); */

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
/* SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable *GeneralizedPropagation::clone)
{
   *valid = true;
   return new GeneralizedPropagation(scip, );
} */

/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method, which copies a constraint from one SCIP data structure into a other
 *  SCIP data structure.
 */
SCIP_RETCODE GeneralizedPropagation::createGenPropagationCons(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name)
{
   SCIP_CONSHDLR *conshdlr;
   SCIP_CONSDATA *consdata;

   // find the generalized propagation constraint handler
   conshdlr = SCIPfindConshdlr(scip, "GPC");
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("GPC constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   //create constraint
   SCIP_CALL(
       SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, TRUE, TRUE, TRUE,
                      TRUE, FALSE, FALSE, TRUE, FALSE, FALSE));

   return SCIP_OKAY;
}