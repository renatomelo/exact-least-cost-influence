/**@file   GeneralizedPropagationCons.cpp
 * @brief  Generalized propagation constraint handler for GLCIP problems
 */
#include "generalizedpropagationcons.h"
#include <stack>

typedef lemon::Dijkstra<Digraph, ArcValueMap> SptSolver;

struct SCIP_ConsData
{
   Digraph *graph;
   SCIP_ROW *row;
};

//peforms a DFS to verify if there is cycle in the graph suport of the feasible solution
SCIP_Bool findDirectedCycle(
    SCIP *scip,
    SCIP_SOL *sol,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z)
{
   DNodeBoolMap visited(instance.g);
   DNodeBoolMap onStack(instance.g);

   //initialization
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      visited[v] = FALSE;
      onStack[v] = FALSE;
   }

   stack<DNode> s;
   for (DNodeIt u(instance.g); u != INVALID; ++u)
   {
      if (visited[u])
         continue;

      s.push(u);
      while (!s.empty())
      {
         DNode v = s.top();

         if (!visited[v])
         {
            visited[v] = TRUE;
            onStack[v] = TRUE;
         }
         else
         {
            onStack[v] = FALSE;
            s.pop();
         }

         for (OutArcIt a(instance.g, v); a != INVALID; ++a)
         {
            if (SCIPisEQ(scip, SCIPgetSolVal(scip, sol, z[a]), 0.0))
               continue;

            DNode w = instance.g.target(a);
            if (!visited[w])
            {
               s.push(w);
            }
            else if (onStack[w])
            {
               return TRUE;
            }
         }
      }
   }
   return FALSE;
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
    set<DNode> generalizedSet,
    DNode k,
    SCIP_Bool lifting)
{
   //cout << "ADDING CONSTRAINT" << endl;

   // add a constraint for each vertex k in the generilized set (k \in X)
   SCIP_ROW *row;

   if (lifting)
      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "GPC_separation", 1.0, SCIPinfinity(scip), FALSE, TRUE, TRUE));
   else
      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, "GPC_separation", 0.0, SCIPinfinity(scip), FALSE, TRUE, TRUE));

   SCIP_CALL(SCIPcacheRowExtensions(scip, row));

   Phi gpcrow;
   gpcrow.k = k;
   if (!lifting)
   {
      SCIPaddVarToRow(scip, row, x[k], -1.0);
   }
   /*  cout << "elements in X : ";
   for (DNode v : generalizedSet)
   {
      cout << instance.nodeName[v] << " ";
   }
   cout << endl; */

   //SCIP_CALL(SCIPwriteTransProblem(scip, "glcip_transformed.lp", "lp", FALSE));

   //add valid influencing-set variables to the row
   for (DNode v : generalizedSet)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (!GLCIPBase::intersects(generalizedSet, infSet[v][i].nodes))
         {
            SCIPaddVarToRow(scip, row, infSet[v][i].var, 1.0);
         }
      }
      gpcrow.generalizedSet.insert(v);
   }

   SCIP_CALL(SCIPflushRowExtensions(scip, row));

   // add cut
   //if (SCIPisCutEfficacious(scip, sol, row))
   {
      SCIP_Bool infeasible;
      SCIP_CALL(SCIPaddRow(scip, row, TRUE, &infeasible));
      //SCIPprintRow(scip, row, NULL);

      if (infeasible)
         *result = SCIP_CUTOFF;
      else
      {
         //save the added constraint to use in the price
         gpcrow.row = row;
         //if (gpcrow.k != INVALID) // save only the rows that have an associated k
         gpcrows.push_back(gpcrow);
      }
   }
   //SCIP_CALL(SCIPreleaseRow(scip, &row));

   return SCIP_OKAY;
}

SCIP_RETCODE GeneralizedPropagation::greedSetExtensionHeur(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr, //the constraint handler itself
    SCIP_SOL *sol,           //primal solution that should be separated
    SCIP_RESULT *result      //pointer to store the result of the separation call
)
{
   assert(result != NULL);
   *result = SCIP_DIDNOTFIND;

   set<DNode> generalizedSet;
   int count = 0;
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cout << "count " << count++ << endl;
      DNode k = INVALID;
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         //get the empty influensin-set var
         if (infSet[v][i].nodes.empty())
         {
            double value = SCIPgetSolVal(scip, sol, infSet[v][i].var);
            // x[v] > 0
            if ((SCIPgetSolVal(scip, sol, x[v]) > SCIPepsilon(scip)) && SCIPisLT(scip, value, SCIPgetSolVal(scip, sol, x[v])))
            {
               //add v to X
               generalizedSet.insert(v);
               k = v;
               cout << "k is the vertex: " << instance.nodeName[k] << endl;
               //break;
            }
         }
      }

      //build a candidate list
      set<DNode> candidates;
      for (DNodeIt u(instance.g); u != INVALID; ++u)
      {
         if (generalizedSet.count(u) == 0)
         {
            //add to the candidate list every vertex that has arc to a element of X
            for (OutArcIt a(instance.g, u); a != INVALID; ++a)
            {
               DNode w = instance.g.target(a);
               if (generalizedSet.count(w) != 0)
               {
                  candidates.insert(u);
                  //cout << "vertex " << instance.nodeName[u] << " added to the candidate list\n";
                  break;
               }
            }
         }
      }

      //find node with the least contribution to LHS of GPC
      double min = SCIPinfinity(scip);
      DNode l = INVALID;
      for (DNode j : candidates)
      {
         double sum = 0;
         for (unsigned int p = 0; p < infSet[j].size(); p++)
         {
            if (!GLCIPBase::intersects(generalizedSet, infSet[j][p].nodes))
            {
               sum += SCIPgetSolVal(scip, sol, infSet[j][p].var);
            }
         }

         //if (sum < min)
         if (SCIPisLT(scip, sum, min))
         {
            min = sum;
            l = j;
         }
      }

      //add to X the candidate node which leads to the smallest LHS of GPC
      if (l != INVALID && k != INVALID)
      {
         cout << "adding node " << instance.nodeName[l] << " to X\n";
         generalizedSet.insert(l);

         //check whether GPC for set X and node k is violated, in which case we add it to the model
         double sum = 0;
         for (DNode j : generalizedSet)
         {
            for (unsigned int p = 0; p < infSet[j].size(); p++)
            {
               if (!GLCIPBase::intersects(generalizedSet, infSet[j][p].nodes))
               {
                  sum += SCIPgetSolVal(scip, sol, infSet[j][p].var);
               }
            }
         }
         cout << "sum = " << sum << ", value of x[k] = " << SCIPgetSolVal(scip, sol, x[k]) << endl;

         if (SCIPisLT(scip, sum, SCIPgetSolVal(scip, sol, x[k])))
         {
            *result = SCIP_SEPARATED;

            cout << "adding violated inequality\n";
            addGeneralizedPropCons(scip, conshdlr, sol, result, generalizedSet, k, FALSE);

            // node has become infeasible
            if (*result == SCIP_CUTOFF)
            {
               cout << "infeasible node\n";
               return SCIP_OKAY;
            }
         }
      }
   }

   return SCIP_OKAY;
}

SCIP_RETCODE GeneralizedPropagation::sepaGeneralizedPropCons(
    SCIP *scip,
    SCIP_CONSHDLR *conshdlr, //the constraint handler itself
    SCIP_SOL *sol,           //primal solution that should be separated
    SCIP_RESULT *result      //pointer to store the result of the separation call
)
{
   //cout << "HEURISITC SEPARATION()" << endl;

   assert(result != NULL);
   *result = SCIP_DIDNOTFIND;

   ArcValueMap weight(instance.g);

   // construct edges weights according to w[(i, j)] = x_i - z_{i,j}
   for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      weight[a] = SCIPgetSolVal(scip, sol, x[instance.g.source(a)]) - SCIPgetSolVal(scip, sol, z[a]);
   }

   for (DNodeIt u(instance.g); u != INVALID; ++u)
   {
      if (SCIPisGT(scip, SCIPgetSolVal(scip, sol, x[u]), 0.0))
      {
         SptSolver spt(instance.g, weight);
         spt.run(u);

         for (InArcIt a(instance.g, u); a != INVALID; ++a)
         {
            DNode v = instance.g.source(a);

            //TODO add the constraint only the cycle that violates the GPC
            // in order to violate cycle inequalities we need that L_{u, v} + w_{v, u} < x_u
            //if (spt.reached(v) && (spt.dist(v) + weight[a] < SCIPgetSolVal(scip, sol, x[u])))
            if (spt.reached(v))
            {
               //cout << "cycle found: ";
               set<DNode> cycle;
               for (DNode w = v; w != u; w = spt.predNode(w))
               {
                  if (w != INVALID && spt.reached(w))
                  {
                     cycle.insert(w);
                     //cout << instance.nodeName[w] << ", ";
                  }
               }
               cycle.insert(u);
               //cout << instance.nodeName[u] << endl;

               //choose the maximum x[k]
               double max_k = -1;
               int idOfNodek = -1;
               for (DNode i : cycle)
               {
                  //cout << "value of x_v: " << SCIPgetSolVal(scip, sol, x[v]) << endl;
                  if (SCIPgetSolVal(scip, sol, x[i]) > max_k)
                  {
                     idOfNodek = instance.g.id(i);
                     max_k = SCIPgetSolVal(scip, sol, x[i]);
                     if (SCIPisEQ(scip, max_k, 1))
                     {
                        break;
                     }
                  }
               }
               assert(idOfNodek != -1);
               DNode k = instance.g.nodeFromId(idOfNodek);

               //check whether GPC for the cycle and node k is violated, in which case we add it to the model
               double sum = 0;
               for (DNode j : cycle)
               {
                  for (unsigned int p = 0; p < infSet[j].size(); p++)
                  {
                     if (!GLCIPBase::intersects(cycle, infSet[j][p].nodes))
                     {
                        sum += SCIPgetSolVal(scip, sol, infSet[j][p].var);
                     }
                  }
               }
               //cout << "sum = " << sum << ", value of x[k] = " << SCIPgetSolVal(scip, sol, x[k]) << endl;

               if (SCIPisLT(scip, sum, SCIPgetSolVal(scip, sol, x[k])))
               {
                  *result = SCIP_SEPARATED;

                  //cout << "adding violated inequality\n";
                  addGeneralizedPropCons(scip, conshdlr, sol, result, cycle, k, FALSE);

                  // node has become infeasible
                  if (*result == SCIP_CUTOFF)
                  {
                     cout << "infeasible node\n";
                     return SCIP_OKAY;
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

void printFractionalSol(SCIP *scip,
                        GLCIPInstance &instance,
                        SCIP_SOL *sol,
                        DNodeSCIPVarMap &x,
                        ArcSCIPVarMap &z,
                        DNodeInfSetsMap &infSet)
{
   cout << "\nFractional solution: \n";
   /*  cout << "x variables\n";
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cout << "x[" << instance.nodeName[v] << "] = " << SCIPgetSolVal(scip, sol, x[v]) << endl;
   } */

   /* cout << "lambda variables\n";
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         cout << SCIPvarGetName(infSet[v][i].var) << " = " << SCIPgetSolVal(scip, sol, infSet[v][i].var) << endl;
      }
   } */

   cout << "LHS of the new constraint: \n";
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      double sum = 0;
      for (InArcIt a(instance.g, v); a != INVALID; ++a)
      {
         DNode u = instance.g.source(a);
         for (unsigned int i = 0; i < infSet[v].size(); i++)
         {
            if (infSet[v][i].nodes.count(u))
            {
               sum += SCIPgetSolVal(scip, sol, infSet[v][i].var);
            }
         }
      }
      InDegMap<Digraph> inDeg(instance.g);

      cout << instance.nodeName[v] << ": degree =  " << inDeg[v] 
            << ", threshold = " << instance.threshold[v] << ", lhs  = " << sum << endl;
   }
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
   assert(result != NULL);
   *result = SCIP_DIDNOTFIND;

   //show fractional solution
   //printFractionalSol(scip, instance, sol, x, z, infSet);

   SCIP *new_scip = NULL;

   // initialize SCIP enviroment
   SCIP_CALL(SCIPcreate(&new_scip));

   SCIP_CALL(SCIPincludeDefaultPlugins(new_scip));
   SCIP_CALL(SCIPsetSeparating(new_scip, SCIP_PARAMSETTING_OFF, TRUE));
   SCIPsetPresolving(new_scip, SCIP_PARAMSETTING_OFF, TRUE);

   // create empty problem
   SCIP_CALL(SCIPcreateProb(new_scip, "GPC_separation", 0, 0, 0, 0, 0, 0, 0));
   SCIP_CALL(SCIPsetIntParam(new_scip, "display/verblevel", 1));

   // declaring variables
   DNodeSCIPVarMap belongsToX(instance.g);  //indicates membership of nodes to set X
   DNodeSCIPVarMap isOnRHS(instance.g);     //equal to one node if i is on the RHS of GPC
   DNodeInfSetsMap validInfSet(instance.g); //equal to one iff node i \in X and U has no
                                            // intersection with X
   //decides whether the lifting to one on the right-hand side is applied or not
   ScipVar *var_lifting = new ScipBinVar(new_scip, "lifting_rhs", -1);
   SCIP_VAR *liftingRHS = var_lifting->var;

   //create "belong to X" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      ScipVar *var = new ScipBinVar(new_scip, "v_" + instance.nodeName[v], 0.0);
      belongsToX[v] = var->var;
   }

   //create "is on the right-hand side" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //get solution value of var x[v]
      double value = SCIPgetSolVal(scip, sol, x[v]);
      ScipVar *var = new ScipBinVar(new_scip, "y_" + instance.nodeName[v], -value);
      isOnRHS[v] = var->var;
   }

   int sum = 0;
   //create "valid influencing-set" variables
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //go through all the influencing-sets of v
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         //get the value of infSet[v][i] variable
         double value = SCIPgetSolVal(scip, sol, infSet[v][i].var);

         //it is sufficient to consider infSet variables greater than zero
         if (SCIPisEQ(new_scip, value, 0.0))
            continue;

         //give a significative name
         std::stringstream stream;
         for (DNode u : infSet[v][i].nodes)
         {
            stream << instance.nodeName[u] + ",";
         }
         string name;
         if (infSet[v][i].nodes.empty())
            name = "u_" + instance.nodeName[v] + "_empty";
         else
            name = "u_" + instance.nodeName[v] + "_{" + stream.str() + "}";

         ScipVar *var = new ScipBinVar(new_scip, name, value);

         InfluencingSet ini;
         ini.cost = value;
         for (DNode u : infSet[v][i].nodes)
         {
            ini.nodes.insert(u);
         }
         ini.var = var->var;

         validInfSet[v].push_back(ini);
      }
      sum += validInfSet[v].size();
   }

   //cout << "number of 'valid inf set' variables = " <<  sum << endl;

   //create constraint to force a minimum size of two for set X
   //or size (1 - alpha) * n if the lifting is applied.
   ScipCons *cons1 = new ScipCons(new_scip, 2.0, SCIPinfinity(new_scip), "size-of-X cons");
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cons1->addVar(belongsToX[v], 1.0);
   }
   cons1->addVar(liftingRHS, -floor((1 - instance.alpha) * instance.n) - 2);
   cons1->commit();

   //constraint to decide for the node or the lifting on the right-hand side
   ScipCons *cons2 = new ScipCons(new_scip, 1.0, 1.0, "decide-lifting cons");
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      cons2->addVar(isOnRHS[v], 1.0);
   }
   cons2->addVar(liftingRHS, 1.0);
   cons2->commit();

   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      ScipCons *cons = new ScipCons(new_scip, -SCIPinfinity(new_scip), 0.0, "node-in-rhs cons");
      cons->addVar(isOnRHS[v], 1.0);
      cons->addVar(belongsToX[v], -1.0);
      cons->commit();
   }

   //constraint to ensure that all variables corresponding to relevant
   //minimal influencing sets are set to one
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      //go through all the influencing-sets of v
      for (unsigned int i = 0; i < validInfSet[v].size(); i++)
      {
         ScipCons *cons = new ScipCons(new_scip, 0.0, SCIPinfinity(new_scip), "valid-inf-sets cons");
         cons->addVar(belongsToX[v], -1.0);

         for (DNode u : validInfSet[v][i].nodes)
         {
            cons->addVar(belongsToX[u], 1.0);
         }

         cons->addVar(validInfSet[v][i].var, 1.0);
         cons->commit();
      }
   }

   SCIP_CALL(SCIPsolve(new_scip));

   //get the vertices that belgong to X
   SCIP_SOL *localSol = SCIPgetBestSol(new_scip);

   set<DNode> generalizedSet;
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      double value = SCIPgetSolVal(new_scip, localSol, belongsToX[v]);
      if (value > 0.5)
      {
         generalizedSet.insert(v);
      }
   }
   //cout << "lifting the RHS = " << SCIPgetSolVal(new_scip, localSol, liftingRHS) << "\n";

   //check if solution value is negative (and hence a violated constraint was found)
   if (SCIPisNegative(new_scip, SCIPgetPrimalbound(new_scip)))
   {
      *result = SCIP_SEPARATED;

      for (DNode j : generalizedSet)
      {
         cout << instance.nodeName[j] << " ";
      }
      cout << endl;

      // in case the lifting isn't done find the vertex to be in the RHS
      if (SCIPgetSolVal(new_scip, localSol, liftingRHS) < 0.5)
      {
         for (DNodeIt v(instance.g); v != INVALID; ++v)
         {
            if (SCIPgetSolVal(new_scip, localSol, isOnRHS[v]) > .5)
            {
               //cout << instance.nodeName[v] << " is on the RHS\n";
               addGeneralizedPropCons(scip, conshdlr, sol, result, generalizedSet, v, FALSE);
            }
         }
      }
      else
         addGeneralizedPropCons(scip, conshdlr, sol, result, generalizedSet, INVALID, TRUE);

      // node has become infeasible
      if (*result == SCIP_CUTOFF)
      {
         cout << "node has become infeasible\n";
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

// check if current solution respect cycle constraints
bool isValid(SCIP *scip, SCIP_SOL *sol, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z)
{
   ArcValueMap weight(instance.g);

   // construct edges weights according to w[(i, j)] = x_i - z_{i,j}
   for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      // compute weight in support graph
      weight[a] = SCIPgetSolVal(scip, sol, x[instance.g.source(a)]) - SCIPgetSolVal(scip, sol, z[a]);
   }

   // for each node u in the graph, we will run dijkstra to get the shortest path from u to its neighbors
   for (DNodeIt u(instance.g); u != INVALID; ++u)
   {
      if (SCIPgetSolVal(scip, sol, x[u]) > 1 - SCIPepsilon(scip))
      {
         // run dijkstra starting from u
         SptSolver spt(instance.g, weight);
         spt.run(u);

         // for each neighbor of u, we check if it violates the cycle inequality
         for (InArcIt a(instance.g, u); a != INVALID; ++a)
         {
            DNode v = instance.g.source(a);

            // in order to violate cycle inequalities we need that L_{u, v} + w_{v, u} < x_u
            if (spt.reached(v) && spt.dist(v) + weight[a] < SCIPgetSolVal(scip, sol, x[u]))
            {
               //cout << "not valid!" << endl;
               return false;
            }
         }
      }
   }

   //cout << "valid!" << endl;
   return true;
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
   //cout << "CONSTRANS()\n";

   /* create target constraint */
   SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, NULL,
                            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
                            SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
                            SCIPconsIsStickingAtNode(sourcecons)));

   return SCIP_OKAY;
}

SCIP_RETCODE printRows(SCIP *scip)
{
   /* get data */
   SCIP_ROW **rows;
   int nrows;

   SCIP_CALL(SCIPgetLPRowsData(scip, &rows, &nrows));
   assert(nrows > 0);
   assert(rows != NULL);
   cout << "number of rows after add a cut: " << nrows << endl;

   string name = "GPC_separation";
   for (int i = 0; i < nrows; i++)
   {
      if (SCIProwGetName(rows[i]) == name)
      {
         SCIPprintRow(scip, rows[i], NULL);
      }
   }

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
   cout << "(GPC) CONSSEPALP()" << endl;

   //implement the fischetti's model for separation
   SCIP_CALL(exactSeparation(scip, conshdlr, NULL, result));
   //SCIP_CALL(sepaGeneralizedPropCons(scip, conshdlr, NULL, result));

   //printRows(scip);

   //cout << "--------END OF SEPARATION---------\n";
   return SCIP_OKAY;
}

/** TODO maybe not necessary but verify 
 * separation method of constraint handler for arbitrary primal solution
 *
 *  Separates all constraints of the constraint handler. 
 *  The method is called outside the LP solution loop (e.g., by
 *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *  Instead, the method should produce cuts that separate the given solution.
 *
 */
SCIP_DECL_CONSSEPASOL(GeneralizedPropagation::scip_sepasol)
{
   cout << "CONSSEPASOL()" << endl;
   // heuristic separation method for an primal solution
   SCIP_CALL(exactSeparation(scip, conshdlr, sol, result));

   return SCIP_OKAY;
}

void getSuportGraph(
    SCIP *scip,
    GLCIPInstance &instance,
    SCIP_SOL *sol,
    ArcSCIPVarMap &z,
    Digraph &new_graph)
{
   DNodeDNodeMap nodeRef(new_graph);
   ArcArcMap arcRef(new_graph);
   digraphCopy(instance.g, new_graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

   for (ArcIt a(new_graph); a != INVALID; ++a)
   {
      if (SCIPisEQ(scip, SCIPgetSolVal(scip, sol, z[arcRef[a]]), 0))
      {
         new_graph.erase(a);
      }
   }

   //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph", nodeRef);
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

   //construct the suport graph
   Digraph new_graph;

   //getSuportGraph(scip, instance, NULL, z, new_graph);
   DNodeDNodeMap nodeRef(new_graph);
   ArcArcMap arcRef(new_graph);
   digraphCopy(instance.g, new_graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

   for (ArcIt a(new_graph); a != INVALID; ++a)
   {
      if (SCIPisEQ(scip, SCIPgetVarSol(scip, z[arcRef[a]]), 0))
      {
         new_graph.erase(a);
      }
   }

   //mostrar variaveis lambdas
   /* for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      for (unsigned int i = 0; i < infSet[v].size(); i++)
      {
         if (SCIPisPositive(scip, SCIPgetVarSol(scip, infSet[v][i].var)))
         {
            cout << SCIPvarGetName(infSet[v][i].var) << "\t = "
                 << SCIPgetVarSol(scip, infSet[v][i].var) << endl;
         }
      }
   } */

   //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph", nodeRef);

   // if a cycle is found, the solution must be infeasible
   if (!dag(new_graph))
   {
      *result = SCIP_INFEASIBLE;
      //cout << "(SUPORT GRAPH) solution has a cycle\n";

      SCIP_CALL(sepaGeneralizedPropCons(scip, conshdlr, NULL, result));

      if (*result == SCIP_DIDNOTFIND)
      {
         cout << "SCIP_DIDNOTFIND\n";
         *result = SCIP_INFEASIBLE;
         //SCIP_CALL(SCIPwriteTransProblem(scip, "glcip_transformed.lp", "lp", FALSE));
      }
      //cout << "number of pseudo candidates = " << SCIPgetNPseudoBranchCands(scip) << endl;
      /* if (SCIPgetNPseudoBranchCands(scip) == 0)
      {
         // if you just return INFEASIBLE without any branching candidates available,
         // SCIP will have no chance to resolve this infeasibility. Therefore, you
         // should return SCIP_CUTOFF, meaning that the node with its current
         // variable bounds (or fixings) contains no feasible solution.
         cout << "SCIP_CUTOFF\n";
         *result = SCIP_CUTOFF;
         //exit(0);
      } */
   }

   /*    if (findDirectedCycle(scip, NULL, instance, x, z))
   {
      //cout << "(DFS) violation: solution has a cycle\n";
      exactSeparation(scip, conshdlr, NULL, result);
      *result = SCIP_INFEASIBLE;
   } */

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

   getSuportGraph(scip, instance, NULL, z, new_graph);

   // if a cycle is found, the solution must be infeasible
   if (!dag(new_graph))
   {
      *result = SCIP_INFEASIBLE;
      //cout << "(COPIED GRAPH) isn't acyclic\n";
   }

   /*  if (findDirectedCycle(scip, NULL, instance, x, z))
   {
      //cout << "(DFS) violation: solution has a cycle\n";
      exactSeparation(scip, conshdlr, NULL, result);
      *result = SCIP_INFEASIBLE;
   } */

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

   //construct the suport graph
   Digraph new_graph;
   getSuportGraph(scip, instance, sol, z, new_graph);

   // if a cycle is found, the solution must be infeasible
   if (!dag(new_graph))
   {
      *result = SCIP_INFEASIBLE;
      cout << "(SUPORT GRAPH) violation: solution has a cycle\n";
   }

   /*  ArcValueMap weight(instance.g);
      for (ArcIt a(instance.g); a != INVALID; ++a)
      {
         weight[a] = SCIPgetSolVal(scip, sol, z[a]);
      }
      GraphViewer::ViewGLCIPFracSolution(instance, weight, "Fractional Sol"); */

   //if (findDirectedCycle(scip, sol, instance, x, z))
   /* if (!isValid(scip, NULL, instance, x, z))
   {
      cout << "(DFS) violation: solution has a cycle\n";
      *result = SCIP_INFEASIBLE;
   } */

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
SCIP_DECL_CONSPROP(GeneralizedPropagation::scip_prop)
{
   cout << "CONSPROP()\n";
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
   //cout << "CONSLOCK\n";

   // round z up may cause cycle inequalities to be infeasible
   for (ArcIt a(instance.g); a != INVALID; ++a)
   {
      SCIP_CALL(SCIPaddVarLocksType(scip, z[a], locktype, nlocksneg, nlockspos));
   }

   // round x down may cause cycle inequalities to be infeasible
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      SCIP_CALL(SCIPaddVarLocksType(scip, x[v], locktype, nlockspos, nlocksneg));
   }

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

SCIP_DECL_CONSPRINT(GeneralizedPropagation::scip_print)
{
   SCIPinfoMessage(scip, file, "cycle of Graph G with %d nodes and %d edges\n", instance.n, instance.m);

   return SCIP_OKAY;
}

/** create genatalized propagation constraint handler */
SCIP_RETCODE GeneralizedPropagation::createGenPropagationCons(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name)
{
   SCIP_CONSHDLR *conshdlr;
   SCIP_CONSDATA *consdata = NULL;

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
                      TRUE, FALSE, TRUE, FALSE, TRUE, FALSE));

   return SCIP_OKAY;
}