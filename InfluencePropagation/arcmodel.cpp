#include "arcmodel.h"

void ArcModel::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    //set initial clock
    double elapsed_time;
    clock_t begin = clock();

    //---------------------------------------------------------------------------
    //SCIP variables and initialization
    SCIP *scip;
    SCIP_CALL(SCIPcreate(&scip));

    //set some parameters
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPsetIntParam(scip, "nodeselection/dfs/stdpriority", 1073741823));
    //SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 5));
    //SCIP_CALL(SCIPsetBoolParam(scip, "display/lpinfo", TRUE));
    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrestarts", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
    SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);
    SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE);
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanupcols", TRUE));
    SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanupcolsroot", TRUE));
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/colagelimit", -1));
    //SCIP_CALL(SCIPsetIntParam(scip, "lp/rowagelimit", -1));
    //SCIP_CALL(SCIPsetIntParam(scip, "separating/cutagelimit",  -1));
    //SCIP_CALL(SCIPsetIntParam(scip, "constraints/agelimit",  -1));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/epsilon", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/feastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/lpfeastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/dualfeastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "separating/minefficacy", 0.001));
    //SCIP_CALL(SCIPsetStringParam(scip, "visual/vbcfilename", "mytree.vbc"));

    // create an empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP Problem", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    // add variables to the model
    DNodeSCIPVarMap x(instance.g);
    ArcSCIPVarMap z(instance.g);

    for(DNodeIt v(instance.g); v != INVALID; ++v){
        ScipVar* var;
        var = new ScipIntVar(scip, "x_" + to_string(instance.g.id(v)), 0.0, 1.0, instance.incentives[v]);
        x[v] = var->var;
    }

    for(ArcIt a(instance.g); a != INVALID; ++a){
        ScipVar* var;
        var = new ScipIntVar(scip, "z_" + to_string(instance.g.id(instance.g.source(a))) + "," + to_string(instance.g.id(instance.g.target(a))),
                             0.0, 1.0, 0.0);
        z[a] = var->var;
    }

    // add threshold constraints
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        ScipCons *cons = new ScipCons(scip, 0.0, SCIPinfinity(scip));

        // \sum_{p \in P_i} p x_{v, p}
        cons->addVar(x[v], instance.incentives[v]);

        // \sum_{(j, i) \in A} d_{j,i} z_{j, i}
        for(InArcIt a(instance.g, v); a != INVALID; ++a){
            cons->addVar(z[a], instance.influence[z]);
        }

        // - h_v x_v
        cons->addVar(x[v], -instance.threshold[v]);

        cons->commit();
    }

    // add z and x coupling constraints
    for(ArcIt a(instance.g); a != INVALID; ++a){
        DNode s = instance.g.source(a);
        DNode t = instance.g.target(a);
        Arc a2 = findArc(g, t, s);

        if(a2 == INVALID){
            ScipCons *cons = new ScipCons(scip, 0.0, SCIPinfinity(scip));

            cons->addVar(x[s], 1.0);
            cons->addVar(z[a], -1.0);

            cons->commit();
        }
    }

    // add alpha constraints
    ScipCons *consAlpha = new ScipCons(scip, ceil(instance.alpha * instance.n), SCIPinfinity(scip));
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        consAlpha->addVar(x[v], 1.0);
    }
    consAlpha->commit();

    //include cycle removal cuts (TODO)
    SCIPCutsCallback cuts = SCIPCutsCallback(scip, instance, x, varPool);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS* cons;
    SCIP_CALL(cuts.SCIPcreateSteinerCuts(scip, &cons, "STPCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));

    // bound the execution time
    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));

    //SCIP tries to solve the LP
    SCIP_CALL(SCIPsolve(scip));
    //SCIP_CALL(SCIPprintStatistics(scip, NULL));
    //SCIP_CALL(SCIPprintOrigProblem(scip, NULL, NULL, FALSE));

    //reached time limit
    if(SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT){
        cout << "reached time limit" << endl;
        return 0;
    }

    //founded optimal solution, now we need to construct the solution
    else{
        // get measures
        SCIP_SOL* sol = SCIPgetBestSol(scip);
        measures.solutionCost = SCIPgetSolOrigObj(scip, sol);

        for(ArcIt a(instance.g); a != INVALID; ++a){
            double aux;
            if(instance.shouldPrice){
                aux = varPool.getArcValue(scip, sol, a);
            }
            else{
                aux = SCIPgetSolVal(scip, sol, x[a]);
            }

            if(aux > 0.1){
                cout << "x[" << instance.g.id(instance.g.source(a)) << "," << instance.g.id(instance.g.target(a)) << "] = " << aux << endl;
                solution[a] = true;
            }
            else{
                solution[a] = false;
            }
        }
    }
}
