#include "GLCIPBase.h"

bool ArcModel::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
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
    DNodeSCIPVarsMap xip(instance.g);
    ArcSCIPVarMap z(instance.g);

    // creates variables x and xip for the incentives of each vertex
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        ScipVar* var  = new ScipBinVar(scip, "x_" + instance.nodeName[v], 0.0);
        x[v] = var->var;

        for(int p = 0; p < instance.incentives[v].size(); p++){
            var  = new ScipBinVar(scip, "x_" + instance.nodeName[v] + "," + to_string(instance.incentives[v][p]), instance.incentives[v][p]);
            xip[v].push_back(var->var);
        }
    }

    // creates variables z for each arc
    for(ArcIt a(instance.g); a != INVALID; ++a){
        ScipVar* var = new ScipIntVar(scip, "z_" + instance.nodeName[instance.g.source(a)] + "," +
            instance.nodeName[instance.g.target(a)], 0.0, 1.0, 0.0);
        z[a] = var->var;
    }

    // add threshold constraints
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        ScipCons *cons = new ScipCons(scip, 0.0, SCIPinfinity(scip));

        // \sum_{p \in P_i} (p - h_v) x_{v, p}
        for (int p = 0; p < instance.incentives[v].size(); ++p){
            cons->addVar(xip[v][p], instance.incentives[v][p]);
        }

        // \sum_{(j, i) \in A} d_{j,i} z_{j, i}
        for(InArcIt a(instance.g, v); a != INVALID; ++a){
            cons->addVar(z[a], instance.influence[a]);
        }

        cons->addVar(x[v], -instance.threshold[v]);
        cons->commit();
    }

    // coupling variables xip and x
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        ScipCons* cons = new ScipCons(scip, 0.0, 0.0);

        for(int p = 0; p < instance.incentives[v].size(); ++p){
            cons->addVar(xip[v][p], 1.0);
        }

        cons->addVar(x[v], -1.0);
        cons->commit();
    }

    // add linking constraints 
    addLinkingConstraints(scip, instance, x, z);
    
    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add small cycle constraints
    addSmallCycleConstraints(scip, instance, x, z);

    //include cycle removal cuts
    //SCIP_CALL( addCuttingPlanes(scip, instance, x, z) );
    
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS* cons;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons, "CycleRemovalCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
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

        for(DNodeIt v(instance.g); v != INVALID; ++v){
            for(int p = 0; p < instance.incentives[v].size(); p++){
                double aux = SCIPgetSolVal(scip, sol, xip[v][p]);

                if(aux > 0.1){
                    cout << "xip[" << instance.nodeName[v] << "," << 
                    instance.incentives[v][p] << "] = " << aux << endl;
                    solution.incentives[v] = instance.incentives[v][p];
                    //cout << "node incentive " << solution.incentives[v] << endl;
                }
                else{
                    solution.incentives[v] = 0.0;
                }
            }
        }
        cout << endl;

        for(ArcIt a(instance.g); a != INVALID; ++a){
            double aux = SCIPgetSolVal(scip, sol, z[a]);

            if(aux > 0.1){
                cout << "z[" << instance.nodeName[instance.g.source(a)] << "," << instance.nodeName[instance.g.target(a)] << "] = " << aux << endl;
                solution.influence[a] = true;
            }
            else{
                solution.influence[a] = false;
            }
        }
    }
}
