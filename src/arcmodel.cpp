#include "arcmodel.h"

// add a cycle founded by the bfs algorithm
void ArcModel::addCycleConstraints(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z, DNodeIntMap &predMap, Arc &backArc){
    DNode s = instance.g.source(backArc);
    DNode t = instance.g.target(backArc);

    // walk from s to t adding corresponding nodes and arcs
    std::vector<Arc> arcs;
    std::vector<DNode> nodes;

    DNode v = s;
    nodes.push_back(v);
    arcs.push_back(backArc);

    do {
        DNode u = instance.g.nodeFromId(predMap[v]);
        Arc a = findArc(instance.g, u, v);
        arcs.push_back(a);
        nodes.push_back(u);
        v = u;
    } while(v != t);

    // now add the corresponding constraints
    for(int i = 0; i < nodes.size(); i++){
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0);

        // add arcs in the cycle
        for(auto a: arcs){
            cons->addVar(z[a], 1.0);
        }

        // add all vertices except one
        for(int j = 0; j < nodes.size(); j++){
            if(j != i){
                DNode v = nodes[j];
                cons->addVar(x[v], -1.0);
            }
        }

        cons->commit();
    }
}

// run a depth first search with maximum height 3
void ArcModel::dfsSmallCycles(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z, DNodeIntMap &colors,
                              DNodeIntMap &predMap, DNode curr, int level, int rootId){
    // visit current node (grey color)
    colors[curr] = 1;

    // pass through current node neighbors
    for(OutArcIt a(instance.g, curr); a != INVALID; ++a){
        DNode next = instance.g.target(a);

        // if next is grey, then we found a cycle
        if(colors[next] == 1){
            addCycleConstraints(scip, instance, x, z, predMap, a);
        }

        // if next is white, then we check the level
        if(colors[next] == 0 && level < 3 && instance.g.id(next) > rootId){
            predMap[next] = instance.g.id(curr);
            dfsSmallCycles(scip, instance, x, z, colors, predMap, next, level + 1, rootId);
        }
    }

    // close this vertex (black color)
    colors[curr] = 2;
}

// add cycle removal constraints for cycles of size up to 4
void ArcModel::addSmallCycleConstraints(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z){
    for(DNodeIt r(instance.g); r != INVALID; ++r){
        // run dfs algorithm for 3 levels
        DNodeIntMap predMap(instance.g);
        DNodeIntMap colors(instance.g);
        mapFill(instance.g, predMap, -1);
        mapFill(instance.g, colors, 0);

        predMap[r] = instance.g.id(r);
        dfsSmallCycles(scip, instance, x, z, colors, predMap, r, 0, instance.g.id(r));
    }
}

/*void ArcModel::addCuttingPlanes(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z){
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS* cons;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons, "CycleRemovalCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));
}*/

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

    // add z and x coupling constraints (z_{s, t} <= x_s)
    for(ArcIt a(instance.g); a != INVALID; ++a){
        DNode s = instance.g.source(a);
        DNode t = instance.g.target(a);
        Arc a2 = findArc(instance.g, t, s);

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

    // add small cycle constraints
    addSmallCycleConstraints(scip, instance, x, z);

    //include cycle removal cuts
    //addCuttingPlanes(scip, instance, x, z);
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
