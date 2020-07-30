#include "GLCIPBase.h"
#include "heur_dualbound.h"
#include "heur_ordering.h"

bool ArcModel::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    //SCIP variables and initialization
    SCIP *scip;
    SCIP_CALL(SCIPcreate(&scip));

    //set some parameters
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    // create an empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP Problem", NULL, NULL, NULL, NULL, NULL, NULL, NULL));
    //SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 0));
    SCIP_CALL(SCIPsetStringParam(scip, "visual/vbcfilename", "branchandbound.vbc"));

    //SCIP_CALL(SCIPsetBoolParam(scip, "lp/presolving", FALSE));
    SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);

    // add variables to the model
    DNodeSCIPVarMap x(instance.g);
    DNodeSCIPVarsMap xip(instance.g);
    ArcSCIPVarMap z(instance.g);

    // creates variables x and xip for the incentives of each vertex
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipVar *var = new ScipBinVar(scip, "x_" + instance.nodeName[v], 0.0);
        x[v] = var->var;

        for (unsigned int p = 0; p < instance.incentives[v].size(); p++)
        {
            var = new ScipBinVar(scip, "x_" + instance.nodeName[v] + "," + to_string(instance.incentives[v][p]), instance.incentives[v][p]);
            xip[v].push_back(var->var);
            //std::cout << "x_" + instance.nodeName[v] + "," + to_string(instance.incentives[v][p]) << endl;
        }
    }

    // creates variables z for each arc
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        ScipVar *var = new ScipIntVar(scip, "z_" + instance.nodeName[instance.g.source(a)] + "," + instance.nodeName[instance.g.target(a)], 0.0, 1.0, 0.0);
        z[a] = var->var;
    }

    // add threshold constraints
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipCons *cons = new ScipCons(scip, 0.0, SCIPinfinity(scip), "threshold cons");

        // \sum_{p \in P_i} (p - h_v) x_{v, p}
        for (unsigned int p = 0; p < instance.incentives[v].size(); p++)
        {
            cons->addVar(xip[v][p], instance.incentives[v][p]);
        }

        // \sum_{(j, i) \in A} d_{j,i} z_{j, i}
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            cons->addVar(z[a], instance.influence[a]);
        }

        cons->addVar(x[v], -instance.threshold[v]);
        cons->commit();
    }

    // coupling variables xip and x
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipCons *cons = new ScipCons(scip, 0.0, 0.0, "linking cons");

        for (unsigned int p = 0; p < instance.incentives[v].size(); p++)
        {
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
    //addSmallCycleConstraints(scip, instance, x, z);
    addAllSmallDirectedCycles(scip, instance, x, z);

    //include cycle removal cuts
    //SCIP_CALL( addCuttingPlanes(scip, instance, x, z) );
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS *cons;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons, "cycle-elimination", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));

    //primal heuristic
    HeurOrdering *ordering = new HeurOrdering(scip, instance, x, z, xip);
    SCIP_CALL(SCIPincludeObjHeur(scip, ordering, TRUE));

    // bound the execution time
    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));

    //SCIP tries to solve the LP
    SCIP_CALL(SCIPsolve(scip));
    //SCIP_CALL(SCIPprintStatistics(scip, NULL));
    //SCIP_CALL(SCIPprintOrigProblem(scip, NULL, NULL, FALSE));

    //reached time limit
    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT)
    {
        //cout << "reached time limit" << endl;
        printf("%.2lf\t%lld\t%d\t%lf\t%lf\t%.2lf\n", SCIPgetSolvingTime(scip), 
                                             SCIPgetNNodes(scip),
                                             SCIPgetNContVars(scip), 
                                             SCIPgetDualbound(scip), 
                                             SCIPgetPrimalbound(scip),
                                             SCIPgetGap(scip));
        return 0;
    }

    //std::cout << SCIPgetSolvingTime(scip) << std::endl;
    //cout << "time \tnodes \tdualbound \tprimalbound \tgap" << endl;
    printf("%.2lf\t%lld\t%d\t%lf\t%lf\t%.2lf\n", SCIPgetSolvingTime(scip), 
                                             SCIPgetNNodes(scip),
                                             SCIPgetNContVars(scip), 
                                             SCIPgetDualbound(scip), 
                                             SCIPgetPrimalbound(scip),
                                             SCIPgetGap(scip));

    //founded optimal solution, now we need to construct the solution
    /* else
    {
        // get measures
        SCIP_SOL *sol = SCIPgetBestSol(scip);

        for (DNodeIt v(instance.g); v != INVALID; ++v)
        {
            solution.incentives[v] = 0.0;
            for (unsigned int p = 0; p < instance.incentives[v].size(); p++)
            {
                double aux = SCIPgetSolVal(scip, sol, xip[v][p]);

                if (aux > 0.1)
                {
                    cout << "xip[" << instance.nodeName[v] << "," << instance.incentives[v][p] << "] = " << aux << endl;
                    solution.incentives[v] = instance.incentives[v][p];
                    //cout << "node incentive " << solution.incentives[v] << endl;
                }
            }
        }
        cout << endl;

        for (ArcIt a(instance.g); a != INVALID; ++a)
        {
            double aux = SCIPgetSolVal(scip, sol, z[a]);

            if (aux > 0.1)
            {
                cout << "z[" << instance.nodeName[instance.g.source(a)] << "," << instance.nodeName[instance.g.target(a)] << "] = " << aux << endl;
                solution.influence[a] = true;
            }
            else
            {
                solution.influence[a] = false;
            }
        }
    } */

    return SCIP_OKAY;
}
