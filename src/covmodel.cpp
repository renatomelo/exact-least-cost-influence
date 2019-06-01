#include "GLCIPBase.h"
#include <map>
#include "pricer_glcip.h"

bool CovModel::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    Digraph &graph = instance.g;
    SCIP* scip = NULL;

    // initialize SCIP enviroment
    SCIP_CALL( SCIPcreate(&scip) );

    // explicitely enables the use of debug solution for this scip instance
    SCIPenableDebugSol(scip);
    
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 5));

    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));

    // create empty problem
    SCIP_CALL( SCIPcreateProb(scip, "GLCIP_ColGeneration", 0, 0, 0, 0, 0, 0, 0) );

    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrestarts", 0));
    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));
    SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);
    SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE);
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanupcols", TRUE));
    SCIP_CALL(SCIPsetBoolParam(scip, "lp/cleanupcolsroot", TRUE));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/epsilon", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/feastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/lpfeastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "numerics/dualfeastol", 0.0001));
    SCIP_CALL(SCIPsetRealParam(scip, "separating/minefficacy", 0.001));

    DNodeSCIPVarMap x(graph); // active-vertex variables
    ArcSCIPVarMap z(graph); // arc-influence variables
    DNodeSCIPVarsMap infSet(graph); // influencing-set variables

    // create variables x
    for (DNodeIt v(graph); v != INVALID; ++v){
        ScipVar* var = new ScipBinVar(scip, "x_" + instance.nodeName[v], 0);  
        x[v] = var->var;
    }

    // create variables z
    for (ArcIt a(graph); a != INVALID; ++a){
        DNode u = graph.source(a);
        DNode v = graph.target(a);
        ScipVar* var = new ScipBinVar(scip, "z_" + instance.nodeName[u] + "" + instance.nodeName[v], 0);
        z[a] = var->var;
    }

    //TODO create influencing-set variables (dynamically)

    // add linking constraints 
    addLinkingConstraints(scip, instance, x, z);
    
    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    addSmallCycleConstraints(scip, instance, x, z);

    // add vertex-coverage constraints
    DNodeConsMap vertCons(graph);
    for (DNodeIt v(graph); v != INVALID; ++v){
        ScipConsPrice* cons = new ScipConsPrice(scip, 0, SCIPinfinity(scip));
        cons->addVar(x[v], -1);
        vertCons[v] = cons->cons;
        cons->commit();
    }

    // add arc-coverage constraints
    ArcConsMap arcCons(graph);
    for (ArcIt a(graph); a != INVALID; ++a){
        ScipConsPrice* cons = new ScipConsPrice(scip, 0, SCIPinfinity(scip));
        cons->addVar(z[a], 1);
        arcCons[a] = cons->cons;
        cons->commit();
    }

    //Pricing to generate the propagation constraints the chosen arcs constraints 
    //TODO decide if we start with a initial solution obtained heuristically
    
    // include pricer
    static const char* PRICER_NAME = "GLCIP_pricer";
    ObjPricerGLCIP* pricer = new ObjPricerGLCIP(scip, PRICER_NAME, instance, z, x, arcCons, vertCons);

    SCIP_CALL( SCIPincludeObjPricer(scip, pricer, TRUE) );

    // activate pricer
    SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME)) );
    std::cout << "\n ACTIVATE PRICER WAS CALLED DHFLKAHFLAKSDHÇFLHSAÇFLÇALSDFLÇASFÇLSAFÇSDFSF" << std::endl;

    //end of pricing

    // add cutting planes
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL( SCIPincludeObjConshdlr(scip, &cuts, TRUE) );

    SCIP_CONS* cons;
    SCIP_CALL( cuts.createCycleCuts(scip, &cons, "CycleRemovalCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE) );
    SCIP_CALL( SCIPaddCons(scip, cons) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons) );

    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );
    SCIP_CALL( SCIPsolve(scip) );

    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT){
        cout << "Reached time limit" << endl;
        return 0;
    } else { // Construct solution
        SCIP_SOL* sol = SCIPgetBestSol(scip);

        for (DNodeIt v(graph); v != INVALID; ++v){
            // get the value of influencing-set variables
            for (unsigned int i = 0; i < infSet[v].size(); i++) {
                double solVal = SCIPgetSolVal(scip, sol, infSet[v][i]);

                // use it to find the amount of incentive paid 
                if (solVal > 0.1) {
                    double cost = 0;//CovModelAllVariables::costInfluencingSet(instance, v, infSets[v][i].elements);

                    std::cout << "lambda[" + instance.nodeName[v] + "," +
                    to_string(cost) + "] = " << solVal << std::endl;
                    solution.incentives[v] = cost;
                } else {
                    solution.incentives[v] = 0;
                }
            }
        }
        std::cout << std::endl;

        for (ArcIt a(graph); a != INVALID; ++a){
            double aux = SCIPgetSolVal(scip, sol, z[a]);

            if(aux > 0.1) {
                DNode u = graph.source(a);
                DNode v = graph.target(a);
                cout << "z[" << instance.nodeName[u] << "," << instance.nodeName[v] << "] = " << aux << endl;
                solution.influence[a] = true;
            } else {
                solution.influence[a] = false;
            }
        }
    }
}
