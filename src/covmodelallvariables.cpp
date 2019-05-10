#include "GLCIPBase.h"

vector<influencingSet> CovModelAllVariables::powerSet(vector<DNode> neighbors){
    unsigned int pSize = pow(2, neighbors.size());

    vector<influencingSet> pSet;

    // run from 000...0 to 111...1
    for (unsigned int i = 0; i < pSize; i++) {
        influencingSet subset;

        for (unsigned int j = 0; j < neighbors.size(); j++) {
            if (i & (1 << j)) {
                subset.elements.insert(neighbors[j]);
            }
        }
        pSet.push_back(subset);
    }
    return pSet;
}
/**
 * Computes the influence of a given set of incoming neigobors of a vertex v
 */
double CovModelAllVariables::costInfluencingSet(GLCIPInstance &instance, DNode v, set<DNode> elements){
    int thr = instance.threshold[v];

    // activation function
    double exertedInfluence = 0;
    for (DNode u: elements){
        Arc a = findArc(instance.g, u, v);
        //std::cout << instance.nodeName[u] + " ";
        exertedInfluence += instance.influence[a];
    }
   // std::cout << ": excerts inflence of " + to_string(exertedInfluence);

    // assuming that the incentives are sorted in an increasing order
    // uses the first incentive that overcomes the threshold of v
    for (int i = 0; i < instance.incentives[v].size(); i++){
        if (exertedInfluence + instance.incentives[v][i] >= thr){
            //std::cout << " with incentive of " + to_string(instance.incentives[v][i]);
            //std::cout << " and threshold of " + to_string(thr) << std::endl;
            return instance.incentives[v][i];
        }
    }
}

/**
 * exactly one influencing set needs to be selected for each active node
 */
void CovModelAllVariables::addPropagationConstraints(SCIP *scip,
                                                     GLCIPInstance &instance, 
                                                     DNodeSCIPVarMap &x, 
                                                     DNodeSCIPVarsMap &infSet, 
                                                     DNodeInfSetsMap &infSets) {
    for (DNodeIt v(instance.g); v != INVALID; ++v){
        ScipCons* cons = new ScipCons(scip, 0, 0);

        // summation of all influencing-set varialbes related to a vertex v
        for (int i = 0; i < infSet[v].size(); i++) {
            cons->addVar(infSet[v][i], 1);
        }
        cons->addVar(x[v], -1);
        cons->commit();
    }
}

/**
 * all arcs on which influence is exerted are chosen as well
 */
void CovModelAllVariables::addChosenArcsConstraints(SCIP *scip,
                                                    GLCIPInstance &instance, 
                                                    ArcSCIPVarMap &z, 
                                                    DNodeSCIPVarsMap &infSet, 
                                                    DNodeInfSetsMap &infSets) {
    for (ArcIt a(instance.g); a != INVALID; ++a) {
        ScipCons* cons = new ScipCons(scip, 0, 0);
        DNode u = instance.g.source(a);
        DNode v = instance.g.target(a);

        // summation
        for (int i = 0; i < infSet[v].size(); i++){
            // add variable if u belongs to the influencing set of v
            if (infSets[v][i].elements.count(u)){
               //std::cout << instance.nodeName[u] + " is in the " + to_string(i) +
               //             "th influencing set of " + instance.nodeName[v] << std::endl;
                cons->addVar(infSet[v][i], 1);
            }
        }

        cons->addVar(z[a], -1);
        cons->commit();
    }
}

bool CovModelAllVariables::run (GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    Digraph &graph = instance.g;
    SCIP* scip = NULL;

    // initialize SCIP enviroment
    SCIP_CALL( SCIPcreate(&scip) );

    // explicitely enables the use of debug solution for this scip instance
    SCIPenableDebugSol(scip);
    
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));

    // create empty problem
    SCIP_CALL( SCIPcreateProb(scip, "GLCIP ColGeneration", 0, 0, 0, 0, 0, 0, 0) );

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

    // create influencing-set variables
    DNodeInfSetsMap infSets(graph); // stores the associated sets
    for (DNodeIt v(graph); v != INVALID; ++v){
        vector<DNode> neighbors;

        for (InArcIt a(graph, v); a != INVALID; ++a) {
            neighbors.push_back(graph.source(a));
        }
        infSets[v] = powerSet(neighbors);

        /*for (int i = 0; i < infSets[v].size(); i++) {
            for(DNode v : infSets[v][i].elements){
                std::cout << instance.nodeName[v] + " ";
            }           
            std::cout << std::endl;
        }*/

        for (int i = 0; i < infSets[v].size(); i++){
            double cost = costInfluencingSet(instance, v, infSets[v][i].elements);
            ScipVar* var = new ScipBinVar(scip, "infSet_" + instance.nodeName[v] 
                                                + "" + to_string(i), cost);
            infSet[v].push_back(var->var);
        }
    }

    // Propagation constraints: 
    addPropagationConstraints(scip, instance, x, infSet, infSets);

    // add the chosen arcs constraints
    addChosenArcsConstraints(scip, instance, z, infSet, infSets);

    // add linking constraints 
    addLinkingConstraints(scip, instance, x, z);
    
    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    addSmallCycleConstraints(scip, instance, x, z);

    // add cutting planes
    //ArcModel::addCuttingPlanes(scip, instance, x, z);
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
                    double cost = costInfluencingSet(instance, v, infSets[v][i].elements);

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