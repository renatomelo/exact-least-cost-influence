#include "GLCIPBase.h"

/**
 * Constructs the solution of the GLCIP
 */
/* void constructSoltion(SCIP *scip, GLCIPInstance &instance, GLCIPSolution &solution, ArcSCIPVarMap& z, DNodeInfSetsMap& infSet)
{
    SCIP_SOL *sol = SCIPgetBestSol(scip);

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        solution.incentives[v] = 0.0;
        // get the value of influencing-set variables
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            double solVal = SCIPgetSolVal(scip, sol, infSet[v][i].var);
            // use it to find the amount of incentive paid
            if (solVal > 0.1)
            {
                solution.incentives[v] = infSet[v][i].cost;
                std::cout << "InfluencingSetVar[" + instance.nodeName[v] + "," << solution.incentives[v] << "]  \t= "
                          << solVal << std::endl;
            }
        }
    }
    std::cout << std::endl;

    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        double aux = SCIPgetSolVal(scip, sol, z[a]);

        if (aux > 0.1)
        {
            DNode u = instance.g.source(a);
            DNode v = instance.g.target(a);
            std::cout << "z[" << instance.nodeName[u] << "," << instance.nodeName[v] << "] = " << aux << std::endl;
            solution.influence[a] = true;
        }
        else
        {
            solution.influence[a] = false;
        }
    }
} */

vector<InfluencingSet> CovModelAllVariables::powerSet(vector<DNode> neighbors)
{
    unsigned int pSize = pow(2, neighbors.size());

    vector<InfluencingSet> pSet;

    // run from 000...0 to 111...1
    for (unsigned int i = 0; i < pSize; i++)
    {
        InfluencingSet subset;

        for (unsigned int j = 0; j < neighbors.size(); j++)
        {
            if (i & (1 << j))
            {
                subset.nodes.insert(neighbors[j]);
            }
        }
        pSet.push_back(subset);
    }
    return pSet;
}
/**
 * Computes the cost paid to activate a vertex v with a given set of incoming neigobors
 */
double CovModelAllVariables::costInfluencingSet(GLCIPInstance &instance, DNode v, set<DNode> nodes)
{
    int thr = instance.threshold[v];
    double cost = 0;
    // activation function
    double exertedInfluence = 0;
    for (DNode u : nodes)
    {
        Arc a = findArc(instance.g, u, v);
        //std::cout << instance.nodeName[u] + " ";
        exertedInfluence += instance.influence[a];
    }

    // assuming that the incentives are sorted in an increasing order
    // uses the first incentive that overcomes the threshold of v
    for (unsigned int i = 0; i < instance.incentives[v].size(); i++)
    {
        if (exertedInfluence + instance.incentives[v][i] >= thr)
        {
            cost = instance.incentives[v][i];
            break;
        }
    }
    return cost;
}

/**
 * exactly one influencing set needs to be selected for each active node
 */
void CovModelAllVariables::addPropagationConstraints(SCIP *scip,
                                                     GLCIPInstance &instance,
                                                     DNodeSCIPVarMap &x,
                                                     DNodeInfSetsMap &infSet)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipCons *cons = new ScipCons(scip, 0, 0);

        // summation of all influencing-set varialbes related to a vertex v
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            cons->addVar(infSet[v][i].var, 1);
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
                                                    DNodeInfSetsMap &infSet)
{
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        ScipCons *cons = new ScipCons(scip, 0, 0);
        DNode u = instance.g.source(a);
        DNode v = instance.g.target(a);

        // summation
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            // add variable if u belongs to the influencing set of v
            if (infSet[v][i].nodes.count(u))
            {
                cons->addVar(infSet[v][i].var, 1);
            }
        }

        cons->addVar(z[a], -1);
        cons->commit();
    }
}

bool CovModelAllVariables::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    Digraph &graph = instance.g;
    SCIP *scip = NULL;

    // initialize SCIP enviroment
    SCIP_CALL(SCIPcreate(&scip));

    // explicitely enables the use of debug solution for this scip instance
    SCIPenableDebugSol(scip);

    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));

    // create empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP ColGeneration", 0, 0, 0, 0, 0, 0, 0));

    DNodeSCIPVarMap x(graph); // active-vertex variables
    ArcSCIPVarMap z(graph);   // arc-influence variables
    //DNodeSCIPVarsMap infSet(graph); // influencing-set variables

    // create variables x
    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        ScipVar *var = new ScipBinVar(scip, "x_" + instance.nodeName[v], 0);
        x[v] = var->var;
    }

    // create variables z
    for (ArcIt a(graph); a != INVALID; ++a)
    {
        DNode u = graph.source(a);
        DNode v = graph.target(a);
        ScipVar *var = new ScipBinVar(scip, "z_" + instance.nodeName[u] + "" + instance.nodeName[v], 0);
        z[a] = var->var;
    }

    // create influencing-set variables
    DNodeInfSetsMap infSet(graph);
    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        vector<DNode> neighbors;

        for (InArcIt a(graph, v); a != INVALID; ++a)
        {
            neighbors.push_back(graph.source(a));
        }
        infSet[v] = powerSet(neighbors);

        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            double cost = costInfluencingSet(instance, v, infSet[v][i].nodes);
            ScipVar *var = new ScipContVar(scip, "infSet_" + instance.nodeName[v] + "" + to_string(i), 0, SCIPinfinity(scip), cost);
            //infSet[v].push_back(var->var);
            infSet[v][i].var = var->var;
            infSet[v][i].cost = cost;
        }
    }

    // Propagation constraints:
    addPropagationConstraints(scip, instance, x, infSet);

    // add the chosen arcs constraints
    addChosenArcsConstraints(scip, instance, z, infSet);

    // add linking constraints
    addLinkingConstraints(scip, instance, x, z);

    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    addSmallCycleConstraints(scip, instance, x, z);

    // add cutting planes
    //ArcModel::addCuttingPlanes(scip, instance, x, z);
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS *cons;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons, "CycleRemovalCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));

    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));
    SCIP_CALL(SCIPsolve(scip));

    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT)
    {
        cout << "Reached time limit" << endl;
        return 0;
    }
    // Construct solution
    CovModel::constructSoltion(scip, instance, solution, z, infSet);

    if (CovModel::isFeasible(instance, solution))
        std::cout << "The solution is feasible" << std::endl;
    else
        std::cout << "The solution is NOT feasible" << std::endl;

    return SCIP_OKAY;
}