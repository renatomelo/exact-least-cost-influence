#include "GLCIPBase.h"
#include "generalizedpropagationcons.h"
#include <stack>

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

vector<InfluencingSet> CovModelAllVariables::powerSet(GLCIPInstance &instance, vector<DNode> neighbors, DNode &v)
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
        //add to the list only if is minimal
        double minCost = costInfluencingSet(instance, v, subset.nodes);
        //cout << "minCost = " << minCost << endl;
        int sum = 0;
        for (DNode u : subset.nodes)
        {
            Arc a = findArc(instance.g, u, v);
            assert(a != INVALID);

            sum += instance.influence[a];
        }
        double slack = (sum + minCost) - instance.threshold[v];
        //cout << "slack = " << slack << endl;

        bool isMinimal = TRUE;
        for (DNode u : subset.nodes)
        {
            Arc a = findArc(instance.g, u, v);
            assert(a != INVALID);

            // cout << "influence weight from (" << instance.nodeName[u] << ") = "
            //     << instance.influence[a] << endl;

            if (instance.influence[a] < slack)
            {
                isMinimal = FALSE;
            }
        }

        if (isMinimal)
            pSet.push_back(subset);
        //else  cout << "influencing set NOT minimal\n";
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
        assert(a != INVALID);
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
 * at least one influencing set needs to be selected for each active node
 */
void CovModelAllVariables::addPropagationConstraints(SCIP *scip,
                                                     GLCIPInstance &instance,
                                                     DNodeSCIPVarMap &x,
                                                     DNodeInfSetsMap &infSet)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipCons *cons = new ScipCons(scip, 0, SCIPinfinity(scip));
        //ScipCons *cons = new ScipCons(scip, 0, 0);

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
 * at most one influencing set needs to be selected for each active node
 */
void unicInfSetConstraints(SCIP *scip,
                           GLCIPInstance &instance,
                           DNodeSCIPVarMap &x,
                           DNodeInfSetsMap &infSet)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 1);
        //ScipCons *cons = new ScipCons(scip, 0, 0);

        // summation of all influencing-set varialbes related to a vertex v
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            cons->addVar(infSet[v][i].var, 1);
        }
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
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0);
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
            // trying to force the model hold when using inequality
            /* if (infSet[v][i].nodes.size() == 0)
            {
                cons->addVar(infSet[v][i].var, 1);
            } */
        }

        cons->addVar(z[a], -1);
        cons->commit();
    }
}

bool CovModelAllVariables::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    SCIP *scip = NULL;

    // initialize SCIP enviroment
    SCIP_CALL(SCIPcreate(&scip));

    // explicitely enables the use of debug solution for this scip instance
    SCIPenableDebugSol(scip);

    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);

    // create empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP_with_all_variables", 0, 0, 0, 0, 0, 0, 0));
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 3));

    // to show the branch and bound tree
    SCIP_CALL(SCIPsetStringParam(scip, "visual/vbcfilename", "branchandbound.vbc"));
    SCIP_CALL(SCIPsetBoolParam(scip, "visual/dispsols", TRUE));

    DNodeSCIPVarMap x(instance.g); // active-vertex variables
    ArcSCIPVarMap z(instance.g);   // arc-influence variables

    // create variables x
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        ScipVar *var = new ScipBinVar(scip, "x_" + instance.nodeName[v], 0);
        x[v] = var->var;
    }

    // create variables z
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        DNode u = instance.g.source(a);
        DNode v = instance.g.target(a);
        ScipVar *var = new ScipBinVar(scip, "z_" + instance.nodeName[u] + "," + instance.nodeName[v], 0);
        z[a] = var->var;
    }

    // create influencing-set variables
    DNodeInfSetsMap infSet(instance.g);
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        vector<DNode> neighbors;

        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            neighbors.push_back(instance.g.source(a));
        }
        infSet[v] = powerSet(instance, neighbors, v);

        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            double cost = costInfluencingSet(instance, v, infSet[v][i].nodes);

            string name;
            if (infSet[v][i].nodes.empty())
                name = "Lambda_" + instance.nodeName[v] + "_empty";
            else
            {
                std::stringstream stream;
                for (DNode u : infSet[v][i].nodes)
                    stream << instance.nodeName[u] << ",";
                name = "Lambda_" + instance.nodeName[v] + "_{" + stream.str() + "}";
            }

            ScipVar *var = new ScipContVar(scip, name, 0, SCIPinfinity(scip), cost);
            infSet[v][i].var = var->var;
            infSet[v][i].cost = cost;
        }
    }

    // Propagation constraints:
    addPropagationConstraints(scip, instance, x, infSet);

    //at most one influencing-set can activate each vertex
    //unicInfSetConstraints(scip, instance, x, infSet);

    // add the chosen arcs constraints
    addChosenArcsConstraints(scip, instance, z, infSet);

    // add linking constraints
    addLinkingConstraints(scip, instance, x, z);

    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    //addSmallCycleConstraints(scip, instance, x, z);

    //cout << "\ntrying a brute force" << endl;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
        {
            DNode u = instance.g.target(a);
            if (findArc(instance.g, u, v) != INVALID)
            {
                //adding inequality
                ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0);

                cons->addVar(z[a], 1);
                cons->addVar(z[findArc(instance.g, u, v)], 1);

                cons->addVar(x[v], -1);
                cons->commit();
            }
            else
            {
                for (OutArcIt b(instance.g, u); b != INVALID; ++b)
                {
                    DNode w = instance.g.target(b);
                    if (findArc(instance.g, w, v) != INVALID)
                    {
                        //adding inequality
                        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0);

                        cons->addVar(z[a], 1);
                        cons->addVar(z[b], 1);
                        cons->addVar(z[findArc(instance.g, w, v)], 1);

                        cons->addVar(x[v], -1);
                        cons->commit();
                    }
                    else
                    {
                        for (OutArcIt c(instance.g, w); c != INVALID; ++c)
                        {
                            DNode y = instance.g.target(c);
                            if (findArc(instance.g, y, v) != INVALID)
                            {
                                //adding inequality
                                ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0);

                                cons->addVar(z[a], 1);
                                cons->addVar(z[b], 1);
                                cons->addVar(z[c], 1);
                                cons->addVar(z[findArc(instance.g, y, v)], 1);

                                cons->addVar(x[v], -1);
                                cons->commit();
                            }
                        }
                    }
                }
            }
        }
    }

    SCIP_CALL(SCIPwriteOrigProblem(scip, "glcip_original.lp", "lp", FALSE));

    //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP");
    //exit(0);
    // add generalized propagation constraints
    GeneralizedPropagation *gpc = new GeneralizedPropagation(scip, instance, x, z, infSet);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, gpc, TRUE));

    SCIP_CONS *cons;
    SCIP_CALL(gpc->createGenPropagationCons(scip, &cons, "GPC"));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));
    // end of GPC

    // add cutting planes
    /* CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS *cons1;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons1, "CycleRemovalCuts", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons1));
    SCIP_CALL(SCIPreleaseCons(scip, &cons1)); */
    //end of cutting planes

    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));
    SCIP_CALL(SCIPsolve(scip));

    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT)
    {
        cout << "Reached time limit" << endl;
        return 0;
    }

    //std::cout << SCIPgetSolvingTime(scip) << std::endl;
    // Construct solution
    //CovModel::constructSoltion(scip, instance, solution, z, infSet);

    /* if (CovModel::isFeasible(instance, solution))
        std::cout << "The solution is feasible" << std::endl;
    else
        std::cout << "The solution is NOT feasible" << std::endl */
    ;

    return SCIP_OKAY;
}