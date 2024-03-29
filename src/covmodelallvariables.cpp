#include "GLCIPBase.h"
#include "generalizedpropagationcons.h"
#include <stack>
#include "heur_mininfluence.h"

vector<InfluencingSet> CovModelAllVariables::powerSet(
    GLCIPInstance &instance, 
    vector<DNode> neighbors, 
    DNode &v)
{
    unsigned int pSize = pow(2, neighbors.size());

    vector<InfluencingSet> pSet;

    // run from 000...0 to 111...1
    for (unsigned int i = 0; i < pSize; i++)
    {
        InfluencingSet subset(instance, v);

        for (unsigned int j = 0; j < neighbors.size(); j++)
        {
            if (i & (1 << j))
            {
                subset.addNode(neighbors[j]);
            }
        }
        //add to the list only if is minimal
        double minCost = costInfluencingSet(instance, v, subset.getNodes());
        //cout << "minCost = " << minCost << endl;
        int sum = 0;
        for (DNode u : subset.getNodes())
        {
            Arc a = findArc(instance.g, u, v);
            assert(a != INVALID);

            sum += instance.influence[a];
        }
        double slack = (sum + minCost) - instance.threshold[v];
        //cout << "slack = " << slack << endl;

        bool isMinimal = TRUE;
        for (DNode u : subset.getNodes())
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
void CovModelAllVariables::addPropagationConstraints(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    DNodeInfSetsMap &infSet)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //ScipCons *cons = new ScipCons(scip, 0, SCIPinfinity(scip), "vertex-coverage cons");
        ScipCons *cons = new ScipCons(scip, 0, 0, "vertex-coverage cons");

        // summation of all influencing-set varialbes related to a vertex v
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            cons->addVar(infSet[v][i].getVar(), 1);
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
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 1, "unic-inf-set cons");
        //ScipCons *cons = new ScipCons(scip, 0, 0);

        // summation of all influencing-set varialbes related to a vertex v
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            cons->addVar(infSet[v][i].getVar(), 1);
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
        //ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0, "arc-coverage cons");
        ScipCons *cons = new ScipCons(scip, 0, 0, "arc-coverage cons");
        DNode u = instance.g.source(a);
        DNode v = instance.g.target(a);

        // summation
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            // add variable if u belongs to the influencing set of v
            if (infSet[v][i].getNodes().count(u))
            {
                cons->addVar(infSet[v][i].getVar(), 1);
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
    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 1));

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
            /* double cost = costInfluencingSet(instance, v, infSet[v][i].getNodes());

            string name;
            if (infSet[v][i].getNodes().empty())
                name = "Lambda_" + instance.nodeName[v] + "_empty";
            else
            {
                std::stringstream stream;

                const char *separator = "";
                for (DNode u : infSet[v][i].getNodes())
                {
                    stream << separator << instance.nodeName[u];
                    separator = ",";
                }

                name = "Lambda_" + instance.nodeName[v] + "_{" + stream.str() + "}";
            } */

            infSet[v][i].setCost(costInfluencingSet(instance, v, infSet[v][i].getNodes()));

            ScipVar *var = new ScipContVar(scip, infSet[v][i].getName(), 0, SCIPinfinity(scip), infSet[v][i].getCost());
            /* infSet[v][i].var = var->var;
            infSet[v][i].cost = cost; */
            infSet[v][i].setVar(var->var);
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
    addAllSmallDirectedCycles(scip, instance, x, z);
    //addSmallCycleConstraints(scip, instance, x, z);

    //SCIP_CALL(SCIPwriteOrigProblem(scip, "glcip_original.lp", "lp", FALSE));

    //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP");
    //exit(0);
    // add generalized propagation constraints
    vector<Phi> gpcrows;
    GeneralizedPropagation *gpc = new GeneralizedPropagation(scip, instance, x, z, infSet, gpcrows);
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

    //SCIP_CALL( SCIPincludeObjHeur(scip, new HeurMinInfluence(scip, instance, x, z, infSet), TRUE) );

    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));
    SCIP_CALL(SCIPsolve(scip));

    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT)
    {
        //cout << "Reached time limit" << endl;
        printf("%.2lf\t%lld\t%d\t%lf\t%lf\t%.2lf\n", SCIPgetSolvingTime(scip), 
                                             SCIPgetNNodes(scip),
                                             SCIPgetNContVars(scip), 
                                             SCIPgetDualbound(scip), 
                                             SCIPgetPrimalbound(scip),
                                             SCIPgetGap(scip));

        return 0;
    }

    printf("%.2lf\t%lld\t%d\t%lf\t%lf\t%.2lf\n", SCIPgetSolvingTime(scip), 
                                             SCIPgetNNodes(scip),
                                             SCIPgetNContVars(scip),
                                             SCIPgetDualbound(scip), 
                                             SCIPgetPrimalbound(scip),
                                             SCIPgetGap(scip));

    //std::cout << SCIPgetSolvingTime(scip) << std::endl;
    // Construct solution
    /* CovModel::constructSoltion(scip, instance, solution, z, infSet);

    if (CovModel::isFeasible(instance, solution))
        std::cout << "The solution is feasible" << std::endl;
    else
        std::cout << "The solution is NOT feasible" << std::endl; */

    return SCIP_OKAY;
}