#include "GLCIPBase.h"
#include <map>
#include "pricer_glcip.h"
#include "branch_glcip.h"
#include "event_glcip.h"
#include "consarcmarker.h"

using namespace arcmarker;

/**
 * Add to the model the influencing-set variables considering no incoming
 *  arcs, that is, paying the integral incentive for every vertex
 */
SCIP_RETCODE incentivesForAll(SCIP *scip,
                              GLCIPInstance &instance,
                              DNodeConsMap &vertCons,
                              DNodeInfSetsMap &infSet)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        double cost = GLCIPBase::cheapestIncentive(instance, v, 0);

        SCIP_VAR *var;
        std::string name = "infSetVar" + instance.nodeName[v] + "_empty";
        SCIP_CALL(SCIPcreateVar(scip, &var,
                                name.c_str(),            // var name
                                0,                       // lower bound
                                SCIPinfinity(scip),      // upper bound
                                cost,                    // coeficient in the objective function
                                SCIP_VARTYPE_CONTINUOUS, // continuous variable
                                FALSE,                   // initial variable
                                FALSE,                   // removable variable
                                NULL, NULL, NULL, NULL, NULL));
        // add new variable to the list of variables to price into LP
        SCIP_CALL(SCIPaddVar(scip, var));

        // add to each vertex constraint
        SCIP_CALL(SCIPaddCoefLinear(scip, vertCons[v], var, 1));

        // data structure to save the variables and associated costs
        InfluencingSet initial;
        initial.var = var;
        initial.cost = cost;
        infSet[v].push_back(initial);

        /*  std::cout << "Adding variable "
                  << "infSetVar" + instance.nodeName[v] + "empty" << std::endl; */
        SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    return SCIP_OKAY;
}

void showActivatedNodes(GLCIPInstance &instance, set<DNode> actives, DNodeInfSetsMap &infSet)
{
    std::cout << "Active vertex and its influencing sets" << std::endl;
    for (DNode v : actives)
    {
        std::cout << instance.nodeName[v] + ": ";
        for (unsigned int i = 0; i < infSet[v].size(); i++)
        {
            for (DNode u : infSet[v][i].nodes)
            {
                std::cout << " " + instance.nodeName[u];
            }
            std::cout << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << " There are " << actives.size() << " active nodes" << std::endl;
}

/**
 * The node with minimal incentive to activate it is chosen
 */
double getMinIncentiveNode(GLCIPInstance &instance, set<DNode> actives, DNode &node)
{
    double minCost = 1e+20;
    //DNode choosed;
    set<DNode> activeNeigbors;

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        activeNeigbors.clear();

        if (!actives.count(v))
        {
            //std::cout << instance.nodeName[v] << " is inactive " << std::endl;
            for (InArcIt a(instance.g, v); a != INVALID; ++a)
            {
                DNode u = instance.g.source(a);

                //std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v]
                //          << " visited now " << std::endl;

                if (actives.count(u))
                {
                    //std::cout << "Node " << instance.nodeName[u] << " is an active neigbor of "
                    //          << instance.nodeName[v] << std::endl;
                    activeNeigbors.insert(u);
                }
            }

            double cost = GLCIPBase::costInfluencingSet(instance, v, activeNeigbors);
            if (cost < minCost)
            {
                //std::cout << "Cost of " + instance.nodeName[v]
                //          << " is smaller than " << minCost << std::endl;

                minCost = cost;
                node = v;

                // if the cost is zero it cannot be improved, then stop
                if (cost == 0)
                    break;
            }
        }
    }

    //std::cout << "Node " + instance.nodeName[node] + " has the minimum incentive: "
            //  << minCost << std::endl;

    return minCost;
}

/** 
* Greedy construction heuristic to obtain feasible solutions to warm-start column 
* generation with an initial set of influensing-set variables:
* 
* At each iteration we activate a not yet active node by paying the minimum available 
* incentive to reach its hurdle, taking into account the current influence coming 
* from already active neighbors.
*/
set<DNode> greedyConstruction(GLCIPInstance &instance, DNodeInfSetsMap &infSet)
{
    set<DNode> actives; // start with a empty solution
    while (actives.size() < instance.alpha * instance.n)
    {
        DNode v;
        double minCost = getMinIncentiveNode(instance, actives, v);
        // save v's influencing-set and activate it
        InfluencingSet ifs;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            DNode u = instance.g.source(a);
            if (actives.count(u))
            {
                ifs.nodes.insert(u);
            }
        }
        ifs.cost = minCost; //GLCIPBase::costInfluencingSet(instance, v, ifs.nodes);
        infSet[v].push_back(ifs);
        actives.insert(v);

        //std::cout << instance.nodeName[v] << " was activated " << std::endl;
    }

    //showActivatedNodes(instance, actives, infSet);
    return actives;
}
/**
 * Add to the model the decision variables for an initial solution obtained from 
 * the greedy heuristic construction method
 */
SCIP_RETCODE addHeurInitialSol(SCIP *scip,
                               GLCIPInstance &instance,
                               ArcConsMap &arcCons,
                               DNodeConsMap &vertCons,
                               DNodeInfSetsMap &infSet)
{
    set<DNode> activated = greedyConstruction(instance, infSet);

    //create a var for each vertex v and add it to the associated vertex constraint
    //for (DNodeIt v(instance.g); v != INVALID; ++v)
    for (DNode v: activated)
    {
        // give a representative name to the variable
        std::string name;
        if (infSet[v][0].nodes.size() > 0)
        {
            // give a significative name for the variable
            std::stringstream stream;
            for (DNode u : infSet[v][0].nodes)
                stream << instance.nodeName[u];
            name = "infSetVar_" + instance.nodeName[v] + "_" + stream.str();
        }
        else
            name = "infSetVar_" + instance.nodeName[v] + "_empty";

        //double cost = GLCIPBase::costInfluencingSet(instance, v, infSet[v][0].nodes);
        SCIP_VAR *var;
        SCIP_CALL(SCIPcreateVar(scip, &var,
                                name.c_str(),            // var name
                                0,                       // lower bound
                                SCIPinfinity(scip),      // upper bound
                                infSet[v][0].cost,       // coeficient in the objective function
                                SCIP_VARTYPE_CONTINUOUS, // continuous variable
                                FALSE,                   // initial variable
                                FALSE,                   // removable variable
                                NULL, NULL, NULL, NULL, NULL));
        // add new variable to the list of variables to price into LP
        SCIP_CALL(SCIPaddVar(scip, var));

        // add to each vertex constraint
        SCIP_CALL(SCIPaddCoefLinear(scip, vertCons[v], var, 1));

        // until this momem each vertex has exactly one assossiated influencing-set
        // add the vars associated with each arc constraint
        for (DNode u : infSet[v][0].nodes)
        {
            //std::cout << instance.nodeName[u] << " ";
            Arc a = findArc(instance.g, u, v);
            assert(a != INVALID);
            SCIP_CALL(SCIPaddCoefLinear(scip, arcCons[a], var, -1.0));
        }

        infSet[v][0].var = var;

        //std::cout << "Adding variable " << name << std::endl;
        SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    //std::cout << "Cost of this solution: " << totalCost << std::endl;
    return SCIP_OKAY;
}

/**
 * Performs a propagation from a set of incentives to verify wheter the instance is 
 * feasible
 */
bool CovModel::isFeasible(GLCIPInstance &instance, GLCIPSolution &solution)
{
    // start wiht empty solution
    Digraph::NodeMap<set<DNode>> influencers(instance.g);
    set<DNode> actives;
    list<DNode> seeds;

    // initialize set with seed nodes
    //std::cout << "Initial adopters:  ";
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        if (solution.incentives[v] >= instance.threshold[v])
        {
            seeds.push_back(v);
            //std::cout << instance.nodeName[v] << " ";
        }
    }
    //std::cout << std::endl;

    // while the seed set is not empty try to activate non active vertices
    while (seeds.size() > 0)
    {
        DNode u = seeds.front();
        seeds.pop_front();
        actives.insert(u);

        for (OutArcIt a(instance.g, u); a != INVALID; ++a)
        {
            DNode v = instance.g.target(a);
            if (!actives.count(v) && solution.influence[a])
            {
                //std::cout << instance.nodeName[v] << " is inactive " << std::endl;
                //std::cout << "arc_" << instance.nodeName[u] << "_"
                //<< instance.nodeName[v] << " is active " << std::endl;
                influencers[v].insert(u);

                double exerterdInfluence = 0;
                for (DNode w : influencers[v])
                {
                    Arc e = findArc(instance.g, w, v);
                    assert(e != INVALID);

                    exerterdInfluence += instance.influence[e];
                }
                //std::cout << " exerterd influence + incentive: "
                //<< (exerterdInfluence + solution.incentives[v]) << std::endl;
                if (exerterdInfluence + solution.incentives[v] >= instance.threshold[v])
                {
                    //std::cout << instance.nodeName[v] << " is inserted in seed set "
                    //<< std::endl;
                    seeds.push_back(v);
                }
            }
        }
    }

    if (actives.size() >= instance.alpha * instance.n)
        return true;

    return false;
}
/**
 * Constructs the solution of the GLCIP
 */
void CovModel::constructSoltion(SCIP *scip,
                                GLCIPInstance &instance,
                                GLCIPSolution &solution,
                                ArcSCIPVarMap &z,
                                DNodeInfSetsMap &infSet)
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
                std::cout << SCIPvarGetName(infSet[v][i].var) << " \t\t= " << solVal << std::endl;
            }
        }
    }
    std::cout << std::endl;

    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        double aux = SCIPgetSolVal(scip, sol, z[a]);

        solution.influence[a] = false;
        if (aux > 0.1)
        {
            DNode u = instance.g.source(a);
            DNode v = instance.g.target(a);
            std::cout << "z[" << instance.nodeName[u] << "," << instance.nodeName[v] << "] = " << aux << std::endl;
            solution.influence[a] = true;
        }
    }
}

bool CovModel::run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit)
{
    Digraph &graph = instance.g;
    SCIP *scip = NULL;

    // initialize SCIP enviroment
    SCIP_CALL(SCIPcreate(&scip));

    // explicitely enables the use of debug solution for this scip instance
    SCIPenableDebugSol(scip);

    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    /*   SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetRealParam(scip, "separating/minefficacy", 0.001)); */

    // create empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP_Column_Generation", 0, 0, 0, 0, 0, 0, 0));

    SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 3));

    // to show the branch and bound tree
    SCIP_CALL(SCIPsetStringParam(scip, "visual/vbcfilename", "branchandbound.vbc"));
    SCIP_CALL(SCIPsetBoolParam(scip, "visual/dispsols", TRUE));
    SCIP_CALL(SCIPsetBoolParam(scip, "visual/realtime", FALSE));

    DNodeSCIPVarMap x(graph);      // active-vertex variables
    ArcSCIPVarMap z(graph);        // arc-influence variables
    DNodeInfSetsMap infSet(graph); // influencing-set variables (used to show the solution)

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
        std::string name = "z_" + instance.nodeName[u] + "" + instance.nodeName[v];
        ScipVar *var = new ScipBinVar(scip, name, 0);
        z[a] = var->var;
    }

    // add vertex-coverage constraints
    DNodeConsMap vertCons(graph);
    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        ScipConsPrice *cons = new ScipConsPrice(scip, 0, SCIPinfinity(scip));
        cons->addVar(x[v], -1);
        vertCons[v] = cons->cons;
        cons->commit();
    }

    // add arc-coverage constraints
    ArcConsMap arcCons(graph);
    for (ArcIt a(graph); a != INVALID; ++a)
    {
        ScipConsPrice *cons = new ScipConsPrice(scip, 0, SCIPinfinity(scip));
        cons->addVar(z[a], 1);
        arcCons[a] = cons->cons;
        cons->commit();
    }
    //Pricing to generate the propagation constraints and the chosen arcs constraints

    // start with a initial solution obtained heuristically
    //incentivesForAll(scip, instance, vertCons, infSet); // construct an initial solution
    //add initial heuristic solution
    SCIP_CALL(addHeurInitialSol(scip, instance, arcCons, vertCons, infSet));
    //exit(0);

    // add linking constraints
    addLinkingConstraints(scip, instance, x, z);

    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    addSmallCycleConstraints(scip, instance, x, z);

    //SCIPwriteOrigProblem(scip, "initial.lp", "lp", FALSE);

    // include pricer
    ArcIntMap isOnSolution(graph);
    static const char *PRICER_NAME = "GLCIP_pricer";
    ObjPricerGLCIP *pricer = new ObjPricerGLCIP(scip, PRICER_NAME, instance, z, x, arcCons, vertCons, infSet, isOnSolution);

    SCIP_CALL(SCIPincludeObjPricer(scip, pricer, TRUE));
    SCIP_CALL(SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME)));
    //end of pricing

   /*  ConshdlrArcMarker *arcMarker = new ConshdlrArcMarker(scip, instance, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, arcMarker, TRUE));

    //include branching rule
    static const char* BRANCH_NAME = "GLCIP_branch";
    ObjBranchruleGLCIP* branch = new ObjBranchruleGLCIP(scip, BRANCH_NAME, instance, x, z, infSet);

    SCIP_CALL(SCIPincludeObjBranchrule(scip, branch, TRUE)); */
    //end of branching rule

    // include event handler pluging
    /* static const char *EVENTHDLR_NAME = "GLCIP_eventhdlr";
    ObjEventhdlrGLCIP *event = new ObjEventhdlrGLCIP(scip, EVENTHDLR_NAME, instance);
    SCIP_CALL(SCIPincludeObjEventhdlr(scip, event, TRUE)); */

    // add cutting planes
    CycleCutsGenerator cuts = CycleCutsGenerator(scip, instance, x, z);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, &cuts, TRUE));

    SCIP_CONS *cons;
    SCIP_CALL(cuts.createCycleCuts(scip, &cons, "CycleRemovalCuts",
                                   FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE));
    SCIP_CALL(SCIPaddCons(scip, cons));
    SCIP_CALL(SCIPreleaseCons(scip, &cons));

    SCIP_CALL(SCIPsetRealParam(scip, "limits/time", timeLimit));
    SCIP_CALL(SCIPsolve(scip));

    if (SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT)
    {
        cout << "Reached time limit" << endl;
        return 0;
    }

    std::cout << SCIPgetSolvingTime(scip) << std::endl;

    // Construct solution
    constructSoltion(scip, instance, solution, z, infSet);

    if (isFeasible(instance, solution))
        std::cout << "The solution is feasible" << std::endl;
    else
        std::cout << "The solution is NOT feasible" << std::endl;

    return SCIP_OKAY;
}
