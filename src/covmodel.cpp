#include "GLCIPBase.h"
#include <map>
#include "pricer_glcip.h"

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
        double cost = 0;
        // uses the first incentive that overcomes the threshold of v
        for (unsigned int i = 0; i < instance.incentives[v].size(); i++)
        {
            if (instance.incentives[v][i] >= instance.threshold[v])
            {
                cost = instance.incentives[v][i];
                break;
            }
        }

        SCIP_VAR *var;
        std::string name = "infSetVar" + instance.nodeName[v] + "_empty";
        SCIP_CALL(SCIPcreateVar(scip, &var,
                                name.c_str(),            // var name
                                0,                       // lower bound
                                SCIPinfinity(scip),      // upper bound
                                cost,                    // coeficient in the objective function
                                SCIP_VARTYPE_CONTINUOUS, // continuous variable
                                TRUE,                    // initial variable
                                FALSE,                    // removable variable
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
        //SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    return SCIP_OKAY;
}

/** 
* greedy construction heuristic to obtain feasible solutions to warm-start column 
* generation with an initial set of influensing-set variables
*/
/* void greedyConstruction(SCIP *scip, GLCIPInstance &instance, ArcConsMap &arcCons, DNodeConsMap &vertCons)
{
   vector<DNode> activeNodes;
   vector<DNode> seedNodes;
   //double paidIncentives[] = new double[instance.n];

   // initialize set with seed nodes - add every node which has no incoming arcs by paying its incentives
   for (DNodeIt v(instance.g); v != INVALID; ++v)
   {
      InDegMap<Digraph> inDeg(instance.g);
      if (inDeg[v] == 0)
      {
         seedNodes.push_back(v);
      }
   }

   while (seedNodes.size() > 0)
   {
      DNode u = seedNodes.back();
   }
} */

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
                //std::cout << "arc_" << instance.nodeName[u] << "_" << instance.nodeName[v] << " is active " << std::endl;
                influencers[v].insert(u);

                double exerterdInfluence = 0;
                for (DNode w : influencers[v])
                {
                    Arc e = findArc(instance.g, w, v);
                    assert(e != INVALID);

                    exerterdInfluence += instance.influence[e];
                }
                //std::cout << " exerterd influence + incentive: " << (exerterdInfluence + solution.incentives[v]) << std::endl;
                if (exerterdInfluence + solution.incentives[v] >= instance.threshold[v])
                {
                    //std::cout << instance.nodeName[v] << " is inserted in seed set " << std::endl;
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
void CovModel::constructSoltion(SCIP *scip, GLCIPInstance &instance, GLCIPSolution &solution, ArcSCIPVarMap &z, DNodeInfSetsMap &infSet)
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
                std::cout << "InfluencingSetVar[" + instance.nodeName[v] + "," 
                          << solution.incentives[v] << "]  \t= " << solVal << std::endl;
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
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE));
    SCIP_CALL(SCIPsetRealParam(scip, "separating/minefficacy", 0.001));

    // create empty problem
    SCIP_CALL(SCIPcreateProb(scip, "GLCIP_Column_Generation", 0, 0, 0, 0, 0, 0, 0));

    //SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);
    /* SCIP_CALL(SCIPsetIntParam(scip, "display/verblevel", 5));
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
     */
    
    DNodeSCIPVarMap x(graph); // active-vertex variables
    ArcSCIPVarMap z(graph);   // arc-influence variables
    DNodeInfSetsMap infSet(graph);  // influencing-set variables (used to show the solution)

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
    incentivesForAll(scip, instance, vertCons, infSet); // construct an initial solution

    // include pricer
    static const char *PRICER_NAME = "GLCIP_pricer";
    ObjPricerGLCIP *pricer = new ObjPricerGLCIP(scip, PRICER_NAME, instance, z, x, arcCons, vertCons, infSet);

    SCIP_CALL(SCIPincludeObjPricer(scip, pricer, TRUE));
    SCIP_CALL(SCIPactivatePricer(scip, SCIPfindPricer(scip, PRICER_NAME)));
    //end of pricing

     // add linking constraints
    addLinkingConstraints(scip, instance, x, z);

    // add coverage constraints:
    addCoverageConstraints(scip, instance, x);

    // add all cycles of size up to 4
    addSmallCycleConstraints(scip, instance, x, z);

    // add cutting planes
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
    constructSoltion(scip, instance, solution, z, infSet);

    if (isFeasible(instance, solution))
        std::cout << "The solution is feasible" << std::endl;
    else
        std::cout << "The solution is NOT feasible" << std::endl;

    return SCIP_OKAY;
}
