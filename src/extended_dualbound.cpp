#include "extended_dualbound.h"

ExtendedDualBound::ExtendedDualBound(
    SCIP *scip,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z
    //DNodeSCIPVarsMap &p_xip
    ) : ObjRelax(scip,
                 "extended-dual-bound",
                 "Extended dual bound for LCIP",
                 -1.0,   //priority of the relaxator (negative: after LP, non-negative: before LP)
                 1,      //frequency for calling relaxator
                 FALSE), //Does the relaxator contain all cuts in the LP?
        instance(p_instance),
        x(p_x),
        z(p_z),
        //xip(p_xip),
        sol_(NULL)
{
    env = new GRBEnv();
}

SCIP_DECL_RELAXFREE(ExtendedDualBound::scip_free)
{
    //cout << "RELAXFREE()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINIT(ExtendedDualBound::scip_init)
{
    //cout << "RELAXEINIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXIT(ExtendedDualBound::scip_exit)
{
    //cout << "RELAXEXIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINITSOL(ExtendedDualBound::scip_initsol)
{
    //cout << "RELAXINITSOL()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXITSOL(ExtendedDualBound::scip_exitsol)
{
    //cout << "RELAXEXITSOL()" << endl;
    return SCIP_OKAY;
}

void ExtendedDualBound::getCondensedGraph(
    Digraph &condensed,
    Digraph &graph,
    DNodeDNodeMap &nodeRef,
    Digraph::NodeMap<int> &components,
    vector<vector<DNode>> &listOfComponents)
{
    //Digraph condensed;
    int nComponents = countStronglyConnectedComponents(graph);
    for (int i = 0; i < nComponents; i++)
    {
        condensed.addNode();
        //DNode v = condensed.addNode();
        //assert(condensed.id(v) == i);
    }

    //compute the components of 'graph'
    stronglyConnectedComponents(graph, components);

    //save each component in a vector of vertices
    for (DNodeIt v(graph); v != INVALID; ++v)
        listOfComponents[components[v]].push_back(nodeRef[v]);
}

void ExtendedDualBound::getCondensedThresholds(
    Digraph &graph,
    vector<vector<DNode>> &listOfComponents,
    int nComponents,
    DNodeValueMap &thr)
{
    //vector<double> thr(nComponents);

    //find the minimum threshold on each component
    for (int i = 0; i < nComponents; i++)
    {
        double minThreshold = 1e+20;
        //printf("size of component %d: %ld\n", i, listOfComponents[i].size());
        for (size_t j = 0; j < listOfComponents[i].size(); j++)
        {
            DNode w = listOfComponents[i][j];
            minThreshold = min(minThreshold, instance.threshold[w]);
            //printf("threshold of %s: %f\n", instance.nodeName[w].c_str(), instance.threshold[w]);
        }
        thr[graph.nodeFromId(i)] = minThreshold;
        //printf("threshold(%d) = %.1f\n", i, minThreshold);
    }
    //return thr;
}

// add the condensed arcs in the condensed graph
void ExtendedDualBound::getCondensedArcWeights(
    Digraph &condensed,
    Digraph &graph,
    ArcArcMap &arcRef,
    ArcValueMap &weights,
    Digraph::NodeMap<int> &components)
{
    //get the cut arcs of the strongly connected components
    Digraph::ArcMap<bool> cutArcs(graph, FALSE);
    stronglyConnectedCutArcs(graph, cutArcs);

    //ArcValueMap condensedInfluence(condensed);
    for (ArcIt a(graph); a != INVALID; ++a)
    {
        //compute the weigh of influence of the arcs
        //each arc receives the total of weights of all arcs from a component to another
        if (cutArcs[a])
        {
            //reference to what component each vertex belong
            int i = components[graph.source(a)];
            int j = components[graph.target(a)];

            if (i == j)
            {
                cout << "something wrong i == j\n";
            }

            Arc b = findArc(condensed, condensed.nodeFromId(i), condensed.nodeFromId(j));
            if (b == INVALID)
            {
                Arc c = condensed.addArc(condensed.nodeFromId(i), condensed.nodeFromId(j));
                weights[c] = instance.influence[arcRef[a]];
            }
            else
                weights[b] += instance.influence[arcRef[a]];
        }
        else
        {
            //reference to what component each vertex belong
            int i = components[graph.source(a)];
            int j = components[graph.target(a)];

            if (i != j)
            {
                cout << "something wrong i != j\n";
            }
        }
    }
}

void ExtendedDualBound::printCondensedArcs(Digraph &condensed, ArcValueMap &condensedInfluence)
{
    for (ArcIt a(condensed); a != INVALID; ++a)
    {
        int i = condensed.id(condensed.source(a));
        int j = condensed.id(condensed.target(a));

        printf("condensed arc: %d - > %d: %f", i, j, condensedInfluence[a]);
    }
}

double ExtendedDualBound::getMinimumThreshold(DNode &node)
{
    double minimum = 1e+09;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        //cout << "threshold of " << instance.nodeName[v] << ": " << instance.threshold[v] << endl;
        double thr = instance.threshold[v];
        if (thr < minimum)
        {
            minimum = thr;
            node = v;

            // stop if thr(v) = 1 because there is no smaller threshold
            if (instance.threshold[v] == 1)
                break;
        }
    }

    return minimum;
}

/* int ExtendedDualBound::getIndexOfChepeastIncentive(GLCIPInstance &instance, DNode &node)
{
    int index = 0;
    for (size_t i = 0; i < instance.incentives[node].size(); i++)
    {
        if (instance.incentives[node][i] >= instance.threshold[node])
        {
            //cout << "incentive paid: " << instance.incentives[node][i] << endl;
            index = i;
            break;
        }
    }
    return index;
} */

void ExtendedDualBound::getSubGraph(
    SCIP *scip,
    Digraph &graph,
    DNodeDNodeMap &nodeRef,
    ArcArcMap &arcRef)
{
    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    for (ArcIt a(graph); a != INVALID; ++a)
    {

        if (SCIPisEQ(scip, SCIPvarGetUbLocal(z[arcRef[a]]), 0))
        {
            //removing arc variables fixed in zero
            graph.erase(a);
        }
    }

    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        if (SCIPisEQ(scip, SCIPvarGetUbLocal(x[nodeRef[v]]), 0))
            graph.erase(v);
    }

    //GraphViewer::ViewGLCIPSupportGraph(instance, graph, "Sub-graph", nodeRef);
}

/**
 * propagate in the topological ordering of condensed graphfor each condensed node, 
 * if the total of influence incident on it is less than the threshold, 
 * pay the difference (needed incentive) find the vertex of smaller threshold or
 * choose the vertices in the component who receives influence of the previous components
 */
double ExtendedDualBound::getCostInTopologicalOrdering(
    Digraph &condensed,
    int nComponents,
    DNodeValueMap &thr,
    ArcValueMap &arcWeight,
    vector<double> incentives)
{

    //linear time algorithm to solve the problem in DAGs
    double total = 0;
    //for (int i = 0; i < nComponents; i++)
    for (DNodeIt v(condensed); v != INVALID; ++v)
    {
        //cout << "visiting node " << i << endl;
        double sum = 0;
        for (InArcIt a(condensed, v); a != INVALID; ++a)
        {
            sum += arcWeight[a];
        }

        //cout << "condensed node id = " << condensed.id(condensed.nodeFromId(i)) << ", index = " << i << endl;

        //cout << "sum = " << sum << " thr = " << thr[i] << endl;
        if (sum < thr[v])
        {
            double p = 0;
            for (double j : incentives)
            {
                if (sum + j >= thr[v])
                {
                    p = j;
                    break;
                }
            }

            total += p;
            //cout << "paying incentive of: " << p << " to " << i << endl;
        }
    }

    return total;
}

//maps of gurobi variables
typedef Digraph::NodeMap<GRBVar> DNodeGRBVarMap;
typedef Digraph::NodeMap<vector<GRBVar>> DNodeGRBVarsMap;

//TODO implemente the same model exactWLCIPonDAG, but now considering only a set of possible incentives given as input
double ExtendedDualBound::exactWLCIPonDAG(
    SCIP *scip,
    Digraph &graph,
    ArcValueMap &arcWeight,
    DNodeValueMap &thr,
    DNodeIntMap &w,
    vector<double> incentives)
{
    double obj = 0; //stores the o objective function value

    //declare the variables
    DNodeGRBVarsMap xvp(graph);   // indicate what of the possible incentives is given to v
    DNodeGRBVarMap active(graph); // indicates if a vertex is active (equivalent to variable x of LCIP formulation)

    try
    {
        //model
        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "exact formulation of WLCIP on DAG");

        //create variables
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            string vname = "a_" + to_string(graph.id(v));
            active[v] = model.addVar(0, 1, 0, GRB_BINARY, vname);

            xvp[v].reserve(incentives.size());
            for (size_t i = 0; i < incentives.size(); i++)
            {
                string yname = "y_" + to_string(graph.id(v)) + to_string(incentives[i]);
                xvp[v][i] = model.addVar(0, 1, incentives[i], GRB_BINARY, yname);
            }
        }

        //the objective is to minimize the total incentive offered
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        //threshold constraints: the total influence comming to a vertex needs to be greater than
        //its threshold minus the incentive
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            GRBLinExpr lhs;

            for (size_t i = 0; i < incentives.size(); i++)
                lhs += incentives[i] * xvp[v][i];

            for (InArcIt a(graph, v); a != INVALID; ++a)
                lhs += arcWeight[a] * active[graph.source(a)];

            model.addConstr(lhs, GRB_GREATER_EQUAL, thr[v] * active[v], "thr-constraint");
        }

        //linking constraints
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            GRBLinExpr lhs;

            for (size_t i = 0; i < incentives.size(); i++)
                lhs += xvp[v][i];

            model.addConstr(lhs, GRB_EQUAL, active[v], "linking-constraint");
        }

        //cover constraints
        GRBLinExpr lhs;
        for (DNodeIt v(graph); v != INVALID; ++v)
            lhs += w[v] * active[v];
        model.addConstr(lhs, GRB_GREATER_EQUAL, ceil(instance.alpha * instance.n), "cover-constraint");

        model.set(GRB_IntParam_OutputFlag, 0);
        model.optimize();

        obj = model.get(GRB_DoubleAttr_ObjVal);

        //getting the active vertices of "graph"
        /* cout << "active vertices:";
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            if (SCIPisPositive(scip, active[v].get(GRB_DoubleAttr_X)))
            {
                cout << " " << to_string(graph.id(v));
            }
        }
        cout << "\n";
        cout << "incentives given:\n";
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            for (size_t i = 0; i < incentives.size(); i++)
            {
                if (SCIPisPositive(scip, xvp[v][i].get(GRB_DoubleAttr_X)))
                {
                    cout << "y_" << to_string(graph.id(v)) << " = " << incentives[i] << endl;
                }
            }
        }

        cout << "obj = " << obj << endl;
        cout << endl; */
    }
    catch (GRBException e)
    {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Error during optimization of GPC's separation" << endl;
    }

    return obj;
}

/**
 *  Integer linear program for the LCIP with weight on vertices (WLCIP) for alpha < 1
 */
double ExtendedDualBound::exactWLCIPonDAG(
    SCIP *scip,
    Digraph &graph,
    ArcValueMap &arcWeight,
    DNodeValueMap &thr,
    DNodeIntMap &w)
{
    double obj = 0; //stores the o objective function value

    //declare the variables
    DNodeGRBVarMap y(graph);      // incentive given to each vertex
    DNodeGRBVarMap active(graph); // indicates if a vertex is active (equivalent to variable x of LCIP formulation)

    try
    {
        //model
        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "exact formulation of WLCIP on DAG");

        //create variables
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            string vname = "a_" + to_string(graph.id(v));
            string yname = "y_" + to_string(graph.id(v));

            active[v] = model.addVar(0, 1, 0, GRB_BINARY, vname);
            y[v] = model.addVar(0, GRB_MAXINT, 1.0, GRB_INTEGER, yname);
        }

        //the objective is to minimize the total incentive offered
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        //threshold constraint: the total influence comming to a vertex needs to be greater than
        //its threshold minus the incentive
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            GRBLinExpr lhs;

            for (InArcIt a(graph, v); a != INVALID; ++a)
                lhs += arcWeight[a] * active[graph.source(a)];
            lhs += y[v];

            model.addConstr(lhs, GRB_GREATER_EQUAL, thr[v] * active[v], "thr-constraint");
        }

        //cover constraint
        GRBLinExpr lhs;
        for (DNodeIt v(graph); v != INVALID; ++v)
            lhs += w[v] * active[v];
        model.addConstr(lhs, GRB_GREATER_EQUAL, ceil(instance.alpha * instance.n), "cover-constraint");

        model.set(GRB_IntParam_OutputFlag, 0);
        model.optimize();

        obj = model.get(GRB_DoubleAttr_ObjVal);

        //getting the active vertices of "graph"
        /* cout << "active vertices:";
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            if (SCIPisPositive(scip, active[v].get(GRB_DoubleAttr_X)))
            {
                cout << " " << to_string(graph.id(v));
            }
        }
        cout << "\n";
        cout << "incentives given:\n";
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            if (SCIPisPositive(scip, y[v].get(GRB_DoubleAttr_X)))
            {
                cout << "y_" << to_string(graph.id(v)) << " = " << y[v].get(GRB_DoubleAttr_X) << endl;
            }
        }

        cout << "obj = " << obj << endl;
        cout << endl; */
    }
    catch (GRBException e)
    {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Error during optimization of GPC's separation" << endl;
    }

    return obj;
}

SCIP_DECL_RELAXEXEC(ExtendedDualBound::scip_exec)
{
    //cout << "RELAXEXEC()" << endl;

    SCIP_Real relaxval;

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    //get the support graph of the current feasible solution
    Digraph graph;

    DNodeDNodeMap nodeRef(graph); //save the reference to the original node
    ArcArcMap arcRef(graph);      //save the reference to the original arc
    getSubGraph(scip, graph, nodeRef, arcRef);

    vector<double> incentives;
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        incentives = instance.incentives[v];
        break; // only one arbitrary vertex is needed
    }

    if (stronglyConnected(graph))
    {
        cout << "assossiated subgraph is strongly connected\n";
        //find minimum threshold vertex
        DNode node = INVALID;
        double minThr = getMinimumThreshold(node);

        relaxval = minThr;
        for (double j : incentives)
        {
            if (j >= minThr)
            {
                relaxval = j;
                break;
            }
        }

        //printf("Heuristic lower bound = %g\n", relaxval);
        /* *lowerbound = relaxval;
        *result = SCIP_SUCCESS; */
    }
    else
    {
        int nComponents = countStronglyConnectedComponents(graph);
        /* cout << "assossiated subgraph isn't strongly connected: ";
        cout << nComponents << " components\n";

        cout << "current node of the three: " << SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) << endl; */

        Digraph condensed;
        vector<vector<DNode>> listOfComponents(nComponents);
        Digraph::NodeMap<int> components(graph);
        ArcValueMap arcWeight(condensed);
        DNodeValueMap thr(condensed);

        getCondensedGraph(condensed, graph, nodeRef, components, listOfComponents);
        getCondensedArcWeights(condensed, graph, arcRef, arcWeight, components);
        getCondensedThresholds(condensed, listOfComponents, nComponents, thr);
        //printCondensedArcs(condensed, arcWeight);

        if (instance.alpha < 1)
        {
            //cout << "alpha < 1 AND G'\n";

            //generate the weight of each condensed vertex
            DNodeIntMap w(condensed);

            for (DNodeIt v(condensed); v != INVALID; ++v)
            {
                w[v] = listOfComponents[condensed.id(v)].size();
                //cout << condensed.id(v) << "'s weight: " << w[v] << endl;
            }

            //implement the ILP model in gurobi to solve the subproblem
            //relaxval = exactWLCIPonDAG(scip, condensed, arcWeight, thr, w);
            relaxval = exactWLCIPonDAG(scip, condensed, arcWeight, thr, w, incentives);

            //TODO compare this relaxval with the obtained by the minimum threshold
            DNode node = INVALID;
            double minThr = getMinimumThreshold(node);

            double m = minThr;
            for (double j : incentives)
            {
                if (j >= minThr)
                {
                    m = j;
                    break;
                }
            }

            /* if (relaxval > m)
                cout << "the solution of the ILP model is greater than the incentive to the minimum threshol vertex\n"; */
            
        }
        else
        {
            //linear time algorithm to solve the problem in DAGs for alpha = 1
            relaxval = getCostInTopologicalOrdering(condensed, nComponents, thr, arcWeight, incentives);
        }
    }

    //printf("Heuristic lower bound = %g\n", relaxval);
    *lowerbound = relaxval;
    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}
