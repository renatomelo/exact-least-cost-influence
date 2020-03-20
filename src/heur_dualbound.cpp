#include "heur_dualbound.h"

HeurDualBound::HeurDualBound(
    SCIP *scip,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z,
    DNodeSCIPVarsMap &p_xip) : ObjRelax(scip,
                                        "heuristic-dual-bound",
                                        "Heuristic dual bound for GLCIP",
                                        -1.0,   //priority of the relaxator (negative: after LP, non-negative: before LP)
                                        10,     //frequency for calling relaxator
                                        FALSE), //Does the relaxator contain all cuts in the LP?
                               instance(p_instance),
                               x(p_x),
                               z(p_z),
                               xip(p_xip),
                               sol_(NULL)
{
}

SCIP_DECL_RELAXFREE(HeurDualBound::scip_free)
{
    //cout << "RELAXFREE()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINIT(HeurDualBound::scip_init)
{
    //cout << "RELAXEINIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXIT(HeurDualBound::scip_exit)
{
    //cout << "RELAXEXIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINITSOL(HeurDualBound::scip_initsol)
{
    //cout << "RELAXINITSOL()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXITSOL(HeurDualBound::scip_exitsol)
{
    //cout << "RELAXEXITSOL()" << endl;
    return SCIP_OKAY;
}

void getCondensedGraph(
    GLCIPInstance &instance,
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

vector<double> getCondensedThresholds(
    GLCIPInstance &instance,
    vector<vector<DNode>> &listOfComponents,
    int nComponents)
{
    vector<double> thr(nComponents);

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
        thr[i] = minThreshold;
        //printf("threshold(%d) = %.1f\n", i, minThreshold);
    }
    return thr;
}

// add the condensed arcs in the condensed graph
void getCondensedArcWeights(
    GLCIPInstance &instance,
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

void printCondensedArcs(Digraph &condensed, ArcValueMap &condensedInfluence)
{
    for (ArcIt a(condensed); a != INVALID; ++a)
    {
        int i = condensed.id(condensed.source(a));
        int j = condensed.id(condensed.target(a));

        printf("condensed arc: %d - > %d: %f", i, j, condensedInfluence[a]);
    }
}

//propagate in the graph using only the incentives selected in lowerBound()
SCIP_RETCODE setRelaxedSolutionSpreading(
    SCIP *scip,
    GLCIPInstance &instance,
    ArcSCIPVarMap &z,
    DNodeSCIPVarMap &x,
    DNode &node)
{
    // start wiht empty solution
    Digraph::NodeMap<set<DNode>> influencers(instance.g);
    set<DNode> actives;
    list<DNode> seeds;
    list<DNode> waiting;

    // initialize set with seed nodes
    seeds.push_back(node);
    actives.insert(node);

    // while the seed set is not empty try to activate non active vertices
    while (seeds.size() > 0)
    {
        DNode u = seeds.front();
        seeds.pop_front();
        //actives.insert(u);

        cout << "u: " << instance.nodeName[u] << endl;

        size_t prevSize = actives.size();

        for (OutArcIt a(instance.g, u); a != INVALID; ++a)
        {
            DNode v = instance.g.target(a);
            if (!actives.count(v) && SCIPisPositive(scip, SCIPgetVarSol(scip, z[a])))
            {
                influencers[v].insert(u);
                SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], 1.0));

                double exerterdInfluence = 0;
                for (DNode w : influencers[v])
                {
                    Arc e = findArc(instance.g, w, v);
                    assert(e != INVALID);

                    exerterdInfluence += instance.influence[e];
                }

                if (exerterdInfluence >= instance.threshold[v])
                {
                    cout << instance.nodeName[v] << " is inserted in seed set " << endl;
                    seeds.push_back(v);
                    actives.insert(v);

                    SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], 1.0));
                }
            }
        }

        if (seeds.empty() && actives.size() == prevSize && actives.size() < instance.n * instance.alpha)
        {
            cout << "no vertex activated in this iteration. Forcing some nodes to be active\n";
            for (OutArcIt a(instance.g, u); a != INVALID; ++a)
            {
                DNode v = instance.g.target(a);
                if (!actives.count(v) && SCIPisPositive(scip, SCIPgetVarSol(scip, z[a])))
                {
                    for (InArcIt b(instance.g, v); b != INVALID; ++b)
                    {
                        DNode y = instance.g.source(b);
                        if (y != u && !actives.count(v) && SCIPisPositive(scip, SCIPgetVarSol(scip, z[b])))
                        {
                            influencers[v].insert(y);
                            SCIP_CALL(SCIPsetRelaxSolVal(scip, z[b], 1.0));

                            actives.insert(y);
                            SCIP_CALL(SCIPsetRelaxSolVal(scip, x[y], 1.0));

                            waiting.push_back(y);

                            double exerterdInfluence = 0;
                            for (DNode w : influencers[v])
                            {
                                Arc e = findArc(instance.g, w, v);
                                assert(e != INVALID);

                                exerterdInfluence += instance.influence[e];
                            }

                            if (exerterdInfluence >= instance.threshold[v])
                            {
                                cout << instance.nodeName[v] << " is inserted in seed set *" << endl;
                                seeds.push_back(v);
                                actives.insert(v);
                                SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], 1.0));
                                break;
                            }
                        }
                    }
                }
            }

            cout << "size of waiting list = " << waiting.size() << endl;
        }
    }

    while (!waiting.empty())
    {
        DNode v = waiting.front();
        waiting.pop_front();

        double exerterdInfluence = 0;
        for (InArcIt a(instance.g, v); a != INVALID; ++a)
        {
            if (SCIPisPositive(scip, SCIPgetVarSol(scip, z[a])))
            {
                DNode u = instance.g.source(a);

                influencers[v].insert(u);
                SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], 1.0));
                cout << instance.nodeName[u] << " is influencer of waiting node " << instance.nodeName[v] << endl;

                exerterdInfluence += instance.influence[a];

                if (exerterdInfluence >= instance.threshold[v])
                {
                    cout << "threshold of waiting node achieved" << endl;
                    break;
                }
            }
        }
    }

    if (actives.size() < instance.alpha * instance.n)
        cout << "not all vertices were activated\n";
    cout << "active nodes: ";
    for (DNode v : actives)
    {
        cout << " " << instance.nodeName[v];
    }
    cout << endl;

    return SCIP_OKAY;
}

double getMinimumThreshold(GLCIPInstance &instance, DNode &node)
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

int getIndexOfChepeastIncentive(GLCIPInstance &instance, DNode &node)
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
}

void getSubGraph(
    SCIP *scip,
    GLCIPInstance &instance,
    Digraph &graph,
    DNodeDNodeMap &nodeRef,
    ArcArcMap &arcRef,
    ArcSCIPVarMap &z)
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

    //GraphViewer::ViewGLCIPSupportGraph(instance, graph, "Sub-graph", nodeRef);
}

/**
 * propagate in the topological ordering of condensed graphfor each condensed node, 
 * if the total of influence incident on it is less than the threshold, 
 * pay the difference (needed incentive) find the vertex of smaller threshold or
 * choose the vertices in the component who receives influence of the previous components
 */
double getCostInTopologicalOrdering(
    Digraph &condensed,
    int nComponents,
    vector<double> thr,
    ArcValueMap &arcWeight)
{
    // start wiht empty solution
    set<DNode> actives;
    /* Digraph::NodeMap<set<DNode>> influencers(condensed);
    list<DNode> seeds;

    //TODO add to seed set every node that has no input arc
    double cost = 0;
    InDegMap<Digraph> inDegree(condensed);
    for (DNodeIt v(condensed); v != INVALID; ++v)
    {
        if (inDegree[v] == 0)
        {
            //in the associated component find the smaller threshold vertex to be the seed
            DNode node = INVALID;
            // find the vertex in component 0 to pay the incentive and start spreading
            for (size_t i = 0; i < listOfComponents[condensed.id(v)].size(); i++)
            {
                DNode w = listOfComponents[condensed.id(v)][i];
                if (thr[condensed.id(v)] == instance.threshold[w])
                {
                    node = w;
                    break;
                }
            }

            assert(node != INVALID);

            //cout << "seed node: " << instance.nodeName[node] << endl;
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

            cost += instance.incentives[node][index];
            //cost += thr[condensed.id(v)];

            seeds.push_back(v);
            originalSeeds.push_back(node);
        }
    }

    //cout << "cost of seed nodes: " << cost << endl;

    // while the seed set is not empty try to activate non active vertices
    while (seeds.size() > 0)
    {
        DNode u = seeds.front();
        seeds.pop_front();
        actives.insert(u);

        for (OutArcIt a(condensed, u); a != INVALID; ++a)
        {
            DNode v = condensed.target(a);
            if (!actives.count(v))
            {
                influencers[v].insert(u);

                double exerterdInfluence = 0;
                for (DNode w : influencers[v])
                {
                    Arc e = findArc(condensed, w, v);
                    assert(e != INVALID);

                    exerterdInfluence += condensedInfluence[e];
                }

                cout << " exerterd influence: " << exerterdInfluence << endl;
                if (exerterdInfluence >= thr[condensed.id(v)])
                {
                    cout << condensed.id(v) << " is inserted in seed set" << endl;
                    seeds.push_back(v);
                }
            }
        }
    }

    //if not all vertices have been reached give more incentives
    if ((int)actives.size() < nComponents)
    {
        cout << "not all vertices have been reached\n";
        //visit all non-active vertices and pay the necessary incentive
        for (DNodeIt v(condensed); v != INVALID; ++v)
        {
            if (!actives.count(v))
            {
                double exerterdInfluence = 0;
                for (DNode w : influencers[v])
                {
                    Arc e = findArc(condensed, w, v);
                    assert(e != INVALID);

                    exerterdInfluence += condensedInfluence[e];
                }
                cout << condensed.id(v) << " not active: ";
                cout << "exerterdInfluence = " << exerterdInfluence << " thr = " << thr[condensed.id(v)] << endl;
                double dif = thr[condensed.id(v)] - exerterdInfluence;
                cout << "difference = " << dif << endl;
                //TODO finish later if needed
            }
        }
        exit(0);
    } */

    //linear time algorithm to solve the problem in DAGs
    vector<double> incentives(nComponents);
    double total = 0;
    for (int i = 0; i < nComponents; i++)
    {
        //cout << "visiting node " << i << endl;
        double sum = 0;
        for (InArcIt a(condensed, condensed.nodeFromId(i)); a != INVALID; ++a)
        {
            sum += arcWeight[a];
        }

        //cout << "condensed node id = " << condensed.id(condensed.nodeFromId(i)) << ", index = " << i << endl;

        //cout << "sum = " << sum << " thr = " << thr[i] << endl;
        if (sum >= thr[i])
        {
            actives.insert(condensed.nodeFromId(i));
            incentives[i] = 0;
        }
        else
        {
            incentives[i] = thr[i] - sum;
            actives.insert(condensed.nodeFromId(i));
            //cout << "paying incentive of: " << incentives[i] << " to " << i << endl;
        }
        total += incentives[i];
    }

    /* if (total > incentives[0])
    {
        //show the diference and study what happens
        for (int i = 0; i < nComponents; i++)
        {
            if (incentives[i] > 0)
            {
                cout << "paying incentive of: " << incentives[i] << " to " << i << endl;
            }
        }
        cout << "total incentives = " << total << endl;
    } */

    //cout << "total incentives = " << total << endl;
    //cout << "size of actives = " << actives.size() << endl;
    //return incentives[0];
    return total;
}

bool isIntegral(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z)
{
    bool integral = TRUE;
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        if (!SCIPisIntegral(scip, SCIPgetVarSol(scip, z[a])))
        {
            integral = FALSE;
            //cout << SCIPvarGetName(z[a]) << " = " << SCIPgetVarSol(scip, z[a]) << endl;
        }
    }

    if (!integral)
    {
        cout << "fractional solution\n";
    }
    return integral;
}

SCIP_RETCODE setRelaxedSol(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    DNodeSCIPVarsMap &xip,
    set<DNode> seeds)
{
    //cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
    SCIP_CALL(SCIPclearRelaxSolVals(scip));

    //it is necessary to set the new solution? I think so
    for (DNode u : seeds)
    {
        SCIP_CALL(SCIPsetRelaxSolVal(scip, x[u], 1.0));

        //int index = getIndexOfChepeastIncentive(instance, u);
        //SCIP_CALL(SCIPsetRelaxSolVal(scip, xip[u][index], 1.0));

        for (DNodeIt v(instance.g); v != INVALID; ++v)
            if (v != u)
                SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], SCIPgetVarSol(scip, x[v])));

        for (ArcIt a(instance.g); a != INVALID; ++a)
        {
            //arcs pointing to the seed node are not selected
            if (instance.g.target(a) != u)
                SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], SCIPgetVarSol(scip, z[a])));
        }
    }
    // propagate in the graph using only the incentives selected in lowerBound()

    //TODO: if we found a strongly connected component add a cuting plane

    //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
    SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
    //cout << "SCIPgetRelaxSolObj(scip) = " << SCIPgetRelaxSolObj(scip) << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    //cout << "RELAXEXEC()" << endl;

    SCIP_Real relaxval;

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    //get the support graph of the current feasible solution
    Digraph graph;

    DNodeDNodeMap nodeRef(graph); //save the reference to the original node
    ArcArcMap arcRef(graph);      //save the reference to the original arc
    getSubGraph(scip, instance, graph, nodeRef, arcRef, z);

    if (stronglyConnected(graph))
    {
        //cout << "support graph is strongly connected\n";
        //find minimum threshold vertex
        DNode node = INVALID;
        relaxval = getMinimumThreshold(instance, node);

        //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
        /* if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            set<DNode> seed;
            seed.insert(node);
            SCIP_CALL(setRelaxedSol(scip, instance, x, z, xip, seed));
        } */
        //printf("Heuristic lower bound = %g\n", relaxval);
        *lowerbound = relaxval;
        *result = SCIP_SUCCESS;
    }
    else
    {
        int nComponents = countStronglyConnectedComponents(graph);
        //cout << "support graph isn't strongly connected: ";
        //cout << nComponents << " components\n";

        Digraph condensed;
        vector<vector<DNode>> listOfComponents(nComponents);
        vector<double> thr(nComponents);
        Digraph::NodeMap<int> components(graph);
        ArcValueMap arcWeight(condensed);

        getCondensedGraph(instance, condensed, graph, nodeRef, components, listOfComponents);
        getCondensedArcWeights(instance, condensed, graph, arcRef, arcWeight, components);
        thr = getCondensedThresholds(instance, listOfComponents, nComponents);
        //printCondensedArcs(condensed, arcWeight);

        //linear time algorithm to solve the problem in DAGs
        relaxval = getCostInTopologicalOrdering(condensed, nComponents, thr, arcWeight);

        //printf("Heuristic lower bound = %g\n", relaxval);
        *lowerbound = relaxval;
        *result = SCIP_SUCCESS;
    }

    return SCIP_OKAY;
}