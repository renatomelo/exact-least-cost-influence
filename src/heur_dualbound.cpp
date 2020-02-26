#include "heur_dualbound.h"

HeurDualBound::HeurDualBound(
    SCIP *scip,
    GLCIPInstance &p_instance,
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z,
    DNodeSCIPVarsMap &p_xip) : ObjRelax(scip,
                                        "heuristic-dual-bound",
                                        "Heuristic dual bound for GLCIP",
                                        -1.0,  //priority of the relaxator (negative: after LP, non-negative: before LP)
                                        10,    //frequency for calling relaxator
                                        TRUE), //Does the relaxator contain all cuts in the LP?
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

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    //cout << "RELAXEXEC()" << endl;

    //SCIP_Bool *cutoff = 0;
    //cout << "SCIPconstructLP(scip, cutoff) = " << SCIPconstructLP(scip, cutoff) << endl;
    /* cout << "SCIPgetDualbound(scip) = " << SCIPgetDualbound(scip) << endl;
    cout << "SCIPgetRelaxSolObj(scip) = " << SCIPgetRelaxSolObj(scip) << endl;
    cout << "SCIPgetFocusNode(scip) = " << SCIPnodeGetNumber(SCIPgetFocusNode(scip)) << endl;
    cout << "SCIPgetNodeDualbound() = " << SCIPgetNodeDualbound(scip, SCIPgetFocusNode(scip)) << endl;
    cout << "SCIPgetDualboundRoot(scip) = " << SCIPgetDualboundRoot(scip) << endl;
    cout << "*lowerbound = " << *lowerbound << endl; */

    SCIP_Real relaxval;

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    double minCost = SCIPinfinity(scip);

    //get the support graph of the current feasible solution
    Digraph graph;
    //GLCIPBase::getSuportGraph(scip, instance, NULL, z, graph);

    DNodeDNodeMap nodeRef(graph); //save the reference to the original node
    ArcArcMap arcRef(graph);      //save the reference to the original arc
    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    for (ArcIt a(graph); a != INVALID; ++a)
    {
        if (SCIPisEQ(scip, SCIPgetVarSol(scip, z[arcRef[a]]), 0))
        {
            graph.erase(a);
        }
    }

    //GraphViewer::ViewGLCIPSupportGraph(instance, graph, "Support Graph", nodeRef);

    if (stronglyConnected(graph))
    {
        //cout << "support graph is strongly connected\n";
        //find minimum threshold vertex
        DNode node = INVALID;
        for (DNodeIt v(instance.g); v != INVALID; ++v)
        {
            //double cost = GLCIPBase::cheapestIncentive(instance, v, 0);
            //cout << "threshold of " << instance.nodeName[v] << ": " << instance.threshold[v] << endl;
            //cout << "cost of " << instance.nodeName[v] << ": " << cost << endl;
            double thr = instance.threshold[v];
            if (thr < minCost)
            {
                minCost = thr;
                node = v;

                // stop if thr(v) = 1 because there is no smaller threshold
                if (instance.threshold[v] == 1)
                    break;
            }
        }

        relaxval = minCost;
        //relaxval = instance.threshold[node];

        //cout << "selected node: " << instance.nodeName[node] << endl;
        /* int index = 0;
        for (size_t i = 0; i < instance.incentives[node].size(); i++)
        {
            if (instance.incentives[node][i] >= instance.threshold[node])
            {
                //cout << "incentive paid: " << instance.incentives[node][i] << endl;
                index = i;
                break;
            }
        } */

        //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
        /* if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            //cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
            SCIP_CALL(SCIPclearRelaxSolVals(scip));

            //it is necessary to set the new solution? I think so
            //SCIP_CALL(SCIPsetRelaxSolVal(scip, x[node], 1.0));
            //SCIP_CALL(SCIPsetRelaxSolVal(scip, xip[node][index], 1.0));

            for (DNodeIt v(instance.g); v != INVALID; ++v)
                //if (v != node)
                SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], SCIPgetVarSol(scip, x[v])));

            for (ArcIt a(instance.g); a != INVALID; ++a)
            {
                //arcs pointing to the seed node are not selected
                //if (instance.g.target(a) != node)
                    SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], SCIPgetVarSol(scip, z[a])));
            } */

            // propagate in the graph using only the incentives selected in lowerBound()

            // start wiht empty solution
            /* Digraph::NodeMap<set<DNode>> influencers(instance.g);
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
            cout << endl; */

            //TODO: if we found a strongly connected component add a cuting plane

            //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
            /* SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
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
        for (int i = 0; i < nComponents; i++)
        {
            condensed.addNode();
            //DNode v = condensed.addNode();
            //cout << "adding condensed node: " << condensed.id(v) << endl;
        }

        //compute the components of 'graph'
        Digraph::NodeMap<int> components(graph);
        stronglyConnectedComponents(graph, components);

        //save each component in a vector of vertices
        vector<vector<DNode>> listOfComponents(nComponents);
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            listOfComponents[components[v]].push_back(nodeRef[v]);
        }

        vector<double> thr(nComponents);
        //find the minimum threshold on each component
        for (int i = 0; i < nComponents; i++)
        {
            double minThreshold = SCIPinfinity(scip);
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

        //get the cut arcs of the strongly connected components
        Digraph::ArcMap<bool> cutArcs(graph, FALSE);
        stronglyConnectedCutArcs(graph, cutArcs);

        // add the condensed arcs in the condensed graph
        ArcValueMap condensedInfluence(condensed);
        for (ArcIt a(graph); a != INVALID; ++a)
        {
            //compute the weigh of influence of the arcs
            //each arc receives the total of weights of all arcs from a component to another
            if (cutArcs[a])
            {
                //reference to what component each vertex belong
                int i = components[graph.source(a)];
                int j = components[graph.target(a)];

                Arc b = findArc(condensed, condensed.nodeFromId(i), condensed.nodeFromId(j));
                if (b == INVALID)
                {
                    Arc c = condensed.addArc(condensed.nodeFromId(i), condensed.nodeFromId(j));
                    condensedInfluence[c] = instance.influence[arcRef[a]];
                }
                else
                    condensedInfluence[b] += instance.influence[arcRef[a]];
            }
        }

        /* for (ArcIt a(condensed); a != INVALID; ++a)
            cout << "condensed arc: " << condensed.id(condensed.source(a)) << " -> "
                 << condensed.id(condensed.target(a)) << ": " << condensedInfluence[a] << endl; */

        //propagate in the topological ordering of condensed graph
        //for each condensed node, if the total of influence incident on
        //it is less than the threshold, pay the difference (needed incentive)
        //find the vertex of smaller threshold or choose the vertices in the
        //component who receives influence of the previous components

        // start wiht empty solution
        Digraph::NodeMap<set<DNode>> influencers(condensed);
        set<DNode> actives;
        list<DNode> seeds;
        list<DNode> originalSeeds;

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
                /* int index = 0;
                for (size_t i = 0; i < instance.incentives[node].size(); i++)
                {
                    if (instance.incentives[node][i] >= instance.threshold[node])
                    {
                        //cout << "incentive paid: " << instance.incentives[node][i] << endl;
                        index = i;
                        break;
                    }
                }

                cost += instance.incentives[node][index]; */
                cost += thr[condensed.id(v)];

                seeds.push_back(v);
                originalSeeds.push_back(node);
            }
        }

        //cout << "cost of seed nodes: " << cost << endl;
        relaxval = cost;

        //linear time algorithm to solve the problem in DAGs
        vector<double> incentives(nComponents);
        double total = 0;
        for (int i = 0; i < nComponents; i++)
        {
            //cout << "visiting node " << i << endl;
            double sum = 0;
            for (InArcIt a(condensed, condensed.nodeFromId(i)); a != INVALID; ++a)
            {
                sum += condensedInfluence[a];
            }

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

        //cout << "total incentives = " << total << endl;
        //cout << "size of actives = " << actives.size() << endl;

        // while the seed set is not empty try to activate non active vertices
        /* while (seeds.size() > 0)
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
                    cout << "exerterdInfluence = " << exerterdInfluence << " thr = "<< thr[condensed.id(v)] << endl;
                    double dif = thr[condensed.id(v)] - exerterdInfluence;
                    cout << "difference = " << dif << endl;
                    //TODO finish later if needed
                }
            }
            exit(0);
        } */

        relaxval = total;

        //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
        /* if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            //cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
            SCIP_CALL(SCIPclearRelaxSolVals(scip));

            //it is necessary to set the new solution? I think so
            for (DNode u : originalSeeds)
            {
                //SCIP_CALL(SCIPsetRelaxSolVal(scip, x[u], 1.0));

                int index = 0;
                for (size_t i = 0; i < instance.incentives[u].size(); i++)
                {
                    if (instance.incentives[u][i] >= instance.threshold[u])
                    {
                        //cout << "incentive paid: " << instance.incentives[u][i] << endl;
                        index = i;
                        break;
                    }
                }

                SCIP_CALL(SCIPsetRelaxSolVal(scip, xip[u][index], 1.0));

                for (DNodeIt v(instance.g); v != INVALID; ++v)
                    //if (v != u)
                        SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], SCIPgetVarSol(scip, x[v])));

                for (ArcIt a(instance.g); a != INVALID; ++a)
                {
                    //arcs pointing to the seed node are not selected
                    //if (instance.g.target(a) != u)
                        SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], SCIPgetVarSol(scip, z[a])));
                }
            }
            // propagate in the graph using only the incentives selected in lowerBound()

            //TODO: if we found a strongly connected component add a cuting plane

            //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
            SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
            //cout << "SCIPgetRelaxSolObj(scip) = " << SCIPgetRelaxSolObj(scip) << endl;
        }
        else
            cout << "relaxed solution didn't improve the corrent relaxation\n"; */

        //printf("Heuristic lower bound = %g\n", relaxval);
        *lowerbound = relaxval;
        *result = SCIP_SUCCESS;
    }

    return SCIP_OKAY;
}