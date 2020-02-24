#include "heur_dualbound.h"

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

/* double getLowerBound(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    DNode &node)
{
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

    //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph", nodeRef);

    if (stronglyConnected(graph))
    {
        cout << "support graph is strongly connected\n";
        //find minimum threshold vertex
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            double cost = GLCIPBase::cheapestIncentive(instance, nodeRef[v], 0);
            if (cost < minCost)
            {
                minCost = cost;
                node = nodeRef[v];

                // stop if thr(v) = 1 because there is no smaller threshold
                if (instance.threshold[nodeRef[v]] == 1)
                    break;
            }
        }
    }
    else
    {
        int nComponents = countStronglyConnectedComponents(graph);
        cout << "support graph isn't strongly connected: ";
        cout << nComponents << " components\n";

        Digraph condensed;
        for (int i = 0; i < nComponents; i++)
        {
            DNode v = condensed.addNode();
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
            printf("size of component %d: %ld\n", i, listOfComponents[i].size());
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
                    //cout << "adding condensed arc: " << i << " -> " << j << ": " << instance.influence[arcRef[a]] << endl;
                }
                else
                    condensedInfluence[b] += instance.influence[arcRef[a]];
            }
        }

        for (ArcIt a(condensed); a != INVALID; ++a)
        {
            cout << "condensed arc: " << condensed.id(condensed.source(a)) << " -> "
                 << condensed.id(condensed.target(a)) << ": " << condensedInfluence[a] << endl;
        }
    }

    //exit(0);
    return minCost;
} */

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    cout << "RELAXEXEC()" << endl;

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

    //GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph", nodeRef);

    if (stronglyConnected(graph))
    {
        cout << "support graph is strongly connected\n";
        //find minimum threshold vertex
        DNode node = INVALID;
        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            double cost = GLCIPBase::cheapestIncentive(instance, nodeRef[v], 0);
            if (cost < minCost)
            {
                minCost = cost;
                node = nodeRef[v];

                // stop if thr(v) = 1 because there is no smaller threshold
                if (instance.threshold[nodeRef[v]] == 1)
                    break;
            }
        }

        relaxval = minCost;

        //cout << "selected node: " << instance.nodeName[node] << endl;
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

        //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
        if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
            SCIP_CALL(SCIPclearRelaxSolVals(scip));

            //it is necessary to set the new solution? I think so
            SCIP_CALL(SCIPsetRelaxSolVal(scip, x[node], 1.0));
            SCIP_CALL(SCIPsetRelaxSolVal(scip, xip[node][index], 1.0));

            for (DNodeIt v(instance.g); v != INVALID; ++v)
                if (v != node)
                    SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], SCIPgetVarSol(scip, x[v])));

            for (ArcIt a(instance.g); a != INVALID; ++a)
            {
                //arcs pointing to the seed node are not selected
                if (instance.g.target(a) != node)
                    SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], SCIPgetVarSol(scip, z[a])));
            }
            // propagate in the graph using only the incentives selected in lowerBound()

            //TODO: if we found a strongly connected component add a cuting plane

            //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
            SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));

            printf("Heuristic lower bound = %g\n", relaxval);
            *lowerbound = relaxval;
            *result = SCIP_SUCCESS;
        }
    }
    else
    {
        int nComponents = countStronglyConnectedComponents(graph);
        cout << "support graph isn't strongly connected: ";
        cout << nComponents << " components\n";

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

        //for (ArcIt a(condensed); a != INVALID; ++a)
        //    cout << "condensed arc: " << condensed.id(condensed.source(a)) << " -> "
        //         << condensed.id(condensed.target(a)) << ": " << condensedInfluence[a] << endl;

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

                cout << "seed node: " << instance.nodeName[node] << endl;
                int index = 0;
                for (size_t i = 0; i < instance.incentives[node].size(); i++)
                {
                    if (instance.incentives[node][i] >= instance.threshold[node])
                    {
                        cout << "incentive paid: " << instance.incentives[node][i] << endl;
                        index = i;
                        break;
                    }
                }

                cost += instance.incentives[node][index];

                seeds.push_back(v);
                originalSeeds.push_back(node);
            }
        }

        //cout << "cost of seed nodes: " << cost << endl;
        relaxval = cost;

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
                    //std::cout << " exerterd influence + incentive: "
                    //<< (exerterdInfluence + solution.incentives[v]) << std::endl;
                    if (exerterdInfluence >= thr[condensed.id(v)])
                    {
                        //cout << condensed.id(v) << " is inserted in seed set" << endl;
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
                    cout << condensed.id(v) << " not active\n";
                    //double dif = thr[condensed.id(v)] - exerterdInfluence;
                    //TODO finish later if needed
                }
            }
            exit(0);
        }

        //relaxval = minCost;

        //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
        if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
            SCIP_CALL(SCIPclearRelaxSolVals(scip));

            //it is necessary to set the new solution? I think so
            for (DNode u : originalSeeds)
            {
                SCIP_CALL(SCIPsetRelaxSolVal(scip, x[u], 1.0));

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
            cout << "SCIPgetRelaxSolObj(scip) = " << SCIPgetRelaxSolObj(scip) << endl;
            printf("Heuristic lower bound = %g\n", relaxval);
            *lowerbound = relaxval;
            *result = SCIP_SUCCESS;
        }
        else
            cout << "relaxed solution didn't improve the corrent relaxation\n";

        // test fucntions here
    }

    return SCIP_OKAY;
}