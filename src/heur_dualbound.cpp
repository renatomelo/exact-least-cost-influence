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

double getLowerBound(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    DNodeSCIPVarsMap &xip,
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
                    //cout << "adding condensed arc: " << i << " -> " << j << ": " << instance.influence[arcRef[a]] << endl;
                }
                else
                    condensedInfluence[b] += instance.influence[arcRef[a]];
            }
        }

        /* for (ArcIt a(condensed); a != INVALID; ++a)
        {
            cout << "condensed arc: " << condensed.id(condensed.source(a)) << " -> " 
                 << condensed.id(condensed.target(a)) << ": " << condensedInfluence[a] << endl;
        } */
    }

    exit(0);
    return minCost;
}

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    cout << "RELAXEXEC()" << endl;

    SCIP_Real relaxval;
    //SCIP_Bool valid;

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    //call the heuristic here
    DNode node = INVALID;
    relaxval = getLowerBound(scip, instance, x, z, xip, node);
    //cout << "selected node: " << instance.nodeName[node] << endl;
    //cout << "relaxval = " << relaxval << endl;
    int idx = 0;
    for (size_t i = 0; i < instance.incentives[node].size(); i++)
    {
        if (instance.incentives[node][i] >= instance.threshold[node])
        {
            cout << "incentive paid: " << instance.incentives[node][i] << endl;
            idx = i;
            break;
        }
    }

    //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
    /* if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
    {
        cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
        SCIP_CALL(SCIPclearRelaxSolVals(scip));

        //it is necessary to set the new solution? I think so
        SCIP_CALL(SCIPsetRelaxSolVal(scip, x[node], 1.0));
        SCIP_CALL(SCIPsetRelaxSolVal(scip, xip[node][idx], 1.0));

        for (DNodeIt v(instance.g); v != INVALID; ++v)
            if (v != node)
                SCIP_CALL(SCIPsetRelaxSolVal(scip, x[v], SCIPgetVarSol(scip, x[v])));

        for (ArcIt a(instance.g); a != INVALID; ++a)
        {
            //not chosing the arcs that pointing to the seed node
            if (instance.g.target(a) != node)
                SCIP_CALL(SCIPsetRelaxSolVal(scip, z[a], SCIPgetVarSol(scip, z[a])));
        }
        // propagate in the graph using only the incentives selected in lowerBound()

        //TODO: if we found a strongly connected component add a cuting plane

        //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
        SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
    }

    printf("Heuristic lower bound = %g\n", relaxval);
    *lowerbound = relaxval;
    *result = SCIP_SUCCESS; */

    return SCIP_OKAY;
}