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

/* void getSuportGraph(
    SCIP *scip,
    GLCIPInstance &instance,
    SCIP_SOL *sol,
    ArcSCIPVarMap &z,
    Digraph &new_graph)
{
    DNodeDNodeMap nodeRef(new_graph);
    ArcArcMap arcRef(new_graph);
    digraphCopy(instance.g, new_graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    for (ArcIt a(new_graph); a != INVALID; ++a)
    {
        if (SCIPisEQ(scip, SCIPgetSolVal(scip, sol, z[arcRef[a]]), 0))
        {
            new_graph.erase(a);
        }
    }

    GraphViewer::ViewGLCIPSupportGraph(instance, new_graph, "Support Graph", nodeRef);
} */

double getLowerBound(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z)
{
    double minCost = -SCIPinfinity(scip);
    //get the support graph of the current feasible solution
    Digraph graph;
    GLCIPBase::getSuportGraph(scip, instance, NULL, z, graph);

    if(stronglyConnected(graph))
    {
        cout << "support graph is strongly connected\n";
        return 2.0;
    }
    else
    {
        cout << "support graph isn't strongly connected: ";
        cout << countStronglyConnectedComponents(graph) << " components\n";
        Digraph::NodeMap<int> components(graph);
        stronglyConnectedComponents(graph, components);
/*         for (DNodeIt v(graph); v != INVALID; ++v)
        {
            cout << instance.nodeName[v] << " is in component " << components[v] << endl;
        } */

        Digraph condensed;
        for (int i = 0; i < countStronglyConnectedComponents(graph); i++)
        {
            DNode v = condensed.addNode();
            cout << "creating condensed node: " << condensed.id(v) << endl;
        }
        
    }
    
    //exit(0);
    return minCost;
}

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    cout << "RELAXEXEC()" << endl;

    SCIP_Real relaxval;
    SCIP_Bool valid;

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    //call the heuristic here
    relaxval = getLowerBound(scip, instance, x, z);

    //store relaxation solution in original SCIP if it improves the best relaxation solution thus far
    if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
    {
        cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
        SCIP_CALL(SCIPclearRelaxSolVals(scip));

        //it is necessary to set the new solution? I think so

        //TODO: if we found a strongly connected component add a cuting plane

        //mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows
        SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
    }

    printf("Heuristic lower bound = %g\n", relaxval);
    *lowerbound = relaxval;
    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}