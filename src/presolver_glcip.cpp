#include "presolver_glcip.h"

PresolverGLCIP::PresolverGLCIP(
    SCIP *scip,
    GLCIPInstance &p_instance, // problem data
    DNodeSCIPVarMap &p_x,
    ArcSCIPVarMap &p_z
    //DNodeInfSetsMap &p_inf_set // influencing set data structure
    ) : ObjPresol(scip, 
                  "preprocessing-rule",
                  "Preprocessing rules",
                  10,                    //priority of the presolver
                  1,                     //maximal number of presolving rounds the presolver participates in (-1: no limit)
                  SCIP_PRESOLTIMING_FAST //timing mask of the presolver
                  ),
        instance(p_instance),
        x(p_x),
        z(p_z)
        //infSet(p_inf_set)
{
}

PresolverGLCIP::~PresolverGLCIP() {}

void printArc(GLCIPInstance &instance, Arc a)
{
    cout << "(" << instance.nodeName[instance.g.source(a)] << ", " 
         << instance.nodeName[instance.g.target(a)] << ")" << endl;
}

// see later: https://www.scipopt.org/doc/html/presol__trivial_8c_source.php#l00067
SCIP_DECL_PRESOLEXEC(PresolverGLCIP::scip_exec)
{
    //cout << "--------- PRESOLVER -----------------\n";
    int count = 0;
    SCIP_Bool infeasible;
    SCIP_Bool fixed;

    vector<Arc> fixedInOne;
    vector<Arc> fixedInZero;

    assert(result != NULL);

    *result = SCIP_DIDNOTRUN;

    if (instance.alpha < 1)
        return SCIP_OKAY;

    Digraph graph;
    DNodeDNodeMap nodeRef(graph); //save the reference to the original node
    ArcArcMap arcRef(graph);      //save the reference to the original arc

    digraphCopy(instance.g, graph).nodeCrossRef(nodeRef).arcCrossRef(arcRef).run();

    //create a copy of the thresholds
    DNodeValueMap thr(graph);
    for (DNodeIt v(graph); v != INVALID; ++v)
        thr[v] = instance.threshold[nodeRef[v]];

    InDegMap<Digraph> inDegree(graph);
    OutDegMap<Digraph> outDegree(graph);

    vector<DNode> toRemove; //save the vertices to be removed in the preprocessing

    for (DNodeIt v(graph); v != INVALID; ++v)
    {
        if (inDegree[v] == 0 && outDegree[v] == 0)
            toRemove.push_back(v);
        else if (inDegree[v] == 0)
        {
            for (OutArcIt a(graph, v); a != INVALID; ++a)
            {
                fixedInOne.push_back(arcRef[a]);

                DNode w = graph.target(a);
                thr[w] -= instance.influence[arcRef[a]];
                /* cout << "to fix arc: ";
                printArc(instance, arcRef[a]); */
            }
            count++;
            //TODO fix the inf_set var too

            toRemove.push_back(v);
            //cout << "to remove vertex " << instance.nodeName[nodeRef[v]] << endl;
        }
        else if (outDegree[v] == 0)
        {
            //cout << instance.nodeName[v] << "'s (out, in)-degree = "  << outDegree[v] << ", " << inDegree[v] << endl;
            for (InArcIt a(graph, v); a != INVALID; ++a)
            {
                fixedInOne.push_back(arcRef[a]);
                /* cout << "to fix arc: ";
                printArc(instance, arcRef[a]); */
            }
            count++;
            //TODO fix the inf_set var too

            toRemove.push_back(v);
            //cout << "to remove vertex " << instance.nodeName[nodeRef[v]] << endl;
        }
    }

    for (DNode v : toRemove)
    {
        //cout << "removing vertex " << instance.nodeName[nodeRef[v]] << endl;
        graph.erase(v);
    }

    do
    {
        toRemove.clear();

        for (DNodeIt v(graph); v != INVALID; ++v)
        {
            if (thr[v] <= 0)
            {
                for (InArcIt a(graph, v); a != INVALID; ++a)
                {
                    fixedInZero.push_back(arcRef[a]);
                    //cout << "to fix arc in zero: ";
                    //printArc(instance, arcRef[a]);
                }
                for (OutArcIt a(graph, v); a != INVALID; ++a)
                {
                    fixedInOne.push_back(arcRef[a]);

                    DNode w = graph.target(a);
                    thr[w] -= instance.influence[arcRef[a]];

                    /* cout << "to fix arc: ";
                    printArc(instance, arcRef[a]); */
                }
                count++;
                //TODO fix the inf_set var too
                //cout << "to remove vertex " << instance.nodeName[nodeRef[v]] << endl;
                toRemove.push_back(v);
                break;
            }
            else if (outDegree[v] == 0)
            {
                //cout << instance.nodeName[v] << "'s (out, in)-degree = "  << outDegree[v] << ", " << inDegree[v] << endl;
                for (InArcIt a(graph, v); a != INVALID; ++a)
                {
                    fixedInOne.push_back(arcRef[a]);
                    /* cout << "to fix arc: ";
                    printArc(instance, arcRef[a]); */
                }
                count++;
                //TODO fix the inf_set var too

                //cout << "to remove vertex " << instance.nodeName[nodeRef[v]] << endl;
                toRemove.push_back(v);
            }
        }

        for (DNode v : toRemove)
        {
            //cout << "removing vertex " << instance.nodeName[nodeRef[v]] << endl;
            graph.erase(v);
        }

    } while (!toRemove.empty());

    //fix variables
    for (Arc a : fixedInOne)
    {
        //cout << "fixing variable in one " << SCIPvarGetName(z[a]) << endl;
        SCIP_CALL(SCIPfixVar(scip, z[a], 1, &infeasible, &fixed));
        if (infeasible)
        {
            cout << " -> infeasible fixing\n";
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
        }
        assert(fixed);
        (*nfixedvars)++;
        *result = SCIP_SUCCESS;
    }
    for (Arc a : fixedInZero)
    {
        //cout << "fixing variable in zero " << SCIPvarGetName(z[a]) << endl;
        SCIP_CALL(SCIPfixVar(scip, z[a], 0, &infeasible, &fixed));
        if (infeasible)
        {
            cout << " -> infeasible fixing\n";
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
        }
        assert(fixed);
        (*nfixedvars)++;
        *result = SCIP_SUCCESS;
    }

    //cout << "count = " << count << endl;
    //exit(0);

    return SCIP_OKAY;
}