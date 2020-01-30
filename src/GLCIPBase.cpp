#include "GLCIPBase.h"

// add a cycle founded by the bfs algorithm
void GLCIPBase::addCycleConstraints(SCIP *scip,
                                    GLCIPInstance &instance,
                                    DNodeSCIPVarMap &x,
                                    ArcSCIPVarMap &z,
                                    DNodeIntMap &predMap,
                                    Arc &backArc)
{
    DNode s = instance.g.source(backArc);
    DNode t = instance.g.target(backArc);

    // walk from s to t adding corresponding nodes and arcs
    std::vector<Arc> arcs;
    std::vector<DNode> nodes;

    DNode v = s;
    nodes.push_back(v);
    arcs.push_back(backArc);

    do
    {
        DNode u = instance.g.nodeFromId(predMap[v]);
        Arc a = findArc(instance.g, u, v);
        arcs.push_back(a);
        nodes.push_back(u);
        v = u;
    } while (v != t);

    //cout << "small cycle found: ";
    // now add the corresponding constraints
    for (unsigned int i = 0; i < nodes.size(); i++)
    {
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0, "small-cycle cons");

        // add arcs in the cycle
        for (auto a : arcs)
        {
            cons->addVar(z[a], 1.0);
        }

        // add all vertices except one
        for (unsigned int j = 0; j < nodes.size(); j++)
        {
            if (j != i)
            {
                DNode w = nodes[j];
                cons->addVar(x[w], -1.0);
            }
        }
        //cout << instance.nodeName[nodes[i]] << " ";

        cons->commit();
    }
    //cout << endl;
}

// run a depth first search with maximum height 3
void GLCIPBase::dfsSmallCycles(SCIP *scip,
                               GLCIPInstance &instance,
                               DNodeSCIPVarMap &x,
                               ArcSCIPVarMap &z,
                               DNodeIntMap &colors,
                               DNodeIntMap &predMap,
                               DNode curr,
                               int level,
                               int rootId)
{
    // visit current node (grey color)
    colors[curr] = 1;

    // pass through current node neighbors
    for (OutArcIt a(instance.g, curr); a != INVALID; ++a)
    {
        DNode next = instance.g.target(a);

        // if next is grey, then we found a cycle
        if (colors[next] == 1)
        {
            addCycleConstraints(scip, instance, x, z, predMap, a);
        }

        // if next is white, then we check the level
        if (colors[next] == 0 && level < 3 && instance.g.id(next) > rootId)
        {
            predMap[next] = instance.g.id(curr);
            dfsSmallCycles(scip, instance, x, z, colors, predMap, next, level + 1, rootId);
        }
    }

    // close this vertex (black color)
    colors[curr] = 2;
}

// add cycle removal constraints for cycles of size up to 4
void GLCIPBase::addSmallCycleConstraints(SCIP *scip,
                                         GLCIPInstance &instance,
                                         DNodeSCIPVarMap &x,
                                         ArcSCIPVarMap &z)
{
    for (DNodeIt r(instance.g); r != INVALID; ++r)
    {
        // run dfs algorithm for 3 levels
        DNodeIntMap predMap(instance.g);
        DNodeIntMap colors(instance.g);
        mapFill(instance.g, predMap, -1);
        mapFill(instance.g, colors, 0);

        predMap[r] = instance.g.id(r);
        dfsSmallCycles(scip, instance, x, z, colors, predMap, r, 0, instance.g.id(r));
    }
}

void addCycleCons(
    SCIP *scip,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z,
    vector<Arc> arcs,
    vector<DNode> nodes)
{
    for (size_t k = 0; k < nodes.size(); k++)
    {
        //adding inequality
        ScipCons *cons = new ScipCons(scip, -SCIPinfinity(scip), 0.0, "small-dirCicle cons");
        for (Arc a : arcs)
            cons->addVar(z[a], 1);

        for (size_t i = 0; i < nodes.size(); i++)
        {
            if (i != k)
                cons->addVar(x[nodes[i]], -1);
        }
        cons->commit();
        /* SCIPprintCons(scip, cons->cons, NULL);
        printf("\n"); */
    }
}

/**
 * add cycle removal constraints for cycles of size up to 4
 */
void GLCIPBase::addAllSmallDirectedCycles(
    SCIP *scip,
    GLCIPInstance &instance,
    DNodeSCIPVarMap &x,
    ArcSCIPVarMap &z)
{
    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        for (OutArcIt a(instance.g, v); a != INVALID; ++a)
        {
            //cycles of length two
            DNode u = instance.g.target(a);
            if (findArc(instance.g, u, v) != INVALID)
            {
                vector<DNode> nodes{u, v};
                vector<Arc> arcs{a, findArc(instance.g, u, v)};

                addCycleCons(scip, x, z, arcs, nodes);
            }
            else
            {
                for (OutArcIt b(instance.g, u); b != INVALID; ++b)
                {
                    //cycles of length three
                    DNode w = instance.g.target(b);
                    if (findArc(instance.g, w, v) != INVALID)
                    {
                        vector<DNode> nodes{v, u, w};
                        vector<Arc> arcs{a, b, findArc(instance.g, w, v)};

                        addCycleCons(scip, x, z, arcs, nodes);
                    }
                    else
                    {
                        for (OutArcIt c(instance.g, w); c != INVALID; ++c)
                        {
                            //cycles of length four
                            DNode y = instance.g.target(c);
                            if (findArc(instance.g, y, v) != INVALID)
                            {
                                vector<DNode> nodes{v, u, w, y};
                                vector<Arc> arcs{a, b, c, findArc(instance.g, y, v)};

                                addCycleCons(scip, x, z, arcs, nodes);
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * add arc-influence constraints - a vertex v needs to be active to send influence to w
 */
void GLCIPBase::addLinkingConstraints(SCIP *scip,
                                      GLCIPInstance &instance,
                                      DNodeSCIPVarMap &x,
                                      ArcSCIPVarMap &z)
{
    for (ArcIt a(instance.g); a != INVALID; ++a)
    {
        DNode v = instance.g.source(a);
        DNode w = instance.g.target(a);
        Arc back = findArc(instance.g, w, v);

        if (back == INVALID)
        {
            ScipCons *cons = new ScipCons(scip, 0, SCIPinfinity(scip), "linking cons");

            cons->addVar(x[v], 1);
            cons->addVar(z[a], -1);

            cons->commit();
        }
    }
}

/**
 * the number of activated vertices is at least (alpha * num_vertices) 
 */
void GLCIPBase::addCoverageConstraints(SCIP *scip,
                                       GLCIPInstance &instance,
                                       DNodeSCIPVarMap &x)
{
    ScipCons *covConstraint = new ScipCons(scip, ceil(instance.alpha * instance.n), SCIPinfinity(scip), "coverage cons");

    for (DNodeIt v(instance.g); v != INVALID; ++v)
    {
        covConstraint->addVar(x[v], 1);
    }
    covConstraint->commit();
}

/**
 * Computes the cost paid to activate a vertex v with a given weight of influence
 */
double GLCIPBase::cheapestIncentive(const GLCIPInstance &instance,
                                    const DNode &v,
                                    double exertedInfluence)
{
    double cost = 0;
    // assuming that the incentives are sorted in an increasing order
    // uses the first incentive that overcomes the threshold of v
    for (unsigned int i = 0; i < instance.incentives[v].size(); i++)
    {
        if (exertedInfluence + instance.incentives[v][i] >= instance.threshold[v])
        {
            cost = instance.incentives[v][i];
            break;
        }
    }
    return cost;
}

/**
 * Computes the cost paid to activate a vertex v with a given set of incoming neigobors
 */
double GLCIPBase::costInfluencingSet(const GLCIPInstance &instance,
                                     const DNode &v,
                                     const set<DNode> &nodes)
{
    int thr = instance.threshold[v];
    double cost = 0;

    // activation function
    double exertedInfluence = 0;
    for (DNode u : nodes)
    {
        Arc a = findArc(instance.g, u, v);
        exertedInfluence += instance.influence[a];
    }

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

bool GLCIPBase::intersects(set<DNode> set1, set<DNode> set2)
{
    if (set1.empty() || set2.empty())
        return false;

    for (DNode v : set1)
    {
        if (set2.count(v))
            return true;
    }

    return false;
}
