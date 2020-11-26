import snap
import networkx as nx 
import random
import math

#H = snap.GenForestFire(200, 0.9, 0.1)
H = snap.GenCopyModel(200, .5)

#for v in SG.Nodes():
    #print(v.GetId())

""" for EI in G.Edges():
    print("edge: (%d, %d)" % (EI.GetSrcNId(), EI.GetDstNId())) """

G = nx.DiGraph()
for v in H.Nodes():
    G.add_node(v.GetId())

#add arcs and defining the weight of influence
for e in H.Edges():
    w = random.randint(1, 10)
    G.add_edge(e.GetSrcNId(), e.GetDstNId(), weight = w)
    #G.add_edge(e.GetSrcNId(), e.GetDstNId())

#removing selfloops of G
""" for u,v in list(G.edges()):
    if(u==v):
        G.remove_edge(u,v)
        print("removing selfloop from ", u) """

print("nnodes narcs type")
print(len(G.nodes()), len(G.edges()), "digraph")

#TODO verify if G is acyclic
if(nx.is_directed_acyclic_graph(G)):
    print("Is DAG")
else:
    print("NOT a DAG")
exit(0)

#setting the thresholds on vertices
thr = list()
for v in list(G.nodes()):
    #print(v)
    mu = 0.7 * G.in_degree(v)
    sum = 0
    for s,t,d in G.in_edges(v, data=True):
        sum += d['weight']

    if G.in_degree(v) != 0 :
        sigma = sum/G.in_degree(v)
    else :
        sigma = 0
    
    zeta = random.normalvariate(mu, sigma)
    thr.append( math.floor( max(1, min(zeta, sum) ) ) )
    #print("thr = ", max(1, min(zeta, sum)))
    #print("floor of thr = ", math.floor(max(1, min(zeta, sum))))

pos = nx.circular_layout(G)
posx = list()
posy = list()
for item in pos.items():
    #print(item[1][0], item[1][1])
    posx.append(item[1][0])
    posy.append(item[1][1])

max_thr = 0
for i in range(len(thr)) :
    if thr[i] > max_thr :
        max_thr = thr[i]
#print("max thr = ", max_thr)

#printing the node attributes
print("nodename posx posy threshold")
i = 0
for v in list(G.nodes()):
    #print(v+1, posx[i], posy[i], thr[i], "1,2,3")
    print(v, posx[i], posy[i], thr[i])
    i += 1

print("tail head influence")
for u,v,d in list(G.edges(data=True)):
    print(u, v, d['weight'])