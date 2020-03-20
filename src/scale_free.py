import networkx as nx 
import matplotlib.pyplot as plt 
from operator import itemgetter
import random
import math

#G = nx.directed_configuration_model()
M = nx.scale_free_graph(2000, alpha=0.15, beta=0.7, gamma=0.15, delta_in=0.2,
                     delta_out=0, create_using=None, seed=None)
#M = nx.scale_free_graph(20)
G = nx.DiGraph(M) # to simple digraph

print("nnodes narcs type")
print(len(G.nodes()), len(G.edges()), "digraph")

#random.normalvariate(mu, sigma)

#removing selfloops of G
for u,v in list(G.edges()):
    if(u==v):
        G.remove_edge(u,v)
        #print("removing selfloop from ", u)

#defining the weight of influence on arcs
for u,v,d in G.edges(data=True):
    d['weight'] = random.randint(1, 10)

#setting the thresholds on vertices
thr = list()
for v in list(G.nodes()):
    #print(v)
    mu = 0.7 * G.in_degree(v)
    sum = 0
    for s,t,d in G.in_edges(v, data=True):
        #print(s, t, d["weight"])
        sum += d['weight']
    #print("sum = ", sum)

    if G.in_degree(v) != 0 :
        sigma = sum/G.in_degree(v)
    else :
        sigma = 0
    
    #print("mu = ",mu, "sigma = ", sigma)
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
print("nodename posx posy threshold incentives")
i = 0
for v in list(G.nodes()):
    print(v+1, posx[i], posy[i], thr[i], "1,2,3")
    i += 1

print("tail head influence")
for u,v,d in list(G.edges(data=True)):
    print(u + 1, v + 1, d['weight'])

""" H = nx.DiGraph(G)
print("number of selfloops = ", H.number_of_selfloops())
print("number of edges = ", H.number_of_edges())
for u,v,d in list(H.edges(data=True)):
    print(u, v, d['weight']) """

#nx.draw_circular(G, with_labels = True, alpha=.4)
#nx.draw_spring(G, with_labels=True)
#plt.show()
#plt.savefig("scalefree.png")


#plt.subplot(122)
#nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True)
#plt.show()