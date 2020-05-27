import networkx as nx 
import matplotlib.pyplot as plt 
from operator import itemgetter
import random
import math

#fh = open("../in/realworld/wiki-vote.txt", 'rb')
#G = nx.read_edgelist(fh, nodetype=str, create_using=nx.DiGraph)

#fh = open("../in/realworld/adolescent_health.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

#fh = open("../in/realworld/highschool.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

#fh = open("../in/realworld/innovation.txt", 'rb')
#G = nx.read_edgelist(fh, nodetype=str, create_using=nx.DiGraph)

#fh = open("../in/realworld/residense.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

#fh = open("../in/realworld/advogato.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)
#fh = open("../in/realworld/HepTh.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

#fh = open("../in/realworld/dblp.txt", 'rb')
#G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

fh = open("../in/realworld/cora.txt", 'rb')
G = nx.nx.read_weighted_edgelist(fh, comments='%', nodetype=int, create_using=nx.DiGraph)

fh.close()
labels = [1]
nx.set_node_attributes(G, labels, 'labels')
i = 0
for v in G.nodes():
    G.nodes[v]['label'] = i
    i += 1

#for u, v in G.edges():
#    print(G.nodes[u]['label'], G.nodes[v]['label'])

#removing selfloops of G
for u,v in list(G.edges()):
    if(u==v):
        G.remove_edge(u,v)

print("nnodes narcs type")
print(len(G.nodes()), len(G.edges()), "digraph")

#defining the weight of influence on arcs
for u,v,d in G.edges(data=True):
#    d['weight'] = random.randint(1, 10)
    d['weight'] = 1

#setting the thresholds on vertices
thr = list()
for v in list(G.nodes()):
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

pos = nx.circular_layout(G)
posx = list()
posy = list()
for item in pos.items():
    posx.append(item[1][0])
    posy.append(item[1][1])

max_thr = 0
for i in range(len(thr)) :
    if thr[i] > max_thr :
        max_thr = thr[i]

#printing the node attributes
print("nodename posx posy threshold")
i = 0
for v in list(G.nodes()):
    print(G.nodes[v]['label'], posx[i], posy[i], thr[i])
    i += 1

print("tail head influence")
for u,v,d in list(G.edges(data=True)):
    print(G.nodes[u]['label'], G.nodes[v]['label'], d['weight'])