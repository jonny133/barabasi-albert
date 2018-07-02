# -*- coding: utf-8 -*-
"""
Barabasi-Albert network model
@author: jonny133
"""

from random import choice, random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('ggplot')

def draw(G):
    plt.figure()
    nx.draw_networkx(G, pos=nx.spring_layout(G))

def BA(N=10**3, m=3, G0=None):
    
    #################################################
#    Generate graph G/check G0 is acceptable

    if G0 == None:          
        if m == 1 or N < m:
            raise Exception("Invalid parameters for producing initial complete graph, check m!=1")#T>=m and m!=1")
        else:
            G = nx.complete_graph(m) #complete graph m0=m instead of arbitrarily linked..
#        "After t timesteps the Barabási-Albert model generates a network with N = t + m0 nodes and m0 + mt links."

    else:
        if len(G0.nodes()) < 2:
            raise Exception("Invalid graph G0 specified (too few vertices)")
        else:
            G = G0
    
    #############################################
    # Use algorithm to grow network
    
    links = [i for i, j in G.degree for z in xrange(j)] #list of existing vertices repeated according to degree

    size = G.order() #find current size of network to start adding

    for t in xrange(size, N):
#        newnode = size + t

        targets = set()
        while len(targets) < m: #choose m random nodes
            x = choice(links)
            #print "target chosen is", x
            targets.add(x)       
#            targets = random.sample(links,m) 

        G.add_edges_from((t, target) for target in targets)
        
        links.extend([t]*m) #extend links with new node m times and the nodes it connected to
        links.extend(targets)
        
    return G
    
def rando(N = 10**3, m = 3, G0=None): 
    
    #################################################
#    Generate graph G/check G0 is acceptable

    if G0 == None:          
        if m == 1  or N < m:# or  T < m :
            raise Exception("Invalid parameters for producing initial complete graph")#T>=m and m!=1")
        else:
            G = nx.complete_graph(m) #complete graph m0=m instead of arbitrarily linked..
#        "After t timesteps the Barabási-Albert model generates a network with N = t + m0 nodes and m0 + mt links."
#         num = T#-m
    else:
        if len(G0.nodes()) < 2:
            raise Exception("Invalid graph G0 specified (too few vertices)")
        else:
            G = G0
#         num = T 
    
    #############################################
    # Use algorithm to grow network
    
    size = G.order() #find current size of network to start adding

    for t in xrange(size, N):

        targets = set()
        while len(targets) < m: #choose m random nodes
            newTarget = np.random.randint(t)
            #print "new random target", newTarget
            targets.add(newTarget)        

        G.add_edges_from((t, target) for target in targets)
        
    return G
    
def mixed(N = 10**3, m = 3, q=0.5, G0=None): 
    
    #################################################
#    Generate graph G/check G0 is acceptable

    if G0 == None:          
        if m == 1  or N < m:# or  T < m :
            raise Exception("Invalid parameters for producing initial complete graph")#T>=m and m!=1")
        else:
            G = nx.complete_graph(m) #complete graph m0=m instead of arbitrarily linked..
#        "After t timesteps the Barabási-Albert model generates a network with N = t + m0 nodes and m0 + mt links."
#         num = T#-m
    else:
        if len(G0.nodes()) < 2:
            raise Exception("Invalid graph G0 specified (too few vertices)")
        else:
            G = G0
#         num = T 
    
    #############################################
    # Use algorithm to grow network
    links = [i for i, j in G.degree_iter() for z in xrange(j)] #list of existing vertices repeated according to degree
    size = G.order() #find current size of network to start adding
    
    for t in xrange(size, N):
        targets = set()
        while len(targets) < m: #choose m random nodes
            if random()<q:
                x = choice(links)
#                print x, "chosen preferentially"
                targets.add(x)       
            else:
                y = np.random.randint(t)
#                print y, "chosen randomly"
                targets.add(y)
#        print targets
    
        G.add_edges_from((t, target) for target in targets)
        
        links.extend([t]*m) #extend links with new node m times and the nodes it connected to
        links.extend(targets)
    return G

def theopref(x, m, q = None):
    return [float(2*m*(m+1))/(k*(k+1)*(k+2)) for k in x]    
def theorand(x, m,  q = None):
    return [(float(m)/(m+1))**(k-m)/(m+1) for k in x]   
def theomixed(x, m, q = 0.5):
    return [(2/q)*(m*(2-q)/q)**(2/q)*(k+2*m*(1-q)/q)**(-1-2/q) for k in x]

theodists = {
    'pref': theopref, 
    'random': theorand, 
    'mixed': theomixed
    }

def dorepeats(N, m, q=0.5, alg='pref', repeats=20):
    choices = {'pref': BA, 'random': rando, 'mixed': mixed}
    degs = [] # repeated lists all degrees e.g. [3,3,5,3,14,3]
    hists =[] # lists of histogram of degs e.g. [0, 0, 0, 851, 412, 272]
    k1s = [] # list of maximum degree k_1

    for i in xrange(repeats):
        G = choices[alg](N, m, G0=None)      
        degree = G.degree().values()
        hist = nx.degree_histogram(G)
        
        degs.extend(degree)
        hists.append(np.array(hist))
        k1s.append(max(degree))       

    maxlenhist = max(len(i) for i in hists)
    for i in hists:
        i.resize(maxlenhist, refcheck=False)
    ndist = np.sum(hists, axis = 0)
    
    return degs, ndist, k1s
    