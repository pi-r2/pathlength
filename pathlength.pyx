# distutils: language = c++
# cython: profile=False

#imports

cimport cython
from libcpp.algorithm cimport sort as csort
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from sage.sets.disjoint_set import DisjointSet
from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx
from cysignals.memory cimport sig_calloc, sig_free


from itertools import combinations
from sage.graphs.distances_all_pairs import diameter
from sage.graphs.distances_all_pairs import distances_all_pairs, eccentricity, floyd_warshall
from math import log, factorial
from sage.graphs.graph import Graph
from sage.graphs.graph_generators_pyx import RandomGNP
import matplotlib.pyplot as plt
from random import shuffle
import time


#Utilitaires

 



cdef class UnionFind() :
    #A class that implements union find in cython


    def __init__(self,set l, int n) : #set l is a subset of [1:n] (allows us to work with arrays instead of dicts)
        self.n = n 
        self.rank = <int *>sig_calloc(n, sizeof(int)) 
        self.dico = <int *>sig_calloc(n, sizeof(int))
        self.elem = l
        cdef int i
        for i in l : 
            self.dico[i] = i
            self.rank[i] = 0 


    cdef void free(self) : #you need to clean manually !
        sig_free(self.rank)
        sig_free(self.dico)
    
    cdef int find(self, int i) :
        if self.dico[i] != i :
            self.dico[i] = self.find(self.dico[i])
        return self.dico[i]


    cdef void union(self, int i, int j) :
        
        cdef int ri 
        cdef int rj 
        
        ri = self.find(i) 
        rj = self.find(j)

        if ri!= rj :
            if self.rank[ri] < self.rank[rj] :
                self.dico[i] = rj 
                self.rank[rj] += 1
                
            else :
                self.dico[j] = ri 
                self.rank[ri] += 1

    cdef set set_(self, int i) :
        cdef int ri 
        ri = self.find(i)
        ret = set()
        cdef int j 
        for j in self.elem :
            if ri == self.find(j) :
                ret.append(j)
        return ret

    cdef list sets(self) : 
        cdef dict temp = dict() 
        
        cdef int i
        cdef int ri 
        
        for i in self.elem :
            ri = self.find(i)
            if ri not in temp:
                temp[ri] = set()
            temp[ri].add(i)

        
        
        return list(temp.values())



def layout_to_decomp(g,L) : #sage graph g, L permutation of vertices. Returns the bags induced by the layout L
    bags = []
    for i in range(len(L))  :
        v = L[i]
        Bv = [v]
        for u in L[:i] :
            for w in g.neighbors(u) : 
                if w not in L[:i] :
                    Bv.append(u)
                    break
        bags.append(Bv)
        
    return bags

def length2(g, L) : #Given layout L on a graph G, return length of the decomposition induced by the layout
    distance = distances_all_pairs(g) 
    path = layout_to_decomp(g,L)
    diam= 0
    for bag in path :
        for u in bag :
            for v in bag :
                if distance[u][v] > diam :
                    diam = distance[u][v]
    return diam



def is_path_decomposition(g, path) : #checks wether path is a valid path decomposition of g. path is a list of lists of vertices
    V = set(g.vertices())
    X = set([])
    #i) union is all vertices
    for bag in path :
        X = X.union(set(bag))
    if X.cardinality() != V.cardinality() :
        print("rule 1")
        return False
    #ii) each edge appears in a bag 
    flag = False
    for u,v in g.edges(labels=False) :
        for bag in path :
            if (u in bag) and (v in bag) :
                flag = True
                break
    if not flag :
        print("rule 2")
        return False
    
    #iii) for all i<=j<=k, Xi inter Xk included in Xj
    X= map(set,path)
    for i in range(len(X)) :
        for j in range(i,len(X)) :
            for k in range(j, len(X)) :
                if not (X[i].intersection(X[k])).issubset(X[j]) :
                    print("rule 3")
                    print(i,j,k)
                    return False
    return True
    

def decomp_to_layout(g,bags) : #Given a path decomposition of g, returns the induced layout.
    r = len(bags)

    L = list(bags[0])
    V = set(bags[0])
    for i in range(1,r) :
        Y = list(set(bags[i]).difference(V))
        while Y!=[] :
            L.append(Y.pop())
            V.add(L[len(L)-1])
    return L





def longest_isometric_cycle(G):  #longest isometric cycle detection by Lohkstanov, original implementation by Theo Qui

    """
    computes the longest isometric cycle of G
    Algorithm by Daniel Lokshtanov
    """


    def power(G, k): #A remplacer par l'algorithme de folklore du calcul de puissances de graphe
        """
        computes G power k
        """
        """
        G_pow_k = Graph()
        G_pow_k.add_vertices(G.vertices())
        """
        dist = G.distance_all_pairs() #on doit pouvoir faire mieux

        edges  = {}

        for u in G :
            for v in G :
                dist_u_v = dist[u].get(v, None)
                if not u == v and not dist_u_v == None and dist_u_v <= k: # if v isn't reachable from u, we say the distance is one more than the number of edges (~ +inf)
                    #G_pow_k.add_edge(u, v)
                    edges[(u,v)] = 1
                else :
                    edges[(u,v)] = 0

        return edges

    #### the actual algorithm starts here
    ans = 0
    if G.is_tree():
        return 0
    dist = G.distance_all_pairs()
    star = (-1,-1)

    for k in range(3, G.order()+1):
        Vk = set()
        for u in G :
            for v in G :
                if dist[u][v] == k // 2:
                    Vk.add((u, v)) #O(n**2) => OK

        Ek = set()
        for u, v in Vk:
            for w, x in Vk:
                if dist[u][w] == 1 and dist[x][v] == 1:
                    Ek.add(((u,v), (w,x)))  #O(Vk**2) => OK

        Gk = Graph([Vk, Ek], format='vertices_and_edges')

        Gk_pow_k = power(Gk, k//2)

        for u, v in Vk :
            for x in G :
                if not k%2 :
                    if u == x and Gk_pow_k[((u, v), (v, x))] :
                        ans = k
                else :
                    if (dist[x][u] == 1 and (v, x) in Vk) and Gk_pow_k[((u, v), (v, x))]:
                        ans = k
    return ans






        
        
#2-approximation algorithm by **Dragan, Kohler, Leitert**    
def two_approximation(g) : #returns a 2-approximation of the pathlength of g
    d = distances_all_pairs(g) 
    diam = diameter(g) #these two operations can be fused 
    best_bags = []
    for s in g.vertices() :
        current_diam = 0
        L = [[] for _ in g.vertices()]
        Lplus = [[] for _ in g.vertices()]
        Lplus[0]= [s]
        for v in g.vertices() :
            L[d[s][v]].append(v)
        for i in range(1,len(L)) :
            Lplus[i] = list(L[i])
            for v in L[i-1] :
                if set(g.neighbors(v)).intersection(set(L[i])) :
                    Lplus[i].append(v)
            for u in Lplus[i] :
                for v in Lplus[i] :
                    if d[u][v] > current_diam :
                        current_diam = d[u][v]
        if current_diam <= diam :
            diam = current_diam
            best_bags = Lplus

    while best_bags[len(best_bags)-1] == [] :
        best_bags.pop() 
    return best_bags, diam
           
#lower bound by detection of isometric spider graph



def lower_bound2(g) : #a homemade lower bound by spider-graph detection - it is not exhaustive

    d, sps = floyd_warshall(g,paths=True, distances=True)
    ecc = eccentricity(g)
    s = 0
    t = 0
    w = 0
    diam = 1
    m2 = 1
    for i in g :
        for j in g :
            if d[i][j] > diam :
                    diam = d[i][j]
                    s = i
                    t= j
    #s,t achieves a diameter

    diameters = set([(s,t)])
    for i in g : #check all diameters, not just one.
        for j in g :
            if d[i][j] == diam :
                    diameters=diameters.union(set([(i,j)]))

    LB = 1
    star = (0,0,0,-1)
    for s,t in diameters :
        sp_st = [s]
        while sp_st[len(sp_st)-1] != t :
            sp_st.append(sps[t][sp_st[len(sp_st)-1]])
        for v in sp_st[LB:-LB] :
            L = []
            for w in g.vertices() :
                if (d[s][w] == d[s][v] + d[v][w]) and (d[t][w] == d[t][v] + d[v][w]) and (LB < min(d[s][v],d[v][w],d[v][t])) :
                    LB = min(d[s][v],d[v][w],d[v][t])
                    star = (s,t,w,v)
            if d[s][v] < LB or d[v][t] < LB :
                break


    return LB







@cython.boundscheck(False)
@cython.wraparound(False)
cdef vector[int] * neighbors_from_adj(int n, unsigned short ** adj) :
    
    cdef vector[int] * neighbors
    neighbors = <vector[int] *> sig_calloc(n, sizeof(vector[int])) 



    cdef int i 
    cdef int j 
    for i in range(n) :
        for j in range(i) :
            if adj[i][j] == 1 : #allows adj to be a distance matrix since im lazy
                neighbors[i].push_back(j)
                neighbors[j].push_back(i)
    return neighbors


@cython.boundscheck(False)
@cython.wraparound(False)
cdef vector[pair[int,int]] vflatten( vector[vector[pair[int,int]]] L ) : #vector flatten : 2D to 1D 
    cdef vector[pair[int,int]] L2 
    cdef int n 
    n = L.size()
    cdef int i 
    cdef int j
    for i in range(n) :
        for j in range(L[i].size()) :
            L2.push_back(L[i][j])
    return L2 


        
#connected components tools
@cython.boundscheck(False)
@cython.wraparound(False)
cdef list connectedComponentsUF(int n, vector[int] * neighbors, set V_S) : #find connected components of G\S with UF data structure
    cdef UnionFind D
    D = UnionFind(V_S, n)

    cdef int u 
    cdef int v 

    for u in V_S :
        for v in neighbors[u] :
            if v in V_S :
                D.union(u,v)
    
    cdef list ret = D.sets()
    D.free()
    return ret


##Extended greedy Branch and Bound
@cython.boundscheck(False)
@cython.wraparound(False)
cdef vector[int] gBAB2(
        vector[int] * neighbors,
        unsigned short ** distances, 
        int n,
        vector[int] P, 
        set S,
        set V_S, 
        vector[int] P_star, 
        int * UB, 
        int LB, 
        int D, 
        long * coupes, 
        long * count, 
        set prefix_storage,
        int max_prefix_length, 
        int max_prefix_number,
        list cc,
        bint jump_flag,
        int * jumps,
        float t0
        ) :#,
        #bitset_t b_prefix
        #) :
    
    #Implements a step of gBAB for a given prefix
    
    cdef float t 
    t = time.time()

    if t - t0 > 600 :
        return P_star
    count[0]+=1

    if UB[0] == LB :
        return P_star
    

    cdef int r
    r = P.size()
    
    
    if r == n : #The layout is complete
        coupes[0]+=1
        if D < UB[0] :
            UB[0]= D 
            return P
        else : 
            return P_star


    #cdef vector[int] list_V_S 





    
  
    #make sure this prefix is worth considering
    cdef frozenset frozen_prefix
    cdef list fp = []
    frozen_prefix= frozenset()

    cdef int i 

    if r <= max_prefix_length:
        for i in range(r) :
            fp.append(P[i])
        frozen_prefix = frozenset(fp)
        if frozen_prefix in prefix_storage:
            return P_star    

    
    #Determine the list of candidates : it is faster than the extension lemma computation
    #and can be used to prune a branch
    cdef vector[vector[pair[int,int]]] L = vector[vector[pair[int,int]]] (n) #List of potential candidates ordered by increasing new bag diameter

    cdef list F  #border of S : any new bag is exactly F \cup {v}
    F = []
    cdef int diam_F 
    cdef bint flag 
    cdef int d 
    
    #Currently, if bag diameter is equal, they are then ordered by labelling. This depends on the initial ordering
    #and therefore is something that could be enhanced.
    
    diam_F = 0

    cdef int v =  0
    cdef int u =  0
    cdef int w =  0
    cdef list temp

    cdef set prefix_neighbourhood = set()


    for u in P : 
        flag = False
        for w in neighbors[u] : 
            prefix_neighbourhood.add(u)
            if w in V_S : #find the frontier of S : the set of vertices in S that have a neighbour in $V \setminus S$
                F.append(u)
                flag = True
                break 


    for u in F :#naive diam calculation of the frontier
        for v in F :
            if distances[u][v] > diam_F :
                diam_F = distances[u][v]

    if diam_F > UB[0] : #The frontier is contained in any new bags, so its diameter is a lower bound on bag diameter 
        return P_star

    cdef vector[vector[pair[int,int]]] L_prio = vector[vector[pair[int,int]]] (n) #List of potential candidates ordered by increasing new bag diameter




    cdef bint connected_pl = True #Change to false to remove connected heuristic

    for v in V_S :
        d=0
        #diam(bag_v) = max(diam F, max(dist v, s\in F))
        for u in F :
            if d > UB[0] :
                break
            if distances[u][v] > d : 
                d = distances[u][v]
        if d < UB[0] :#insertion baquet
            if connected_pl :
                if v in prefix_neighbourhood :
                    L_prio[d].push_back(pair[int, int](d,v)) 
                else :
#                    continue #COMMENT IF YOU WANT TRUE PATHLENGTH
                    L[d].push_back(pair[int, int](d,v)) 
            else :
                L_prio[d].push_back(pair[int, int](d,v)) 
                
                
    cdef vector[pair[int,int]] L2_prio
    L2_prio = vflatten(L_prio)
    cdef vector[pair[int,int]] L2 
    L2 = vflatten(L)
    
    for i in range(L2.size()):
        L2_prio.push_back(L2[i])
    


    if not L2_prio.size() : #this branch has no interesting leaves
        return P_star
    #else : there might be better candidates in the form of connected components
    

    #Extension Lemma 1: 
    
    # if c is a set on G\S such that 
    # - N(c) is included in S \cup c ( iff c is a connected component of G\S)
    # - diameter(SUC) = diameter(S)
    #Then there exists an optimal layout prefixed by P that is
    #also prefixd by P.C
    


   
    cdef set c 
    cdef list new_sets
    
    if not jump_flag :     #we know precisely which connected component is subject to having been split (the one containing the last vertex in the prefix)
                           #else we know that we already have removed the previous component : nothing to compute 
        u = P.back()
        for i in range(len(cc)) :
            if u in cc[i] :
                c = cc.pop(i)
                break 
        if len(c) > 1 :
            c.discard(u)
            new_sets = connectedComponentsUF(n, neighbors, c) 
            new_sets.extend(cc)
            cc = new_sets                
        #cc.sort(key = len, reverse = True) #How useful is sorting ?    
    
    
    cdef vector[int] P_prime   

    cdef set ens 
    cdef set set_F
    set_F = set(F)
    cdef list temp_cc

    for i in range(len(cc)) : #note pour le futur : on ne peut pas se passer du nonecheck produit par cython, mais étant "unlikely_", il n'est en pratique jamais considéré : donc pas d'overhead.  
        c = cc[i]
        ens = set_F.union(c)

        d=0

        for u,w in combinations(ens, 2) :
            if distances[u][w] > d : 
                d = distances[u][w]
                if d > D :
                    break
                    
        if d <= D :

            for v in c :
                P.push_back(v)
                S.add(v)
                V_S.remove(v)

            jumps[0]+=1
            cc.pop(i)
            return gBAB2(neighbors, distances, n, P, S, V_S, P_star, UB, LB, D, coupes, count, prefix_storage, max_prefix_length, max_prefix_number, cc, True, jumps, t0)#, b_prefix)
    

    #Copies to be used as arguments
    cdef set S2 
    cdef set V_S2

    cdef list cc_copy
    cdef vector[int] Pv
    for i in range(L2_prio.size()) : #Do not explore branches that cannot be improvements to the current upper bound
        
        
        d = L2_prio[i].first
        v = L2_prio[i].second #Cython unpacking of cpp pairs is not optimized (yet) : we cannot use d,v = pair  (see https://github.com/cython/cython/issues/1429)
        if UB[0] == LB :
            return P_star

        if d < UB[0]  :
            
            Pv.clear()
            for u in P :
                Pv.push_back(u)
            Pv.push_back(v)
            S2 = set(S)
            V_S2 = set(V_S)

            S2.add(v)
            V_S2.remove(v)
            cc_copy = [set(s) for s in cc]
            P_star = gBAB2(neighbors, distances, n, Pv, S2, V_S2, P_star, UB, max(LB,d), max(d,D), coupes, count, prefix_storage, max_prefix_length, max_prefix_number, cc_copy, False, jumps, t0)#, b_prefix)
        else : 
            continue
            
        #prefix storage update :
    if D < UB[0] and r <= max_prefix_length and len(prefix_storage) < max_prefix_number : #after having explored the branch
        prefix_storage.add(frozen_prefix)
            
    return P_star
        
        
cdef list pathlength_gBAB2(GenericGraph_pyx  g,
                          int max_prefix_length=20,
                          int max_prefix_number=10**6,
                          bint verbose = True
                          ) : 
    #allocs

    cdef int LB1 
    cdef int LB2 
    cdef int LB  
    cdef unsigned int n 
    cdef long coupes 
    cdef long count     
    cdef int UB 
    cdef int K

    cdef vector[int] P_star 
    cdef vector[int] P
    cdef list L = []
    cdef list ecc
    cdef list bags 

    cdef unsigned short * c_distances
    cdef set prefix_storage

    cdef vector[int] * neighbors
    cdef set S 
    cdef set V_S


    
    #precomputations 

    n = g.order()
    bags, UB = two_approximation(g)
    LB1 = lower_bound2(g)
    LB2 = longest_isometric_cycle(g)//2
    LB = max(LB1, LB2, UB//2)
    P_star = decomp_to_layout(g,bags)
    ecc = eccentricity(g)
    
    c_distances = c_distances_all_pairs(g)
    cdef unsigned short ** distances 
    distances = <unsigned short **>sig_calloc(n, sizeof(unsigned short *))
    

    cdef unsigned short i
    for i in range(n):
        distances[i] = c_distances + i * n
    neighbors = neighbors_from_adj(n, distances)

    def ec(int i) :
        return ecc[i]

    L = sorted(list(range(n)), key = ec, reverse = True) #placeholder
    coupes = 0
    count = 0
    prefix_storage = set()
    
    cdef vector[int] Pv

    K = UB 
    
    if UB == LB : #If we're lucky, seize it
        return [UB, P_star, coupes, count, LB]
    
    cdef int branches 
    branches = 0

    
    experimental = False #Change to true to break after exploring diameters

    cdef int jumps

    t0 = time.time()
    for v in L:
        if UB == LB :
            break 
        if experimental :
            if ecc[v] < ecc[L[0]] :
                print("exp")
                break
        branches +=1
    
        S = set([v])
        V_S = set(range(n))
        V_S.remove(v)
        Pv = gBAB2(neighbors, distances,n, [v], S, V_S, P_star, &UB, LB, 0,&coupes, &count, prefix_storage, max_prefix_length, max_prefix_number, [set(range(n))], False, &jumps, t0)#, b_prefix)
        if UB <= K :
            K = UB
            P_star = Pv
        if verbose :
            print(str(branches)+"/"+str(n)+" branches explorées")
    
    
    t = time.time()
    #frees

    sig_free(c_distances)
    sig_free(distances)
    sig_free(neighbors)
    
    if verbose:
        print("eccentricities :", ecc)
        print("jumps :", jumps)

    if t - t0 > 600 : 
        return [-1, P_star, coupes, count, LB]    
    return [K, P_star, coupes, count, LB] #return a list containing the pathlength, the optimal layout, the number of cuts and explored leaves aswell as the chosen lower bound.




cpdef run(g, V = None) :


    ##the main function for practical use. Takes a graph and an optional 'permutation' argument.
    ##If no permutation is given, the vertices will be shuffled before the graph is processed.

    ##You can pass V = g.vertices() if you wish not to shuffle the graph.



    #avoid preliminary assumptions on vertex order and remove weird labels 
    g.relabel(inplace=True, perm = list(range(g.order())))

    if V == None :
        V = list(g)
        shuffle(V)
    print("perm = ", V)
    g.relabel(inplace = True, perm=V)
    
    rev = [0 for _ in V]
    for i in range(len(V)) : #permutation inverse
        rev[V[i]] = i
    
    p=  g.plot( layout = 'spring')
    p.show()
    
    
    pl, L, c, c2, LB  = pathlength_gBAB2(g)
    
    app = two_approximation(g)[1]
    print("lower bound : ", LB)
    print("two approx : ", app)
    print("number of layouts explored : ", c)
    print("Magnitude of total number of layouts : 10^" + str(round(log(factorial(len(g.vertices())),10))))
    print('pathlength : ', pl)
    print("diameter : ", g.diameter())
    print('layout : ', L)
    print("calls to BAB : ", c2)
    print('length verification')
    print(length2(g, L))
    
    g.relabel(inplace = True, perm=rev)


cpdef timed_run(g) : #a function for benchmarking - no other practical use.
    #avoid preliminary assumptions on vertex order and remove weird labels

    g.relabel(inplace=True, perm = list(range(g.order())))



    #V = list(g)
    #shuffle(V)
    #g.relabel(inplace = True, perm=V)
    t0 = time.time()
    pl, L, c,c2, LB = pathlength_gBAB2(g, max_prefix_length=40, max_prefix_number=10**7, verbose = False)
    dt = time.time() - t0



    UB = two_approximation(g)[1]
    return g.order(), g.size(), dt, pl, LB, UB, g.diameter()












##Simple BaB for comparison purposes



##Extended greedy Branch and Bound
@cython.boundscheck(False)
@cython.wraparound(False)
cdef vector[int] gBAB(
        vector[int] * neighbors,
        unsigned short ** distances, 
        int n,
        vector[int] P, 
        set S,
        set V_S, 
        vector[int] P_star, 
        int * UB, 
        int LB, 
        int D, 
        long * coupes, 
        long * count, 
        set prefix_storage,
        int max_prefix_length, 
        int max_prefix_number,
        float t0
        ) :
    
    #Implements a step of gBAB for a given prefix
    
    cdef float t 
    t = time.time()

    if t - t0 > 600 :
        return P_star
    count[0]+=1

    if UB[0] == LB :
        return P_star
    

    cdef int r
    r = P.size()
    
    
    if r == n : #The layout is complete
        coupes[0]+=1
        if D < UB[0] :
            UB[0]= D 
            return P
        else : 
            return P_star


    #cdef vector[int] list_V_S 





    
  
    #make sure this prefix is worth considering
    cdef frozenset frozen_prefix
    cdef list fp = []
    frozen_prefix= frozenset()

    cdef int i 

    if r <= max_prefix_length:
        for i in range(r) :
            fp.append(P[i])
        frozen_prefix = frozenset(fp)
        if frozen_prefix in prefix_storage:
            return P_star    

    
    #Determine the list of candidates : it is faster than the extension lemma computation
    #and can be used to prune a branch
    cdef vector[vector[pair[int,int]]] L = vector[vector[pair[int,int]]] (n) #List of potential candidates ordered by increasing new bag diameter

    cdef list F  #border of S : any new bag is exactly F \cup {v}
    F = []
    cdef int diam_F 
    cdef bint flag 
    cdef int d 
    
    #Currently, if bag diameter is equal, they are then ordered by labelling. This depends on the initial ordering
    #and therefore is something that could be enhanced.
    
    diam_F = 0

    cdef int v =  0
    cdef int u =  0
    cdef int w =  0
    cdef list temp

    cdef set prefix_neighbourhood = set()


    for u in P : 
        flag = False
        for w in neighbors[u] : 
            prefix_neighbourhood.add(u)
            if w in V_S : #find the frontier of S : the set of vertices in S that have a neighbour in $V \setminus S$
                F.append(u)
                flag = True
                break 


    for u in F :#naive diam calculation of the frontier
        for v in F :
            if distances[u][v] > diam_F :
                diam_F = distances[u][v]

    if diam_F > UB[0] : #The frontier is contained in any new bags, so its diameter is a lower bound on bag diameter 
        return P_star

    cdef vector[vector[pair[int,int]]] L_prio = vector[vector[pair[int,int]]] (n) #List of potential candidates ordered by increasing new bag diameter




    cdef bint connected_pl = True 

    for v in V_S :
        d=0
        #diam(bag_v) = max(diam F, max(dist v, s\in F))
        for u in F :
            if d > UB[0] :
                break
            if distances[u][v] > d : 
                d = distances[u][v]
        if d < UB[0] :#insertion baquet
            if connected_pl :
                if v in prefix_neighbourhood :
                    L_prio[d].push_back(pair[int, int](d,v)) 
                else :
                #    continue #COMMENT IF YOU WANT TRUE PATHLENGTH
                    L[d].push_back(pair[int, int](d,v)) 
            else :
                L_prio[d].push_back(pair[int, int](d,v)) 
                
                
    cdef vector[pair[int,int]] L2_prio
    L2_prio = vflatten(L_prio)
    cdef vector[pair[int,int]] L2 
    L2 = vflatten(L)
    
    for i in range(L2.size()):
        L2_prio.push_back(L2[i])
    


    if not L2_prio.size() : #this branch has no interesting leaves
        return P_star
 
    #temp
    cdef set S2 
    cdef set V_S2

    cdef vector[int] Pv
    for i in range(L2_prio.size()) : #Do not explore branches that cannot be improvements to the current upper bound
        
        
        d = L2_prio[i].first
        v = L2_prio[i].second #Cython unpacking of cpp pairs is not optimized (yet) : we cannot use d,v = pair  (see https://github.com/cython/cython/issues/1429)
        if UB[0] == LB :
            return P_star

        if d < UB[0]  :
            
            Pv.clear()
            for u in P :
                Pv.push_back(u)
            Pv.push_back(v)
            S2 = set(S)
            V_S2 = set(V_S)

            S2.add(v)
            V_S2.remove(v)
            P_star = gBAB(neighbors, distances, n, Pv, S2, V_S2, P_star, UB, max(LB,d), max(d,D), coupes, count, prefix_storage, max_prefix_length, max_prefix_number, t0)#, b_prefix)
        else : 
            continue
            
        #prefix storage update :
    if D < UB[0] and r <= max_prefix_length and len(prefix_storage) < max_prefix_number : #after having explored the branch
        prefix_storage.add(frozen_prefix)
            
    return P_star







           
cdef list pathlength_gBAB(GenericGraph_pyx  g,
                          int max_prefix_length=20,
                          int max_prefix_number=10**6,
                          bint verbose = True
                          ) : 
    #allocs

    cdef int LB1 
    cdef int LB2 
    cdef int LB  
    cdef unsigned int n 
    cdef long coupes 
    cdef long count     
    cdef int UB 
    cdef int K

    cdef vector[int] P_star 
    cdef vector[int] P
    cdef list L = []
    cdef list ecc
    cdef list bags 

    cdef unsigned short * c_distances
    cdef set prefix_storage

    cdef vector[int] * neighbors
    cdef set S 
    cdef set V_S


    
    #precomputations 

    n = g.order()
    bags, UB = two_approximation(g)
    LB1 = lower_bound2(g)
    LB2 = longest_isometric_cycle(g)//2
    LB = max(LB1, LB2, UB//2)
    P_star = decomp_to_layout(g,bags)
    ecc = eccentricity(g)
    
    c_distances = c_distances_all_pairs(g)

    #bitsets for elements of the prefix 
    
    #cdef bitset_s * b_prefix
    #bitset_init(b_prefix, n)

    
    

    cdef unsigned short ** distances 
    distances = <unsigned short **>sig_calloc(n, sizeof(unsigned short *))
    

    cdef unsigned short i
    for i in range(n):
        distances[i] = c_distances + i * n
    neighbors = neighbors_from_adj(n, distances)

    def ec(int i) :
        return ecc[i]

    L = sorted(list(range(n)), key = ec, reverse = True) #placeholder
    coupes = 0
    count = 0
    prefix_storage = set()
    
    cdef vector[int] Pv

    K = UB 
    
    if UB == LB : #If we're lucky, seize it
        return [UB, P_star, coupes, count, LB]
    
    cdef int branches 
    branches = 0

    
    experimental = False 

    cdef int jumps

    t0 = time.time()
    for v in L:
        if UB == LB :
            break 
        if experimental :
            if ecc[v] < ecc[L[0]] :
                print("exp")
                break
        branches +=1
    
        S = set([v])
        V_S = set(range(n))
        V_S.remove(v)
        Pv = gBAB(neighbors, distances,n, [v], S, V_S, P_star, &UB, LB, 0,&coupes, &count, prefix_storage, max_prefix_length, max_prefix_number, t0)#, b_prefix)
        if UB <= K :
            K = UB
            P_star = Pv
        if verbose :
            print(str(branches)+"/"+str(n)+" branches explorées")
    
    
    t = time.time()
    #frees

    sig_free(c_distances)
    sig_free(distances)
    sig_free(neighbors)
    
    if verbose:
        print("eccentricities :", ecc)
        print("jumps :", jumps)

    if t - t0 > 600 : 
        return [-1, P_star, coupes, count, LB]    
    return [K, P_star, coupes, count, LB] #return a list containing the pathlength, the optimal layout, the number of cuts and explored leaves aswell as the chosen lower bound.





cpdef timed_run_naive(g) : #a function for benchmarking - no other practical use.
    #avoid preliminary assumptions on vertex order and remove weird labels

    g.relabel(inplace=True, perm = list(range(g.order())))



    #V = list(g)
    #shuffle(V)
    #g.relabel(inplace = True, perm=V)


    t0 = time.time()
    pl, L, c,c2, LB = pathlength_gBAB2(g, max_prefix_length=40, max_prefix_number=10**7, verbose = False)
    dt = time.time() - t0



    UB = two_approximation(g)[1]
    return g.order(), g.size(), dt, pl, LB, UB, g.diameter()

