from pathlength import timed_run, timed_run_naive, run, two_approximation, lower_bound2, longest_isometric_cycle
import matplotlib.pyplot as plt
import time


def randomOuterPlanar(n) : #Build a random outerplanar graph by cutting a cycle.
    g = graphs.CycleGraph(n)
    g.name("Outerplanar")
    def subdivide(i,j) :
        if i == j :
            return []
        if i == j-1 :
            if randint(i,j) == j :
                return [(i,j)]
            else :
                return []
        """
        a = randint(i,j)
        b = randint(i,j)
        while b == a :
            b = randint(i,j)
        j = max(a,b)
        i = min(a,b)
        """
        
        coinflip= randint(0,1)
        if coinflip : 
            a = randint(i+1,j)
        else :
            a = randint(i,j-1)
        
        #graphs are too dense,add a coinflip to ignore edge
        coinflip= randint(0,1)
        if coinflip : 
            return [(i,j)] + subdivide(i,a) + subdivide(a,j)
        else :
            return subdivide(i,a) + subdivide(a,j)
            
    i = randint(0,n-1)
    j = randint(0, n-1)
    edges = subdivide(0, min(i,j)) + subdivide(min(i,j),max(i,j)) + subdivide(max(i,j),n-1)
    g.add_edges(edges)
    
    return g

def randomConnectedGNP(n) : #randomly generate a GNP graph with high hopes of it being connex, and select the largest connected component.
    def p(n) :
        return 2*log(n)/n
    g = RandomGNP(n,p(n))
    g2 = g.connected_components_subgraphs()[0]
    g2.relabel(inplace=True, perm = list(range(g2.order())))
    return g2


import networkx
"""
g= Graph(networkx.read_edgelist("./graphs/edgelist_rome/grafo6353.39.edgelist"))
print(g.vertices())
print(g.edges())

p =g.plot(layout = "spring")
p.show()
"""
def test(n_max, n_min = 5,generator="gnp") :
    #generator is either "gnp" or "op" 
    if generator == "gnp" : 
        gen = randomConnectedGNP
    else :
        gen = randomOuterPlanar

    for i in range(n_min,n_max) :
        print(i)
        for _ in range(10) :
            g = gen(i)
            perm = shuffle(list(range(i)))
            g.relabel(inplace=True, perm = perm)
            with open("c_results_"+generator, "a") as results :

          
                n, m, temps, pl, LB, UB, diam = timed_run(g)
                results.write(str(n) + " " + str(m) + " " + str(temps) + " " + str(pl) + " " + str(LB) + " " + str(UB) + " " + str(diam)+'\n')

                with open("c_results_naive_"+generator, "a") as results2 :
                    _, _, temps2, _, _, _, _ = timed_run_naive(g)
                    results2.write(str(n) + " " + str(m) + " " + str(temps2) + " " + str(pl) + " " + str(LB) + " " + str(UB) + " " + str(diam)+'\n')
            

            print(str(n) + " " + str(m) + " " + str(temps) + " " + str(pl) + " " + str(LB) + " " + str(UB) + " " + str(diam)+'\n')

#test(50,5,'op')



def test_precomp(path_to_lib, path_to_results) :
    l = os.listdir(path_to_lib)
    l2 = []
    for f in l : 
        f2 = f.split(".")
        l2.append((int(f2[1]),f2[0],f2[2]))
    l2.sort()
    i = 0
    p = len(l)
    for f2 in l2 :
        f = f2[1] + "." + str(f2[0])  + "." + f2[2] 
        i+=1
        print(i,"/",p,  f)

        with open("c_results_LB_"+path_to_results, "a") as results :
            try :
                g = Graph(networkx.read_edgelist(path_to_lib+f))
                n = g.order()
                
                g.relabel(inplace=True,perm=list(range(n)))

                t0 = time.time()
                K = lower_bound2(g)
                temps = time.time() - t0
                results.write(str(n) + " " + str(temps) + " " + str(K) +'\n')
                print(str(n) + " " + str(temps) + " " + str(K) +'\n')
            except :
                pass


def plot_precomp(gen = "gnp") :
    results = open("c_results_IC_"+gen, 'r')
    tab = []
    x = []
    y= []
    line = results.readline()
    while line != '' :
        
        n, temps, K = list(map(float,line[:len(line)-1].split(" ")))
        
        
        x.append(n)
        y.append(temps)
            
        line = results.readline()
        
    results.close()
    plt.plot(x,y,'r.')
    plt.xlabel("Number of vertices")
    plt.ylabel("Execution time (s)")
    plt.show()





def test_lib(path_to_lib, path_to_results) :
    l = os.listdir(path_to_lib)
    l2 = []
    for f in l : 
        f2 = f.split(".")
        l2.append((int(f2[1]),f2[0],f2[2]))
    l2.sort()
    i = 0
    p = len(l)
    for f2 in l2 :
        f = f2[1] + "." + str(f2[0])  + "." + f2[2] 
        i+=1
        print(i,"/",p,  f)
        #if f2[0] <= 39 :
        #    continue
        #if f2[0] ==35 and f2[1] < "grafo3117.38.edgelist" :
        #    continue
        with open("c_results_"+path_to_results, "a") as results :
            g = Graph(networkx.read_edgelist(path_to_lib+f))
            n = f2[0]
            perm = shuffle(list(range(n)))
            g.relabel(inplace = True, perm = perm)
            n, m, temps, pl, LB, UB, diam = timed_run(g)
            results.write(str(n) + " " + str(m) + " " + str(temps) + " " + str(pl) + " " + str(LB) + " " + str(UB) + " " + str(diam)+'\n')
            with open("c_results_naive_"+path_to_results, "a") as results2 :

                _, _, temps2, _, _, _, _ = timed_run_naive(g)
                results2.write(str(n) + " " + str(m) + " " + str(temps2) + " " + str(pl) + " " + str(LB) + " " + str(UB) + " " + str(diam)+'\n')

#test_lib("./graphs/edgelist_rome/", "rome")






def plot_results(gen = "gnp", res = "all") :
    results = open("c_results_"+gen, 'r')
    tab = []
    x = []
    y= []
    line = results.readline()
    i=0
    while line != '' :
        
        n, m, temps, pl, LB, UB, diam = list(map(float,line[:len(line)-1].split(" ")))
        if temps > 600 :
            i+=1
            line = results.readline()

            continue
        if res == "poly" :
            if LB == UB : 
                    
                tab.append( (n,m,temps,pl,LB,UB,diam))
                x.append(n)
                y.append(temps)
        elif res == "exp" :
            if LB != UB : 
                tab.append( (n,m,temps,pl,LB,UB,diam))
                x.append(n)
                y.append(temps)
        else :
            tab.append( (n,m,temps,pl,LB,UB,diam))
            x.append(n)
            y.append(temps)
            
        line = results.readline()
        
    results.close()
    print(len(x), " graphs in total")
    print(i," graphs could not be computed")
    plt.plot(x,y,'r.')
    plt.xlabel("Number of vertices")
    plt.ylabel("Execution time (s)")
    plt.show()

def plot_results_compare(gen = 'op') :
    results = open("c_results_"+gen, 'r')
    results2 = open("c_results_naive_"+gen, 'r')
    tab = []
    x = []
    y= []
    y2= []
    line = results.readline()
    line2 = results2.readline()
    i=0
    while line != '' and line2 != '':
        
        n, m, temps, pl, LB, UB, diam = list(map(float,line[:len(line)-1].split(" ")))
        _, _, temps2, _, _, _, _ = list(map(float,line2[:len(line2)-1].split(" ")))
        if temps > 600 :
            i+=1
            line = results.readline()
            line2 = results2.readline()

            continue
        
        
        x.append(n)
        y.append(temps)
        y2.append(temps2)    
        line = results.readline()
        line2 = results2.readline()




    results.close()
    print(len(x), " graphs in total")
    print(i," graphs could not be computed")
    plt.plot(x,y,'rx')
    plt.plot(x,y2,'b+')
    plt.xlabel("Number of vertices")
    plt.ylabel("Execution time (s)")
    plt.show()

plot_results_compare('op')
plot_results_compare('rome')     


#print(g.edges(labels=False))

"""edges = [(0, 1), (0, 6), (0, 35), (1, 2), (2, 3), (2, 4), (2, 6), (3, 4), (4, 5), (4, 6), (5, 6), (6, 7), (6, 8), (6, 23), (6, 24), (6, 28), (7, 8), (8, 9), (8, 14), (8, 16), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (17, 21), (18, 19), (18, 20), (19, 20), (20, 21), (21, 22), (21, 23), (22, 23), (23, 24), (24, 25), (24, 26), (25, 26), (26, 27), (27, 28), (28, 29), (28, 31), (29, 30), (29, 31), (30, 31), (31, 32), (31, 33), (31, 34), (31, 35), (32, 33), (33, 34), (34, 35)]


g = Graph(edges)

run(g, [20, 10, 4, 30, 29, 18, 34, 26, 23, 9, 7, 15, 32, 24, 12, 2, 16, 33, 19, 25, 21, 31, 35, 8, 0, 17, 13, 28, 11, 14, 3, 27, 22, 1, 6, 5])"""