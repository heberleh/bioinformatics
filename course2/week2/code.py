import sys
from itertools import product

def composition_k(text, k):
    kmers = set()
    for i in range(len(text) - k + 1):
        kmers.add(text[i:i+k])
    return list(kmers)


def combine(kmers):
    combined = kmers[0]
    for i in range(1,len(kmers)):
        kmer = kmers[i]
        for j in range(len(kmer),0,-1):
            if combined.endswith(kmer[0:j]):
                combined = combined+kmer[j:len(kmer)]
                break
    return combined

def prefix_sufix(kmer1, kmer2):
    return kmer1[1:len(kmer1)] == kmer2[0:len(kmer2)-1]

def adjacency_list(kmers):
    adjacency_list = {}
    for kmer1 in kmers:
        if not kmer1 in adjacency_list:
            adjacency_list[kmer1] = set()           
    for i in range(len(kmers)):
        for j in range(len(kmers)):
            if i!=j and prefix_sufix(kmers[i], kmers[j]):
                adjacency_list[kmers[i]].add(kmers[j])
    return adjacency_list

def deBruijn_k(text, k):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])

    adjacency_list = {}
    for kmer1 in kmers:
        kmer = kmer1[0:len(kmer1)-1]
        if not kmer in adjacency_list:
            adjacency_list[kmer] = []    

    for i in range(len(kmers)):
        kmer = kmers[i]        
        adjacency_list[kmer[0:len(kmer)-1]].append(kmer[1:len(kmer)])
    return adjacency_list

def deBruijn_kmers(kmers):
    adjacency_list = {}
    
    for kmer in kmers:  
        prefix = kmer[0:len(kmer)-1]
        sufix = kmer[1:len(kmer)]
        
        if not prefix in adjacency_list:
            adjacency_list[prefix] = []

        adjacency_list[prefix].append(sufix)        

    return adjacency_list


def print_graph(graph):
    for source in graph:
        line = source + " -> "   
        for target in graph[source]:
            line += target +','
        if (len(graph[source])>0):
            print(line[:-1])


def rm_e(g, src, trg):
    g[src].remove(trg)
    if len(g[src]) == 0:
        g.pop(src)

def next(g, src):
    return g[src][0]

def get_start(g, unused, c):
    for s in c:
        if s in unused:
            return s
    for s in unused:
        for t in unused[s]:
            if t in c:
                return t
    
    print("error")

def traverse(c, start):    
    c2 = []
    c = c[1:]
    trg = None    

    size = len(c)
    index = 0   
    cs = c[index] 
    while cs != start:  
        index = (index + 1) % size      
        cs = c[index] 
        
    c2.append(cs)   
    for i in range(len(c)-1):
        index = (index + 1) % size        
        c2.append(c[index])    
    return c2        


# EulerianCycle(Graph)
#     form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
#     while there are unexplored edges in Graph
#         select a node newStart in Cycle with still unexplored edges
#         form Cycle' by traversing Cycle (starting at newStart) and then randomly walking 
#         Cycle = Cycle'
#     return Cycle
def eulerian_cycle(g):
    unused = {}
    for src in g:
        unused[src] = [trg for trg in g[src]]

    c = []
    start = list(unused)[0]
    now = unused[start][0] 
    rm_e(unused, start, now)
    c.append(start) 
    c.append(now)      
    while c[0] != c[-1]:        
        last = now      
        now = unused[now][0] 
        rm_e(unused, last, now)
        c.append(now)           

    while len(list(unused)) > 0:
        start = get_start(g, unused, c)  

        c2 = traverse(c, start)               
        now = start
        first = True
        while c2[0] != c2[-1]:
            if first:
                c2.append(now)
                first = False
            last = now            
            now = next(unused, now)
            rm_e(unused, last, now)            
            c2.append(now)           
        c = c2
    return c

def find_unbalanced_nodes(graph):
    ind = {}
    outd = {}
    for source in graph:
        for target in graph[source]:
            ind[target] = 0
            outd[target] = 0
    for source in graph:
        ind[source] = 0
        outd[source] = 0

    for source in graph:
        outd[source] += len(graph[source])
    
    for source in graph:
        for target in graph[source]:
            ind[target]+=1

    source = None
    for node in ind:
        if ind[node] > outd[node]:
            source = node
        
    target = None
    for node in ind:
        if ind[node] < outd[node]:
            target = node

    if source is None and list(graph)[0] != target:
        source = list(graph)[0]   
    elif source is None:
        source = list(graph)[1]

    return {"source":source, "target":target}

def eulerian_path(g):
    unbalanced = find_unbalanced_nodes(g)
    source = unbalanced['source']
    target = unbalanced['target']
    
    unused = {}
    for src in g:
        unused[src] = [trg for trg in g[src]]
    if source in g:
        g[source].append(target)
        unused[source].append(target)
    else:
        g[source] = [target]
        unused[source] = [target]
    
    c = []
    start = list(unused)[0]
    now = unused[start][0] 
    rm_e(unused, start, now)
    c.append(start) 
    c.append(now)      
    while c[0] != c[-1]:        
        last = now      
        now = unused[now][0] 
        rm_e(unused, last, now)
        c.append(now)           

    while len(list(unused)) > 0:
        start = get_start(g, unused, c)  

        c2 = traverse(c, start)               
        now = start
        first = True
        while c2[0] != c2[-1]:
            if first:
                c2.append(now)
                first = False
            last = now            
            now = next(unused, now)
            rm_e(unused, last, now)            
            c2.append(now)           
        c = c2

    c = traverse(c, target)
    
    return c

def join_genome(path):
    result = path[0]
    for i in range(1,len(path)):
        result += path[i][-1]
    return result

def binary_kmers(k):    
    return [''.join(map(str,l)) for l in list(product([0,1],repeat=k))]

        

lines = sys.stdin.read().splitlines()

graph = {}
# for line in lines:
#     l = line.split(" -> ")
#     graph[l[0]] = []
#     for trg in l[1].split(','):
#         graph[l[0]].append(trg)

# ex1
#result = eulerian_cycle(graph)

# ex2
# result = eulerian_path(graph)

# str_result = result[0]
# for i in range(1,len(result)):
#     str_result+="->"+result[i]
# print(str_result)

# #ex3
# kmers = [lines[i] for i in range(1,len(lines))]
# graph = deBruijn_kmers(kmers)
# result = eulerian_path(graph)
# print(join_genome(result))


#ex 4
# k = int(lines[0])
# kmers = binary_kmers(k)
# graph = getDeBrujin(kmers)
# path = EulerianCycle(graph)
# genome = join_genome(path[:-(k-1)])
# print(genome)