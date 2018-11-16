import sys, random

def find_degrees(graph):
    in_degree = {}
    out_degree = {}
    for source in graph:
        for target in graph[source]:
            in_degree[target] = 0
            out_degree[target] = 0

    for source in graph:
        in_degree[source] = 0
        out_degree[source] = 0

    for source in graph:
        out_degree[source] = len(graph[source])
    
    for source in graph:
        for target in graph[source]:
            in_degree[target]+=1

    return {"in":in_degree, "out":out_degree}

def prefixPair(pair):
    return "|".join([p[:-1] for p in pair.split("|")])

def sufixPair(pair):
    return "|".join([p[1:] for p in pair.split("|")])

def deBrujinPairs(pairs):
    adj_list = {prefixPair(pair): [] for pair in pairs}
    for pair in pairs:
        adj_list[prefixPair(pair)].append(sufixPair(pair))
    return adj_list

def deBrujinKmers(kmers):
    adj_list = {kmer[:-1]: [] for kmer in kmers}
    for kmer in kmers:
        adj_list[kmer[:-1]].append(kmer[1:])
    return adj_list    

def eulerianPath(graph):
    degree = find_degrees(graph)
    start = None
    for key in degree['in']:
        if degree['out'][key] - degree['in'][key] > 0:
            start = key
            break
    for key in degree['in']:
        if degree['in'][key] - degree['out'][key] > 0:
            end = key
            break   

    if end in graph:
        graph[end].append(start)
    else:
        graph[end] = [start]

    stack = [start]
    path = []
    while stack:		
        if graph[stack[0]]:
            w = random.choice(graph[stack[0]])
            graph[stack[0]].remove(w)
            stack.insert(0,w)
        else:
            a=stack[0]
            path.append(stack[0])
            stack.remove(a)	
    path.reverse()
    return path[:-1]

def concat_pairs(s1, s2, k):
    #GTG|GTG->TGG|TGA
    splited1 = s1.split("|")
    splited2 = s2.split("|")
    return splited1[0][:k-1]+splited2[0][-1]+'|'+splited1[1][:k-1]+splited2[1][-1]

def path_convert(path, k):
    new_path = []
    for i in range(len(path)-1):
        new_path.append(concat_pairs(path[i], path[i+1], k))
    return new_path

def string_reconstruction(pairs, k, d):
    s = ""
    for i in range(0,len(pairs)-1):
        s += pairs[i][0]    

    for i in range(k):
        s+=pairs[len(pairs)-1][i]    
   
    start = (len(pairs)-1)-d
    end = len(pairs)-1

    for i in range(start, end):
        s+=pairs[i][len(pairs[i])-k]  # [_ _ _ | _ _ _] 0 1 2 |3| 4 5 6 -> 7-3 = 4   

    for i in range(k+1,(2*k)+1):
        s+=pairs[len(pairs)-1][i]    

    return s

"""
 MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
"""

def rm_e(g, src, trg):
    g[src].remove(trg)
    if len(g[src]) == 0:
        g.pop(src)

def degreesOne(degree, node):
    return degree['in'][node] == 1 and degree['out'][node] == 1

def maximalNonBranchingPaths(graph):
    paths = []
    degree = find_degrees(graph)
    used_nodes = set()
    for source in graph:
        if degree['in'][source] != 1 or degree['out'][source] != 1:
            if degree['out'][source] > 0:
                for target in graph[source]:
                    used_nodes.add(source)
                    nonBranchingPath = [source, target]    
                    current = target                                    
                    while degreesOne(degree, current):
                        used_nodes.add(current)
                        nonBranchingPath.append(graph[current][0])
                        current = graph[current][0]
                    paths.append(nonBranchingPath)
    
    for source in graph:
        if source not in used_nodes and degreesOne(degree, source):
            used_nodes.add(source)
            current = graph[source][0]
            nonBranchingPath = [source]
            while degreesOne(degree, current):
                used_nodes.add(current)
                nonBranchingPath.append(current)
                current = graph[current][0]
                if current == source:
                    nonBranchingPath.append(current)
                    paths.append(nonBranchingPath)
                    break

    return paths

def concat_kmers(s1, s2):
    #GTG + TGG = GTGG
    return s1+s2[-1]

def contigsKmers(kmers):
    graph = deBrujinKmers(kmers)
    paths = maximalNonBranchingPaths(graph)
    contigs = [] 
    for path in paths:
        word = path[0]
        for i in range(1, len(path)):
            word = concat_kmers(word, path[i])
        contigs.append(word)
    return contigs

def join_genome(path):
    result = path[0]
    for i in range(1,len(path)):
        result += path[i][-1]
    return result

lines = sys.stdin.read().splitlines()
# k = int(lines[0].split(' ')[0])
# d =  int(lines[0].split(' ')[1])
# pairs = []
# for i in range(1,len(lines)):
#     pairs.append(lines[i])

# graph = deBrujinPairs(pairs)
# path = eulerianPath(graph)
# path = path_convert(path,k)
# s = string_reconstruction(path, k, d)
# print(s)

#print(string_reconstruction(lines[1:],k,d))


# reading graph
# graph = {}
# for line in lines:
#     l = line.split(" -> ")
#     graph[l[0]] = []
#     for trg in l[1].split(','):
#         graph[l[0]].append(trg)
# for path in maximalNonBranchingPaths(graph):
#     print(" -> ".join(path))


# print(' '.join(sorted(contigsKmers(lines))))

# quiz ex 1
# kmers = lines
# graph = deBrujinKmers(kmers)
# result = eulerianPath(graph)
# print(join_genome(result))

# quiz ex 3
k=3
d=1
graph = deBrujinPairs(lines)
path = eulerianPath(graph)
path = path_convert(path,k)
s = string_reconstruction(path, k, d)
print(s)

True or False: every Eulerian path in the paired de Bruijn graph constructed from a (k, d)-mer composition must spell out a solution to the String Reconstruction from Read-Pairs Problem.
False

True or False: read breaking can transform reads with imperfect coverage into reads with perfect coverage.
True