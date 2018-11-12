import sys, random

def find_unbalanced_nodes(graph):
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

def prefix(pair):
    return "|".join([p[:-1] for p in pair.split("|")])

def sufix(pair):
    return "|".join([p[1:] for p in pair.split("|")])

def deBrujin(pairs):
    adj_list = {prefix(pair): [] for pair in pairs}
    for pair in pairs:
        adj_list[prefix(pair)].append(sufix(pair))
    return adj_list

def eulerianPath(graph):
    degree = find_unbalanced_nodes(graph)
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

lines = sys.stdin.read().splitlines()
k = int(lines[0].split(' ')[0])
d =  int(lines[0].split(' ')[1])
pairs = []
for i in range(1,len(lines)):
    pairs.append(lines[i])

# graph = deBrujin(pairs)
# path = eulerianPath(graph)
# path = path_convert(path,k)
# s = string_reconstruction(path, k, d)
# print(s)

print(string_reconstruction(lines[1:],k,d))