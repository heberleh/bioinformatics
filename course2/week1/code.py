
import sys

def composition_k(text, k):
    kmers = set()
    for i in range(len(text) - k + 1):
        kmers.add(text[i:i+k])
    return kmers


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


lines = sys.stdin.read().splitlines()

# #ex1
# for kmer in composition_k(lines[1], int(lines[0])):
#     print(kmer)

# #ex2
# print(combine(lines))

#ex3
graph = adjacency_list(lines)
for source in graph:
    line = source + " -> "   
    for target in graph[source]:
        line += target +','
    if (len(graph[source])>0):
        print(line[:-1])