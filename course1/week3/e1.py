
import sys

def GreedyMotifSearch(Dna, k, t):
    bestMotifs = []
    for i in range(0, t):
        bestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        motifs = []
        motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    return (bestMotifs)

def Score(motifs):
    count = 0
    k = len(motifs[0])
    t = len(motifs)
    consensusMotif = Consensus(motifs)
    for i in range(t):
        for j in range(k):
            if motifs[i][j] != consensusMotif[j]:
                count += 1
    return count

def Count(motifs):
    count = {}
    k = len(motifs[0])

    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    return count

def Profile(motifs):
    profile = {}
    t = len(motifs)
    k = len(motifs[0])
    countMotifs = Count(motifs)

    for symbol in "ACGT":
        profile[symbol] = []

    for x in countMotifs:
        for y in countMotifs[x]:
            z = y/float(t)
            profile[x].append(z)

    return profile

def Consensus(motifs):
    k = len(motifs[0])
    count = Count(motifs)

    consensus = ""
    for j in range(k):
        M = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > M:
                M = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def ProfileMostProbableKmer(text, k, profile):
    mostProbVal = -1
    mostProbKmer = ''
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probKmerVal = Pr(kmer, profile)
        if probKmerVal > mostProbVal:
            mostProbVal = probKmerVal
            mostProbKmer = kmer
    return mostProbKmer

def Pr(text, profile):
    P = 1
    for i in range(len(text)):
        P = P * profile[text[i]][i]
    return P


lines = sys.stdin.read().splitlines()

line0 = lines[0].split(" ")
#dna = lines[0]
k = int(line0[0])
#d = int(line0[1])
#k = int(lines[1])
t = int(line0[1])

dna_list = []
for i in range(1,len(lines)):
    dna_list.append(lines[i])

#print(medianString(dna,k))
#print(' '.join([str(s) for s in MotifEnumeration(dna,k,d)]))
print('\n'.join(GreedyMotifSearch(dna_list, k, t)))