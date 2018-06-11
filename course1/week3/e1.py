
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
            count[symbol].append(1)

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


def distance(pattern, dna_list):
    k = len(pattern)
    dist = 0
    for dna in dna_list:        
        hamdist = float("inf")
        for i in range(len(dna)-k+1):
            pattern2 = dna[i:i+k]            
            hamdist_aux = hamming(pattern2, pattern)
            if hamdist_aux < hamdist:
                hamdist = hamdist_aux
        dist = dist + hamdist
    return dist

def hamming(str1, str2):
    h = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            h+=1
    return h


def medianString(dna_list, k):
    dist = float("inf")
    medians = []
    for i in range(0,(4**k)):
        pattern = numberToPattern(i,k)
        dist_aux = distance(pattern, dna_list)
        if dist_aux <= dist:
            dist = dist_aux
            median = pattern
            medians.append(pattern)
    #print(medians)
    return median

nucleotides = ['A','C','G','T']

def id(n):
    if n == 'A':
        return 0
    elif n == 'C':
        return 1
    elif n == 'G':
        return 2
    else:
        return 3

def patternToNumber(pattern):
    if len(pattern) ==  0:
        return 0
    else:
        return 4 * patternToNumber(pattern[:-1]) + id(pattern[-1])

def numberToPattern(number,k):
    if k == 1:
        return nucleotides[number]
    else:        
        return numberToPattern(number // 4, k-1) + nucleotides[number % 4]





#dna = """CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC
# GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC
# GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"""
# dna = dna.split('\n')
# print(medianString(dna,7))
profile = {'A': [0.4,0.3,0.0,0.1,0.0,0.9],
'C': [0.2,0.3,0.0,0.4,0.0,0.1],
'G': [0.1,0.3,1.0,0.1,0.5,0.0],
'T': [0.3,0.1,0.0,0.4,0.5,0.0]}
print(Pr('AAGTTC', profile))

# lines = sys.stdin.read().splitlines()

# #line0 = lines[0].split(" ")
# pattern = lines[0]
# #k = int(line0[0])
# #d = int(line0[1])
# #k = int(lines[1])
# #t = int(line0[1])

# dna_list = lines[1].split(' ')
# # for i in range(1,len(lines)):
# #     dna_list.append(lines[i])

# #print(medianString(dna,k))
# #print(' '.join([str(s) for s in MotifEnumeration(dna,k,d)]))
# #print('\n'.join(distance(pattern, dna_list)))
# print(distance(pattern, dna_list))