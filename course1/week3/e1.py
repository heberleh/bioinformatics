
import sys # you must import "sys" to read from STDIN
#import math 

def skew(genome):    
    c = 0
    v = []
    for i in range(len(genome)):
        if genome[i] == 'C':
            c-=1
        elif genome[i] == 'G':
            c+=1
        v.append(c)
    min_c = min(v)
    result = []
    for i in range(len(genome)):
        if v[i] == min_c:
            result.append(i+1)
    return result

def hamming(str1, str2):
    h = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            h+=1
    return h

def sufix(word):
    return word[1:]

nucleotides = ['A','C','G','T']

def neighbors(word, d):
    if d == 0:
        return set([word])
    if len(word) == 1:
        return set([nucleotides])
    neighborhood = set()
    sufix_word = sufix(word)
    for suf in neighbors(sufix_word, d):
        if hamming(sufix_word, suf) < d:
            for nt in nucleotides:
                neighborhood.add(nt+suf)
        else:
            neighborhood.add(word[0]+suf)
    return neighborhood

def findAproximateMatches(genome, word, k):
    s = len(word)
    l = []
    k += 1
    for i in range(len(genome)-s+1):
        if hamming(genome[i:i+s], word) < k:
            l.append(i)
    return l

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

def computingFrequency(text, k, d):
    freq = [0] * (4**k) 
    for i in range(len(text)-k+1):
        for neighbor in neighbors(text[i:i+k],d): #text[i:i+k]                        
            freq[patternToNumber(neighbor)] += 1
    return freq

drev = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def reverse_complement(dna):
    return "".join([drev[b] for b in dna])[::-1]

def computingFrequencyRegRev(text, k, d):
    freq = [0] * (4**k)     
    for i in range(len(text)-k+1):
        # neighborhood = neighbors(text[i:i+k],d)
        # for neighbor in neighborhood: #text[i:i+k]                        
        #     freq[patternToNumber(neighbor)] += 1
        # neighborhood = neighbors(reverse_complement(text[i:i+k]),d)
        # for neighbor in neighborhood: #text[i:i+k]                        
        #     freq[patternToNumber(neighbor)] += 1
        
        for neighbor in neighbors(text[i:i+k],d): #text[i:i+k]                        
            freq[patternToNumber(neighbor)] += 1      
            freq[patternToNumber(reverse_complement(neighbor))] += 1
    return freq

def freqWords(text, k, d):
    words = set()
    freq = computingFrequency(text, k, d)
    max_freq = max(freq)
    for i in range(4**k):
        if freq[i] == max_freq:
            words.add(numberToPattern(i, k))
    return words

def freqWordsRegRev(text, k, d):
    words = set()
    freq1 = computingFrequency(text, k, d)
    freq2 = computingFrequency(reverse_complement(text), k, d)
    freq = [freq1[i]+freq2[i] for i in range(len(freq1))]
    max_freq = max(freq)
    for i in range(4**k):
        if freq[i] == max_freq:
            words.add(numberToPattern(i, k))
    return words

def MotifEnumeration(dna, k, d):
    patterns = set()
    for i in range(len(dna[0])-k+1):
        for neighbor in neighbors(dna[0][i:i+k],d): #text[i:i+k] 
            all_sub_dna_has = True            
            for sub_dna in dna:
                appear_in_sub_dna = False      
                for neighbor2 in neighbors(neighbor,d):
                    if neighbor2 in sub_dna:
                        appear_in_sub_dna = True
                        break
                if appear_in_sub_dna == False:
                    all_sub_dna_has = False
                    break
            if all_sub_dna_has:
                patterns.add(neighbor)
    return patterns

def distance(pattern, dna_list):
    k = len(pattern)
    dist = 0
    for dna in dna_list:        
        hamdist = float("inf")
        for i in range(len(dna)-k+1):
            pattern2 = dna[i:i+k]
            # for neighbor in neighbors(pattern2, k):
            hamdist_aux = hamming(pattern2, pattern)
            if hamdist_aux < hamdist:
                hamdist = hamdist_aux
        dist = dist + hamdist
    return dist

def medianString(dna_list, k):
    dist = float("inf")
    for i in range(0,(4**k)):
        pattern = numberToPattern(i,k)
        dist_aux = distance(pattern, dna_list)
        if dist_aux < dist:
            dist = dist_aux
            median = pattern
    return median



def ProfileMostProbableKMerProblem(dna, k, probs_matrix):
    idx = {'A':0, 'C': 1, 'G': 2, 'T':3}
    max_prob = 0
    selected_word = ""
    for i in range(len(dna)-k+1):
        word = dna[i:i+k]
        prob = 0
        for w in range(k):
            prob += float(probs_matrix[idx[word[w]]][w])
        if prob > max_prob:
            max_prob = prob
            selected_word = word
    return selected_word




lines = sys.stdin.read().splitlines()

#line0 = lines[0].split(" ")
dna = lines[0]
#k = int(line0[0])
#d = int(line0[1])
k = int(lines[1])
probs_matrix = []
for i in range(2,len(lines)):
    probs_matrix.append(lines[i].split(' '))

#print(medianString(dna,k))
#print(' '.join([str(s) for s in MotifEnumeration(dna,k,d)]))
print(ProfileMostProbableKMerProblem(dna,k, probs_matrix))