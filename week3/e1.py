
import sys # you must import "sys" to read from STDIN

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
        return word
    if len(word) == 1:
        return nucleotides
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


def motifEnumeration(dna, k, d):
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

lines = sys.stdin.read().splitlines()

line0 = lines[0].split(" ")
k = line0[0]
d = line0[1]
dna = []
for i in range(1,len(lines)):
    dna.append(lines[i])

print(' '.join([str(s) for s in motifEnumeration(dna,int(k),int(d))]))