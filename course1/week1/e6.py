
import sys # you must import "sys" to read from STDIN

d = ['A','C','G','T']

def id(n):
    if n == 'A':
        return 0
    elif n == 'C':
        return 1
    elif n == 'G':
        return 2
    else:
        return 3

#id = {'A':0,'C':1,'G':2,'T':3}

def patternToNumber(pattern):
    if len(pattern) ==  0:
        return 0
    else:
        return 4 * patternToNumber(pattern[:-1]) + id(pattern[-1])

def numberToPattern(number,k):
    if k == 1:
        return d[number]
    else:        
        return numberToPattern(number // 4, k-1) + d[number % 4]

def computingFrequency(text, k):
    freq = [0] * (4**k)
    for i in range(len(text)-k+1):
        freq[patternToNumber(text[i:i+k])] += 1
    return freq

def freqWords(text, k):
    words = set()
    freq = computingFrequency(text, k)
    max_freq = max(freq)
    for i in range(4**k):
        if freq[i] == max_freq:
            words.add(numberToPattern(i, k))
    return words

def clumpFinding(genome, k, L, t):
    freq_words = set()
    t = t-1
    for i in xrange(len(genome)-L+1):
        window = genome[i:i+L]
        counts = {}
        for j in xrange(len(window)-k+1):
            word = window[j:j+k]
            if word not in counts:
                counts[word] = 0
            counts[word] += 1
        for kmer in counts:
            if counts[kmer] > t:
                freq_words.add(kmer)
    return freq_words


lines = sys.stdin.read().splitlines()
par = lines[1].split(' ')

#print (' '.join([str(s) for s in clumpFinding(lines[0],int(par[0]),int(par[1]),int(par[2]))]))
print(len(clumpFinding(lines[0],int(par[0]),int(par[1]),int(par[2]))))
