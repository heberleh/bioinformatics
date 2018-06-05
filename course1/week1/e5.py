
import sys # you must import "sys" to read from STDIN


kmer = []
d = ['A','C','G','T']
id = {'A':0,'C':1,'G':2,'T':3}

freq = []
k = 2
for i in range(4**k):
    freq.append(0)

def patternToNumber(pattern):
    if len(pattern) ==  0:
        return 0
    else:
        return 4 * patternToNumber(pattern[:-1]) + id[pattern[-1]]


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




def match(text,word):    
    s = len(word)
    l = []
    for i in range(len(text)-s+1):
        if text[i:i+s] == word:
            l.append(i)
    return l

lines = sys.stdin.read().splitlines()

print (' '.join([str(s) for s in match(lines[0], "CTTGATCAT")]))
# print (len([str(s) for s in match(lines[1],lines[0])]))
# print ("dna len", len(lines[1]))