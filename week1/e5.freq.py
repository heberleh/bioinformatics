
import sys # you must import "sys" to read from STDIN
kmer = []
d = ['A','C','G','T']
id = {'A':0,'C':1,'G':2,'T':3}

def patternToNumber(pattern):
    if len(pattern) ==  0:
        return 0
    else:
        return 4 * patternToNumber(pattern[:-1]) + id[pattern[-1]]

def computingFrequency(text, k):
    freq = [0] * (4**k)
    for i in range(len(text)-k+1):
        freq[patternToNumber(text[i:i+k])] += 1
    return freq

lines = sys.stdin.read().splitlines()
print (' '.join([str(s) for s in computingFrequency(lines[0], int(lines[1]))]))