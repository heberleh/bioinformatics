
import sys # you must import "sys" to read from STDIN

def match(text,word):    
    s = len(word)
    l = []
    for i in range(len(text)-s+1):
        if text[i:i+s] == word:
            l.append(i)
    return l

lines = sys.stdin.read().splitlines()

print (' '.join([str(s) for s in match(lines[1],lines[0])]))
# print (len([str(s) for s in match(lines[1],lines[0])]))
# print ("dna len", len(lines[1]))
#CTTGATCAT