
import sys # you must import "sys" to read from STDIN
id = {'A':0,'C':1,'G':2,'T':3}

def numberToPattern(number,k):
    if k == 1:
        return d[number]
    else:        
        return numberToPattern(number // 4, k-1) + d[number % 4]


print(patternToNumber(sys.stdin.read().splitlines()[0]))