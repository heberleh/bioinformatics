
import sys # you must import "sys" to read from STDIN
id = {'A':0,'C':1,'G':2,'T':3}

def patternToNumber(pattern):
    if len(pattern) ==  0:
        return 0
    else:
        return 4 * patternToNumber(pattern[:-1]) + id[pattern[-1]]

print(patternToNumber(sys.stdin.read().splitlines()[0]))