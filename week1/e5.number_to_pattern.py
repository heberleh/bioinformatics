
import sys # you must import "sys" to read from STDIN

d = ['A','C','G','T']
def numberToPattern(number,k):
    if k == 1:
        return d[number]
    else:        
        return numberToPattern(number // 4, k-1) + d[number % 4]
lines = sys.stdin.read().splitlines()
print(numberToPattern(int(lines[0]),int(lines[1])))