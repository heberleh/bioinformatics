
import sys # you must import "sys" to read from STDIN

def patternCount(text, code):
    c = 0
    k = len(code)
    for i in range(len(text)-k+1):        
        if text[i:i+k] == code:
            c += 1
    return c


lines = sys.stdin.read().splitlines() # read in the input from STDIN

code = lines[0]
kmer = lines[1]
solution = None
if len(lines) > 2:
    solution = lines[2]

k = len(kmer)

result = patternCount(code,kmer)
if solution == None or str(result) == solution:
    print(result)
else:
    print("Error.",result)