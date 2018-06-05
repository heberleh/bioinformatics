
import sys # you must import "sys" to read from STDIN

d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def reverse_complement(dna):
    return "".join([d[b] for b in dna])[::-1]

dna = sys.stdin.read().splitlines()[0]
print(reverse_complement(dna))
