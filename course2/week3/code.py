import sys

amino_hash = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

codon_hash = {}
for codon in amino_hash:
    amino = amino_hash[codon]
    if amino in codon_hash:
        codon_hash[amino].append(codon)
    else:
        codon_hash[amino] = [codon]

def proteinToRNAs(amino_seq):
    sequences = []
    for amino in amino_seq:
        codons = codon_hash[amino]
        if len(sequences) == 0:
            sequences = codons
        else:
            new_sequences = []
            for seq in sequences:
                for codon in codons:            
                    new_sequences.append(seq+codon)
            sequences = new_sequences
    return sequences


def translate(rna):
    i = 0
    protein = ""
    while i < len(rna):
        amino = amino_hash[rna[i:i+3]]
        if amino != "*":
            protein += amino
        i += 3
    return protein

drev = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def reverse_complement(dna):
    return "".join([drev[b] for b in dna])[::-1]

def dna2rna(dna):
    return dna.replace('T','U')

def rna2dna(dna):
    return dna.replace('U','T')

def findEncodingSubstrings(genome, prot):
    rnas = proteinToRNAs(prot)
    dnas = []
    for rna in rnas:
        dna = rna2dna(rna)
        rdna = reverse_complement(dna)
        dnas.append(dna)
        dnas.append(rdna)
    return [dna for dna in dnas if dna in genome]


lines = sys.stdin.read().splitlines()
#print(translate(lines[0]))


dna = lines[0]
prot = lines[1]

#print(proteinToRNAs(prot))
print('\n'.join(findEncodingSubstrings(dna,prot)))