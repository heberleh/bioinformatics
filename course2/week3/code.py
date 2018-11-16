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



class bidict(dict):
    def __init__(self, *args, **kwargs):
        super(bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.iteritems():
            self.inverse.setdefault(value,[]).append(key) 

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key) 
        super(bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value,[]).append(key)        

    def __delitem__(self, key):
        self.inverse.setdefault(self[key],[]).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]: 
            del self.inverse[self[key]]
        super(bidict, self).__delitem__(key)


mass = {'G': 57,'A': 71,'S': 87,'P': 97,'V': 99,'T': 101,'C': 103,'I': 113,'L': 113,'N': 114,'D': 115,'K': 128,'Q': 128,'E': 129,'M': 131,'H': 137,'F': 147,'R': 156,'Y': 163,'W': 186}


aav = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
massv = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

aav_uniq = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
massv_uniq = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

pepdic = {}
for i in range(len(massv_uniq)):
    pepdic[massv_uniq[i]] = aav_uniq[i]

# # removed 'I': 113,'Q': 128
# mass_uniq = bidict({'G': 57,'A': 71,'S': 87,'P': 97,'V': 99,'T': 101,'C': 103,'L': 113,'N': 114,'D': 115,'K': 128, 'E': 129,'M': 131,'H': 137,'F': 147,'R': 156,'Y': 163,'W': 186})



class Memoize:
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args):
        if args not in self.memo:
            self.memo[args] = self.fn(*args)
        return self.memo[args]


def linearSpectrum(pep):
    prefix_mass = [0]    
    for aa in pep:
        prefix_mass.append(prefix_mass[-1]+mass[aa])

    linear_spectrum = [0]
    for i in range(len(pep)):
        for j in range(i+1, len(pep)+1):
            linear_spectrum.append(prefix_mass[j]-prefix_mass[i])    

    return sorted(linear_spectrum)

def linearSpectrumInt(pep):
    prefix_mass = [0]    
    for aa in pep:
        prefix_mass.append(prefix_mass[-1]+aa)

    linear_spectrum = [0]
    for i in range(len(pep)):
        for j in range(i+1, len(pep)+1):
            linear_spectrum.append(prefix_mass[j]-prefix_mass[i])    

    return sorted(linear_spectrum)    

def cyclicSpectrum(pep):
    prefix_mass = [0]    
    for aa in pep:
        prefix_mass.append(prefix_mass[-1]+mass[aa])
    
    peptideMass = prefix_mass[len(pep)]
    cyclic_spectrum = [0]
    for i in range(len(pep)):
        for j in range(i+1, len(pep)+1):
            cyclic_spectrum.append(prefix_mass[j]-prefix_mass[i])
            if i > 0 and j < len(pep):
                cyclic_spectrum.append(peptideMass-(prefix_mass[j]-prefix_mass[i]))            
    return sorted(cyclic_spectrum)


def countingPeptides(target):
    @Memoize
    def countingPeptidesAux(m):
        count = 0
        for aa_mass in massv_uniq:
            totalmass = aa_mass + m
            if totalmass == target:
                count+=1
            elif totalmass < target:                
                count+=countingPeptidesAux(totalmass)               
        return count
    count = 0    
    for m in massv_uniq:        
        if m < target:
            count+=countingPeptidesAux(m)
        elif m == target:
            count+=1
    return count


#    CyclopeptideSequencing(Spectrum)
#         Peptides ← a set containing only the empty peptide
#         while Peptides is nonempty
#             Peptides ← Expand(Peptides)
#             for each peptide Peptide in Peptides
#                 if Mass(Peptide) = ParentMass(Spectrum)
#                     if Cyclospectrum(Peptide) = Spectrum
#                         output Peptide
#                     remove Peptide from Peptides
#                 else if Peptide is not consistent with Spectrum
#                     remove Peptide from Peptides

def expandPeptides(peptides, maximum, masses):
    new_set = []
    for pep in peptides:
        for aa in masses:
            if sum(pep)+aa <= maximum:
                new_pep = pep[:]
                new_pep.append(aa) 
                new_set.append(new_pep)
    return new_set


def consistent(speccheck, specref):
    for c in speccheck:
        if specref.count(c) < speccheck.count(c):
            return False
    return True

def pepNumberToString(pep_masses):
    seq = ""
    for m in pep_masses:
        seq+= pepdic[m]
    return seq

def cycloPeptideSequencing(spectrum):
    spectrum = sorted(spectrum)
    final_peptides = []
    peptides = [[]]
    masses = massv_uniq #[m for m in massv_uniq if m in spectrum]
    while len(peptides) > 0:
        peptides = expandPeptides(peptides, spectrum[-1], masses)
        good_peptides = []      
        for pep in peptides:
            if sum(pep) == spectrum[-1]:                        
                if cyclicSpectrum(pepNumberToString(pep)) == spectrum:
                    final_peptides.append(pep)                                    
            elif consistent(linearSpectrumInt(pep) , spectrum):                        
                good_peptides.append(pep)
        peptides = good_peptides     
    return final_peptides

lines = sys.stdin.read().splitlines()
#print(translate(lines[0]))


# dna = lines[0]
# prot = lines[1]
#print(proteinToRNAs(prot))
# print('\n'.join(findEncodingSubstrings(dna,prot)))

#pep = lines[0]
#print(' '.join(map(str,linearSpectrum(pep))))
#print(' '.join(map(str,cyclicSpectrum(pep))))

#target = int(lines[0])
#print(countingPeptides(target))

# target = 0
# import math
# target = 800
# x1 = 1024#target
# y1 = 14712706211#countingPeptides(target)
# lastc = y1+100000
# currenty = y1
# i = 0
# while i < 10:
#     i+=1
#     target = target + 600
#     x2 = 1307#target
#     y2 = 34544458837656#countingPeptides(x2)
#     currenty = y2
#     logc = math.log10(y1/y2)/(x1-x2)
#     lastc = 10**logc
    
#     print("solution: ")
#     print(lastc)
   
# spectrum = [int(m) for m in lines[0].split(' ')]
# print(' '.join(sorted(['-'.join(map(str,pep)) for pep in cycloPeptideSequencing(spectrum)], reverse=True)))

#print(proteinToRNAs("PRTEIN"))
# print(len(proteinToRNAs("MASS")))
#print('\n'.join(findEncodingSubstrings(dna,prot)))

# peps = cycloPeptideSequencing([0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416])
# for pep in peps:
#     print(pepNumberToString(pep))

# spectrum = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
# for pep in ['TMIA', 'IAMT', 'MLAT', 'MTAL', 'MAIT', 'TAIM']:
#     print(pep)
#     c = cyclicSpectrum(pep)
#     print(c)
#     print(consistent(c,sorted(spectrum)))
#     print("------")


spectrum = [0,71,99,101,103,128,129,199,200,204,227,230,231,298,303,328,330,332,333]
for pep in ['QCV', 'CTV', 'ETC', 'TCE', 'CTQ', 'AQV']:
    print(pep)
    c = linearSpectrum(pep)
    print(c)
    print(consistent(c,sorted(spectrum)))
    print("------")