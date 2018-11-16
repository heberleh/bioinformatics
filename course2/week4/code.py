import sys


aav_uniq = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
massv_uniq = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

getChar = {}
for i in range(len(massv_uniq)):
    getChar[massv_uniq[i]] = aav_uniq[i]



aav = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
massv = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

getInt = {}
for i in range(len(aav)):
    getInt[aav[i]] = massv[i]


def stringPepToIntList(string_peptide):
    return [getInt[aa] for aa in string_peptide]



def cyclicSpectrumInt(pep):
    prefix_mass = [0]    
    for aa in pep:
        prefix_mass.append(prefix_mass[-1]+aa)
    peptideMass = prefix_mass[len(pep)]
    cyclic_spectrum = [0]
    for i in range(len(pep)):
        for j in range(i+1, len(pep)+1):
            cyclic_spectrum.append(prefix_mass[j]-prefix_mass[i])
            if i > 0 and j < len(pep):
                cyclic_spectrum.append(peptideMass-(prefix_mass[j]-prefix_mass[i]))            
    return sorted(cyclic_spectrum)        

def linearSpectrumInt(pep):
    prefix_mass = [0]    
    for aa in pep:
        prefix_mass.append(prefix_mass[-1]+aa)
    linear_spectrum = [0]
    for i in range(len(pep)):
        for j in range(i+1, len(pep)+1):
            linear_spectrum.append(prefix_mass[j]-prefix_mass[i])    
    return sorted(linear_spectrum)

def consistent(speccheck, specref):
    for c in speccheck:
        if specref.count(c) < speccheck.count(c):
            return False
    return True

def score(peptide, spectrum, fn=cyclicSpectrumInt):
    count = 0
    pep_spectrum = fn(peptide)
    for mass in set(pep_spectrum):
        count += min([spectrum.count(mass),pep_spectrum.count(mass)])
    return count


def expandPeptides(peptides, maximum, masses):
    new_set = []
    for pep in peptides:
        for aa in masses:
            if sum(pep)+aa <= maximum:
                new_pep = pep[:]
                new_pep.append(aa) 
                new_set.append(new_pep)
    return new_set

def trim(peptides, spectrum, n):
    if len(peptides) > n:
        tuples = []
        for pep in peptides:
            tuples.append((pep, score(pep, spectrum, fn=linearSpectrumInt)))
        tuples = sorted(tuples, key=lambda tup: tup[1], reverse=True)
        cutoff = tuples[n-1][1]
        for i in range(n+1,len(tuples)):
            if tuples[i][1] < cutoff:
                return [tup[0] for tup in tuples[:i]]
    return peptides


def leaderboardCyclePeptideSequencing(spectrum, n, masses):
    spectrum = sorted(spectrum)
    leaderboard = [[]]
    leader_peptides = []    
    leader_score = 0
    while len(leaderboard) > 0:
        good_peptides = expandPeptides(leaderboard, spectrum[-1], masses)          
        for pep in leaderboard:
            if sum(pep) == spectrum[-1]:                       
                pep_score = score(pep, spectrum)
                if pep_score > leader_score:
                    leader_peptides = [pep]
                    leader_score = pep_score
                elif pep_score == leader_score:
                    leader_peptides.append(pep)                    
        leaderboard = trim(good_peptides, spectrum, n)
    return leader_peptides


lines = sys.stdin.read().splitlines()
#pep_string = lines[0]
#spectrum = [int(m) for m in lines[1].split(' ')]

# peptide = stringPepToIntList(pep_string)
# cy_score = score(peptide, spectrum, fn=cyclicSpectrumInt)
# print(cy_score)

n = int(lines[0])
spectrum = [int(m) for m in lines[1].split(' ')]
final_peps = leaderboardCyclePeptideSequencing(spectrum, n, masses=range(57,201))
print(' '.join(sorted(['-'.join(map(str,pep)) for pep in final_peps], reverse=True)))