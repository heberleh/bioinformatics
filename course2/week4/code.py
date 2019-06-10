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
        for pep in good_peptides:
            if sum(pep) == spectrum[-1]:                                       
                pep_score = score(pep, spectrum)
                if pep_score > leader_score:
                    leader_peptides = [pep]
                    leader_score = pep_score
                elif pep_score == leader_score:
                    leader_peptides.append(pep)                    
        leaderboard = trim(good_peptides, spectrum, n)
    print(leader_score)
    return leader_peptides

def convolution(spectrum, m):
    if 0 not in spectrum:
        spectrum.append(0)
    count = {}
    allmasses = sorted(spectrum, reverse=True)
    for i in range(len(allmasses)):
        for j in range(i+1,len(allmasses)):
            diff = abs(allmasses[i]-allmasses[j])    
            if diff > 56 and diff < 201:          
                if diff in count:
                    count[diff]+=1
                else:
                    count[diff]=1
    tuples = sorted([(key, count[key]) for key in count], key=lambda tup: tup[1], reverse=True)
    cutoff = tuples[m-1][1]
    for i in range(m+1,len(tuples)):
            if tuples[i][1] < cutoff:
                return sorted(list(set([tup[0] for tup in tuples[:i]])))
    return sorted(list(set([tup[0] for tup in tuples])))
    

def convolutionAll(spectrum):
    count = {}    
    allmasses = sorted(spectrum, reverse=True)
    for i in range(len(allmasses)):
        for j in range(i+1,len(allmasses)):
            diff = abs(allmasses[i]-allmasses[j])    
            if diff != 0:          
                if diff in count:
                    count[diff]+=1
                else:
                    count[diff]=1
    return sorted([(key, count[key]) for key in count], key=lambda tup: tup[1], reverse=True)


def pepNumberToString(pep_masses):
    seq = ""
    for m in pep_masses:
        seq+= getChar[m]
    return seq






lines = sys.stdin.read().splitlines()

#pep_string = lines[0]
#spectrum = [int(m) for m in lines[1].split(' ')]

# peptide = stringPepToIntList(pep_string)
# cy_score = score(peptide, spectrum, fn=cyclicSpectrumInt)
# print(cy_score)

# n = int(lines[0])
# spectrum = [int(m) for m in lines[1].split(' ')]
# final_peps = leaderboardCyclePeptideSequencing(spectrum, n, masses=range(57,201))
# print(' '.join(sorted(['-'.join(map(str,pep)) for pep in final_peps], reverse=True)))

# result = convolutionAll(spectrum)
# string = ""
# for tup in result:
#     string = ' '.join([string]+[str(tup[0])]*tup[1])
# print(string)

# m = int(lines[0])
# n = int(lines[1])
# spectrum = [int(m) for m in lines[2].split(' ')]
# final_peps = leaderboardCyclePeptideSequencing(spectrum, n, masses=convolution(spectrum,m))
# print(' '.join(sorted(['-'.join(map(str,pep)) for pep in final_peps], reverse=True)))

# N 114.04293
# Q 128.05858
# Y 163.06333
# V 99.06841
# K 128.09496
# L 113.08406
# G 57.02146
# F 147.06841
# P 97.05276
# W 186.07931
# F 147.06841

# masses = [getInt[c] for c in 'VKLFPWFNQY']

# spectrum = [float(m)-1 for m in lines[2].split(' ')]
# # spectrum = [int(round(x-0.5)) for x in spectrum]
# # spectrum += [int(round(x+0.5)) for x in spectrum]
# spectrum = [int(round((x- 1.007) * (1 - 0.0004522))) for x in spectrum]  
# spectrum = sorted(list(spectrum))
# m = int(lines[0])
# n = int(lines[1])

# #masses += convolution(spectrum, m)
# masses = sorted(list(set(masses)))

# final_peps = leaderboardCyclePeptideSequencing(spectrum, n, masses=masses)
# print('\n'.join(sorted([' '.join(map(str,map(int,pep))) for pep in final_peps], reverse=True)))

# spec = [0,97,97,129,129,194,203,226,226,258,323,323,323,355,403,452]
# pep_string = "PEEP"
# peptide = stringPepToIntList(pep_string)
# cy_score = score(peptide, spec, fn=linearSpectrumInt)
# print(cy_score)

spec = [0,57,118,179,236,240,301]
print(convolutionAll(spec))