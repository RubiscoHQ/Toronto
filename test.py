import pandas as pd

f = open('/Users/rubisco/Desktop/Toronto/Mutation_analyse/input/ExAC.all.hg19_multianno.new.vcf.hq.input.narm.txt')
VI = []
V = []
for line in f:
    ll = line.strip().split('\t')
    if ll[-1] == 'VI':
        VI.append(ll[5])
    elif ll[-1] == 'V':
        V.append(ll[5])


vi = pd.Series(VI, dtype=float)
print vi.mean(), 'VI'

v = pd.Series(V, dtype=float)
print v.mean(), 'V'

