from hasselib import *

S = set([0,3,5,7])
length = len(S)
P = powerset(S)
H = [None]*( length + 1 )

for i in range ( length + 1 ):
    H[i] = [j for j in P if len(j) == i]


#Print partition of powerset by cardinality of the subsets
'''
for h in H:
    print h
'''

for H_level in H:
    for elem in H_level:
        print elem.calc_greater_neighbours(H)
