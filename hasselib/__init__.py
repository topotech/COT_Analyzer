# *-* coding: utf-8 *-*

from itertools import chain, combinations

class HasseElem:
    def __init__(self,iterable):
        self.set = set(iterable)
        self.level = len(self.set)
        self.greater = None

        self.isClosed = False
        self.isSSM = False
        self.isSM = False
        self.isOrg = False

    def __repr__(self):
        return str(self.set)

    def __len__(self):
        return len(self.set)

    '''H(asse) list is needed to search neighbour'''
    def calc_greater_neighbours(self,H):
        ''' If this is the top element, break'''
        if self.level == len(H)-1:
            return self.greater

        self.greater = list()
        j = self.level + 1
        for h in H[j]:
            if (self.set == self.set & h.set):
                self.greater.append(h)
        return self.greater


def powerset(iterable):
    s = list(iterable)
    L = list( chain.from_iterable(set(combinations(s, r)) for r in range(len(s)+1)) )
    return [ HasseElem(i) for i in L ]