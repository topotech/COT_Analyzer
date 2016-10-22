# *-* coding: utf-8 *-*

from itertools import chain, combinations

class HasseElem:
    def __init__(self,iterable):
        self.set = set(iterable)
        self.level = len(self.set)
        self.greater = None #Sets which contain current set (reference)
        self.lesser = [] #Sets contained in current set (reference)
        self.closure = None #The minimal closure (reference)


        self.isClosed = None
        self.isSSM = None
        self.isSM = None
        self.isOrg = None
        self.isVisited = False


    def __repr__(self):
        return str(self.set)

    def __len__(self):
        return len(self.set)

    '''H(asse) list is needed to search neighbour'''
    def calc_neighbours(self,H):
        ''' If this is the top element, break'''
        if self.level == len(H)-1:
            return self.greater

        self.greater = list()
        j = self.level + 1
        for h in H[j]:
            if (self.set == self.set & h.set):
                self.greater.append(h)
                h.lesser.append(self)
        if len(self.lesser) == 0:
            self.lesser = None
        return self.greater


def powerset(iterable):
    s = list(iterable)
    L = list( chain.from_iterable(set(combinations(s, r)) for r in range(len(s)+1)) )
    return [ HasseElem(i) for i in L ]