# *-* coding: utf-8 *-*
"""
COTlib
======

Python package for complex systems analysis through Chemical Organization
Theory.

This package provides three classes needed to build a reaction network:
    1.- species:    Species
    2.- reaction:   Reactions
    3.- renet:      Reaction Networks

This package also provides two global sets containing all the species
and reactions created, called M and R respectively.

"""

import cotlib
import cvxopt as co
import numpy as np
import re
import hasselib

"""
Global variables and functions
==============================
"""

# Global sets
M = list()
R = list()

# Display all elements in M
def show_M():
    for i, elem in enumerate(M):
        print "["+str(i)+"]: "+str(elem)

# Display all reactions in R
def show_R():
    for i,reac in enumerate(R):
        print "["+str(i)+"]: "+str(reac)

# Remove species from M
def rm_species(name):
    if not(isinstance(name,str)):
        raise TypeError("rm_species() expects a string")
    for elem in M:
        if name == elem.name:
            M.remove(elem)
            return True
    return False

# Verify of species exists in M: return boolean and index (-1 if false)
def in_M(name):
    if not(isinstance(name,str)):
        raise TypeError("in_M() expects a string")
    for i,elem in enumerate(M):
        if name == elem.name:
            return True, i
    return False, -1

"""
Classes and their methods
=========================
"""


# Class species
class species:
    def __init__(self, name):
        if not(isinstance(name,str)):
            raise TypeError("species constructor expects a string")

        for elem in M:
            if name == elem.name:
                raise ValueError("there can't be two elements with the same name")
        self.myReactions = list()
        self.name = name
        M.append(self)

    def __repr__(self):
        return self.name

    def ins_reaction(self,myreaction):
        if not(isinstance(myreaction,reaction)):
            raise TypeError("ins_reaction() expects an instance of reaction class")
        self.myReactions.append(myreaction)


# Class reaction
class reaction:

    def __init__(self, repr):
        if not(isinstance(repr,str)):
            raise TypeError("reaction constructor expects a string")

        self.repr = repr.replace(" ", "")

        pattern = re.compile(r"(.*)(?:->|→|⇒)(.*)")
        match = re.match(pattern, self.repr)

        if not(match):
            raise ValueError("The reaction string was badly formatted")
        if len(match.groups()) != 2:
            raise ValueError("The reaction string was badly formatted")

        str_reactants = match.groups()[0].split("+")
        str_products  = match.groups()[1].split("+")

        self.reactants = set()
        self.products= set()
        pattern = re.compile(r"([0-9]*)(\w*)")

        if len(match.groups()[0]) > 0:
            for string in str_reactants:
                match = re.match(pattern, string)
                if not(match):
                    raise ValueError("The reactants were badly formatted")
                if len(match.groups()[1]) <= 0:
                    raise ValueError("The reactants were badly formatted")

                if len(match.groups()[0]) == 0:
                    self.reactants.add((1,match.groups()[1]))
                else:
                    self.reactants.add((int(match.groups()[0]),match.groups()[1]))

        if len(match.groups()[1]) > 0:
            for string in str_products:
                match = re.match(pattern, string)
                if not(match):
                    raise ValueError("The products were badly formatted")
                if len(match.groups()[1]) <= 0:
                    raise ValueError("The products were badly formatted")

                if len(match.groups()[0]) == 0:
                    self.products.add((1,match.groups()[1]))
                else:
                    self.products.add((int(match.groups()[0]),match.groups()[1]))

        reac_list = [item[1] for item in self.reactants]
        prod_list = [item[1] for item in self.products]

        if len(reac_list) != len(set(reac_list)):
            raise ValueError("The reactants were badly formatted: there are repeated species")

        if len(prod_list) != len(set(prod_list)):
            raise ValueError("The products were badly formatted: there are repeated species")

        for reac in R:
            if self.reactants == reac.reactants and self.products == reac.products:
                raise ValueError("there can't be two equal reactions")

        for tup in self.reactants.union(self.products):
            boolean,index = in_M(tup[1])
            if boolean:
                M[index].ins_reaction(self)
            else:
                species(tup[1])
                boolean,index = in_M(tup[1])
                M[index].ins_reaction(self)

        R.append(self)

    def __repr__(self):
        return self.repr


# Class reaction network
class renet:
    S_x = np.array([])

    def __init__(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("renet constructor expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("renet constructor expects a list of integers")

        self.X = list()
        self.X_names = set()
        self.R_x = list()

        for i in indexes:
            self.X.append(M[i])
            self.X_names.add(M[i].name)

        for elem in self.X:
            for reaction in elem.myReactions:
                self.ins_reaction(reaction)

        self.cln_Rx()
        self.createA()
        self.createB()

        self.S = self.B - self.A


    # This let the user to add reactions with empty reactants outside X
    def add_species(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("add_species() expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("add_species() expects a list of integers")


        for i in indexes:
            if not (M[i].name in self.X_names):
                self.X.append(M[i])
                self.X_names.add(M[i].name)

        for elem in self.X:
            for reaction in elem.myReactions:
                self.ins_reaction(reaction)

        self.cln_Rx()
        self.createA()
        self.createB()

        self.S = self.B - self.A


    # This let the user to add reactions with empty reactants outside X
    def add_inflow(self, indexes):
        if not(isinstance(indexes,list)):
            raise TypeError("add_inflow() expects a list of integers")
        if not (all(isinstance(x, int) for x in indexes)):
            raise TypeError("add_inflow() expects a list of integers")

        for i in indexes:
            if len(R[i].reactants) == 0:
                self.ins_reaction(R[i])
            else:
                raise ValueError("Malformed input. All reactant sets need to be empty.")

        self.cln_Rx()
        self.createA()
        self.createB()

        self.S = self.B - self.A

    #cln_Rx(): Clean reaction: just leave the reactions which can be fired by X
    def cln_Rx(self):
        aux = list()
        for myreaction in self.R_x:
            reactants = set([item[1] for item in myreaction.reactants])
            if not (reactants.issubset(self.X_names)):
                aux.append(myreaction)
        for item in aux:
            self.R_x.remove(item)



    def ins_reaction(self,r):
        for reac in self.R_x:
            if (r.reactants == reac.reactants) and (r.products == reac.products):
                return False
        self.R_x.append(r)
        return True

    #createA(): Create reactants matrix
    def createA(self):

        self.A = np.zeros(  (len(self.X) , len(self.R_x) ) )

        for i,r in enumerate(self.R_x):

            for j,e in enumerate(self.X):

                for m in r.reactants:
                    if m[1] == e.name:
                        self.A[j][i] = m[0]


    #createB(): Create products matrix
    def createB(self):

        self.B = np.zeros(  (len(self.X) , len(self.R_x) ) )

        for i,r in enumerate(self.R_x):

            for j,e in enumerate(self.X):

                for m in r.products:
                    if m[1] == e.name:
                        self.B[j][i] = m[0]



    #close(): Turn X into it's closure (add all the species produced by R_x into X)
    def close(self):
        if not(self.isClosed()):
            indexes = list()
            producedSet = set()

            for r in self.R_x:
                aux = set([item[1] for item in r.products])
                producedSet = producedSet.union(aux)

            for name in producedSet:
                indexes.append( in_M(name)[1] )

            #Reinitializing the reaction network
            for i in indexes:
                if not (M[i].name in self.X_names):
                    self.X.append(M[i])
                    self.X_names.add(M[i].name)

            for elem in self.X:
                for reaction in elem.myReactions:
                    self.ins_reaction(reaction)

            self.cln_Rx()
            self.createA()
            self.createB()

            self.S = self.B - self.A



    #show_x(): Display a list with all elements in X and its indexes
    def show_X(self):
        for i, elem in enumerate(self.X):
            print "["+str(i)+"]: "+str(elem)

    #show_Rx(): Display all reactions in R_x
    def show_Rx(self):
        for i,reac in enumerate(self.R_x):
            print "["+str(i)+"]: "+str(reac)

    #isSM(): Verify if this reaction network is overproduced
    def isSM(self):
        #if X is empty, then it's SM
        if len(self.X) == 0:
            return True

        #vars is the numbers of columns of S. Represents
        # the  number of components of the unknow flux vector 'v'
        vars = np.shape(self.S)[1]

        #eqs1 is the numbers of rows of S. Represents
        # the number of inequations introduced in the LP problem.
        eqs1 = np.shape(self.S)[0]

        #eqs2 is the numbers of columns of S. Represents
        # the lower bound inequations introduced in the LP problem.
        eqs2 = vars

        I = np.identity(vars)
        A = self.S
        A = np.concatenate((A, I), axis=0)
        A = -A
        b0 =  np.zeros(eqs1)
        b1 = -np.ones(eqs2)
        b = np.concatenate((b0, b1), axis=0)

        c = np.ones(vars)

        A = co.matrix(A)
        b = co.matrix(b)
        c = co.matrix(c)

        co.solvers.options['show_progress']=False
        sol=co.solvers.lp(c,A,b)

        return bool(sol['x'])



    #isOverproduced(i): Verify if species m_i is overproduced. Use show_X() to list all indexes.
    def isOverproduced(self,i):

        if i >= len(self.X) or i < 0:
            raise ValueError("Index out of bound. A non-negative integer less than len(X) expected.")

        #vars is the numbers of columns of S. Represents
        # the  number of components of the unknow flux vector 'v'
        vars = np.shape(self.S)[1]

        #eqs1 is the numbers of rows of S. Represents
        # the number of inequations introduced in the LP problem.
        eqs1 = np.shape(self.S)[0]

        #eqs2 is the numbers of columns of S. Represents
        # the lower bound inequations introduced in the LP problem.
        eqs2 = vars

        I = np.identity(vars)
        A = self.S
        A = np.concatenate((A, I), axis=0)
        A = -A

        b =  np.zeros(eqs1+eqs2)
        b[i] = -1

        c = np.ones(vars)

        A = co.matrix(A)
        b = co.matrix(b)
        c = co.matrix(c)

        co.solvers.options['show_progress']=False
        sol=co.solvers.lp(c,A,b)

        return bool(sol['x'])

    #show_x(): Display a list with all overproduced elements in X
    def showOverproduced(self):
        boolean = True
        for i, elem in enumerate(self.X):
            if(self.isOverproduced(i)):
                print "["+str(i)+"]: "+str(elem)
                boolean = False

        if (boolean): print "None"


    #isClosed(): Return True if X is closed (False if not).
    def isClosed(self):
        producedSet = set()
        for r in self.R_x:
            aux = set([item[1] for item in r.products])
            producedSet = producedSet.union(aux)
        #if the elements produced by R_x is subset of X
        return producedSet.issubset(set(self.X_names))

    #isSSM(): Return True if X is semi self-maintaining (False if not).
    def isSSM(self):
        producedSet = set()
        for r in self.R_x:
            aux = set([item[1] for item in r.products])
            producedSet = producedSet.union(aux)
        #if X is subset of the elements produced by R_x

        return set(self.X_names).issubset(producedSet)

    # hierarchy_noDecomp() Returns the reference to Hasse graph and a list of organizations.
    def hierarchy_noDecomp(self):

        # Save global indexes to inflow (empty) reactions.
        R_inflow = list()
        for i,my_reac in enumerate(self.R_x):
            if len(my_reac.reactants) == 0:
                for j,global_reac in enumerate(R):
                    if my_reac.reactants == global_reac.reactants and my_reac.products == global_reac.products:
                        R_inflow.append(j)

        # Defining set with global indexes
        # (a.k.a. not using the local indexes of this network)
        S = []
        for elem in self.X:
            S.append(in_M(elem.name)[1])
        S = set(S)


        # Constructing the graph (adding the reference
        # to the neighbours of each element in the list H.
        length = len(S)
        P = hasselib.powerset(S)
        H = [None]*( length + 1 ) #Defining size of list according
                                  # to the height of Hasse Diag.
        for i in range ( length + 1 ):
            H[i] = [j for j in P if len(j) == i]
        for H_level in H:
            for elem in H_level:
                elem.calc_neighbours(H)

        '''
        Visit_queue = list([H[0][0]]) #Enqueuing the first (empty) subset
        '''
        Visit_queue = P

        org_list = []
        while (Visit_queue):
            curr_elem = Visit_queue.pop(0)
            if(curr_elem.isVisited):
                continue
            else:
                curr_elem.isVisited = True

            curr_renet = renet(list(curr_elem.set))
            curr_renet.add_inflow(R_inflow)

            if(curr_renet.isClosed()):
                curr_elem.isClosed = True
                curr_elem.closure = curr_elem
            else:
                del curr_renet
                self.seek_closure(curr_elem,length,R_inflow)
                continue

            if(curr_renet.isSSM()):
                curr_elem.isSSM = True
            else:
                del curr_renet
                curr_elem.isSSM = False
                curr_elem.isSM = False
                curr_elem.isOrg = False
                continue

            if(curr_renet.isSM()):
                curr_elem.isSM = True
            else:
                del curr_renet
                curr_elem.isSM = False
                curr_elem.isOrg = False
                continue

            if( curr_elem.isClosed and curr_elem.isSSM and curr_elem.isSM ):
                curr_elem.isOrg = True
                del curr_renet
                org_list.append(curr_elem)


        return H, org_list
    # Search for minimal closure for NON-CLOSED sets
    def seek_closure(self, hasse_element, max_cardinality, R_inflow):
        height = len(hasse_element.set) # 'Level' of initial element in Hasse diagram by cardinality
        curr_level = height+1 #Setting the next level as the first to vist
        curr_greater = hasse_element.greater
        FoundIt=False
        FoundClosure = None
        breadcrumb = list([hasse_element])

        while (curr_level <= max_cardinality):   # Si aún no llega al tope
            next_greater = set()                 # Crea una lista auxiliar
            if curr_greater != None:
                for e in curr_greater:               # Recorre todos los vecinos próximos
                    breadcrumb.append(e)
                    if e.greater != None:
                        next_greater.update(set(e.greater))  # Va creando la nueva lista de vecinos
                    if(e.isVisited):
                        continue
                    else:
                        e.isVisited = True

                    if FoundIt:
                        e.isClosed = False

                    curr_renet = renet(list(e.set))  # Instancia red reacciones
                    curr_renet.add_inflow(R_inflow)  #Agrega el inflow

                    if(curr_renet.isClosed()):       #Si es cerrado, flagear
                        FoundIt = True
                        e.isVisited = False
                        e.isClosed = True
                        e.closure = e
                        FoundClosure = e
                    else:
                        e.isClosed = False
                    del curr_renet                   #Limpiar referencia

                if FoundIt:
                    break
                curr_level += 1                     # Pasa al siguiente nivel
                curr_greater = list(next_greater)   # El auxiliar pasa a prinicipal

        if FoundClosure != None:
            for e in breadcrumb:
                if e.set.issubset(FoundClosure.set):
                    e.closure = FoundClosure

    # Test simple connection between species from Self (or 'X') with m_i
    def IsConnected (self,  i):
        m = M[i]
        sR_x = set(self.R_x)
        sR = set(R)
        sR_unused =  sR - sR_x
        sR_m = set(m.myReactions)
        sR_xUm = sR_x.union(sR_m)
        sxUm = self.X_names.union(m.name)

        aux = list()

        for myreaction in sR_xUm:
            reactants = set([item[1] for item in myreaction.reactants])
            if not (reactants.issubset(sxUm)):
                aux.append(myreaction)
        for item in aux:
            sR_xUm.remove(item)

        sR_new = sR_xUm - sR_x

        for myreaction in sR_new:
            reactants = set([item[1] for item in myreaction.reactants])
            products = set([item[1] for item in myreaction.products])

        if ( reactants.union(products).isdisjoint(self.X_names) ):
            return False
        else:
            return True