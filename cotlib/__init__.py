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
        return producedSet.issubset(set(self.X))

    #isSSM(): Return True if X is semi self-maintaining (False if not).
    def isSSM(self):
        producedSet = set()
        for r in self.R_x:
            aux = set([item[1] for item in r.products])
            producedSet = producedSet.union(aux)
        #if X is subset of the elements produced by R_x

        return set(self.X).issubset(producedSet)