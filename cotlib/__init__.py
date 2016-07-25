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


M = set()
R = set()


class species:
    myReactions = set()
    def __init__(self, name):
        if isinstance(name,str):
            self.name = name
        else:
            raise TypeError

class reaction:
    MyOwners = set()

class renet:
    X = set()
    R_x = set()
    S_x = np.array([])

    def isSM(self, S):
        #M should be np.array
        vars = np.shape(S)[1]
        eqs1 = np.shape(S)[0]
        eqs2 = vars
        I = np.identity(vars)
        A = S
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



