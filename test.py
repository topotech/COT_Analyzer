# *-* coding: utf-8 *-*

import cotlib
import cvxopt as co
import numpy as np

def isSM(S):
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



'''
Ejemplo ificc
-----
S=     [ 1, -1,  0,  0]
       [ 0,  1, -1,  1]
       [ 0,  0,  1, -2]


'''

'''
Ejemplo correo
==============

S = [-1   1]
....[ 1  -2]

'''


S0 = np.array([[-1  , 1], [ 1 , -2]])
S1= np.array([[ 1, -1,  0,  0],
       [ 0,  1, -1,  1],
       [ 0,  0,  1, -2]])


print isSM(S0)
print isSM(S1)



c = cotlib.species("H2O")

print len(c.myReactions)

c.myReactions.add("C")
c.myReactions.add("D")

print c.myReactions
print len(c.myReactions)


