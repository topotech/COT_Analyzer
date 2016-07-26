# *-* coding: utf-8 *-*

import cotlib
import numpy as np

S0 = np.array([[-1  , 1], [ 1 , -2]])
S1= np.array([[ 1, -1,  0,  0],
       [ 0,  1, -1,  1],
       [ 0,  0,  1, -2]])


e = list()
e.append(cotlib.species("H2O"))
e.append(cotlib.species("CO2"))
e.append(cotlib.species("NaCL"))
e.append(cotlib.species("CHO3"))



cotlib.rm_species("NaCL")

print cotlib.in_M("CO2")
print cotlib.in_M("H2O")



cotlib.show_M()


A = cotlib.reaction('2 Arbol+ 3 H2O + 5 velas -> Hola + Mundo')
A = cotlib.reaction('CO2 + H2O -> H2CO3')
A = cotlib.reaction(' -> H2CO3')
print ""
cotlib.show_M()
print ""

A=cotlib.renet([0])
B=cotlib.renet([0,1])
C=cotlib.renet([0,1,4,6])

D=cotlib.renet([7])


print A.R_x
print B.R_x
print C.R_x
print D.R_x