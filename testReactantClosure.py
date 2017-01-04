# *-* coding: utf-8 *-*

import cotlib
import copy


cotlib.reaction("a+c+d+e -> c+d")
cotlib.reaction("f+g -> a+g+e")
cotlib.reaction("f+e -> d+b")
cotlib.reaction("b+c+g -> f")
cotlib.reaction("h+i+k+a+c+f -> f+a")
cotlib.reaction("h+i+f -> f+a+d+r+w+t+b")
cotlib.reaction("h+a+i+f -> f+a+d+r+w+t+b+z")
cotlib.reaction("->a")

'''
cotlib.reaction("a+b->2b")
cotlib.reaction("a+b->2a")

cotlib.reaction("a+c->2c")
cotlib.reaction("a+c->2a")
'''
'''
cotlib.reaction("a+c+d+e -> c+d")
cotlib.reaction("f+g -> a+g+e")
cotlib.reaction("f+e -> d+b")
cotlib.reaction("b+c+g -> f")
cotlib.reaction("h+i+k+a+c+f -> f+a")
cotlib.reaction("h+i+f -> f+a+d+r+w+t+b")
cotlib.reaction("h+a+i+f -> f+a+d+r+w+t+b+z")
cotlib.reaction("z1 -> c+d")
cotlib.reaction("a1+c1+d1+e1 -> c1+d1")
cotlib.reaction("f1+g1 -> a1+g1+e1")
cotlib.reaction("f1+e1 -> d1+b1")
cotlib.reaction("b1+c1+g1 -> f1")
cotlib.reaction("h1+i1+k1+a1+c1+f1 -> f1+a1")
cotlib.reaction("h1+i1+f1 -> f1+a1+d1+r1+w1+t1+b1")
cotlib.reaction("h1+a1+i1+f1 -> f1+a1+d1+r1+w1+t1+b1+z1")
cotlib.reaction("f12+e12 -> d12+b12")
cotlib.reaction("b12+c12+g12 -> f12")
cotlib.reaction("h12+i12+k12+a12+c12+f12 -> f12+a12")
cotlib.reaction("h12+i12+f12 -> f12+a12+d12+r12+w12+t12+b12")
cotlib.reaction("h12+a12+i12+f12 -> f12+a12+d12+r12+w12+t12+b12+z12")
'''


#Displaying all species in M (global species set)
print "M := \n"
cotlib.show_M()

#Displaying all reactions in R (global reactions set)
print "\nR := \n"
cotlib.show_R()


RN1 = cotlib.renet(range(len(cotlib.M)))



'''

#Playing with RN1

print "\nReaction Network #1:"
print "=====================\n"
print "Set X:"
RN1.show_X()
print "\nSet R_x:"
RN1.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN1.isClosed())
print "Is SSM?:"
print "\t"+str(RN1.isSSM())
print "Is SM?:"
print "\t"+str(RN1.isSM())
print "\nList of all overproduced species in X:"
RN1.showOverproduced()

'''
print "\n"

RN1.reactant_closure(verbose = False, add_separables=True)
RN1.reactant_closure(verbose = False)