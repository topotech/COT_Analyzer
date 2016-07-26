# *-* coding: utf-8 *-*

import cotlib


#Creating species. Note: it's case sensitive.
cotlib.species("A")
cotlib.species("B")
cotlib.species("C")

cotlib.species("a")
cotlib.species("b")


#Creating reactions
cotlib.reaction("-> A")
cotlib.reaction("A  -> B")
cotlib.reaction("B  -> C")
cotlib.reaction("2 C  -> B")

cotlib.reaction(" a → b")
cotlib.reaction("2b → a")

cotlib.reaction("  -> Hello + World")

#Displaying all species in M (global species set)
print "M := \n"
cotlib.show_M()

#Displaying all reactions in R (global reactions set)
print "\nR := \n"
cotlib.show_R()


RN1 = cotlib.renet([0,1,2])
RN2 = cotlib.renet([3,4])

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

print"\n\n"

print "\nReaction Network #2:"
print "=====================\n"
print "Set X:"
RN2.show_X()
print "\nSet R_x:"
RN2.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN2.isClosed())
print "Is SSM?:"
print "\t"+str(RN2.isSSM())
print "Is SM?:"
print "\t"+str(RN2.isSM())
print "\nList of all overproduced species in X:"
RN2.showOverproduced()