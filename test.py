# *-* coding: utf-8 *-*

import cotlib
import copy


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


#RN1 = cotlib.renet([0,1,2])
RN2 = cotlib.renet([1])


'''

#Playing with RN1

print "\nReaction Network #1:"
print "=====================\n"
print "Set X:"
RN1.show_X()
print "\nSet R_x:"
RN1.show_Rx()

RN1.close()
'''
print "______"
RN2.show_X()
print "----"
RN2.show_Rx()
print RN2.isClosed()

'''
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

RN3= copy.deepcopy(RN1)
RN1.add_inflow([0,6])

print"\n\n"
print "\nReaction Network #1` (w/inflow):"
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

RN1.close()

print"\n\n"
print "\nReaction Network #1`` (closed):"
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
print "\nReaction Network #3`` (copia):"
print "=====================\n"
print "Set X:"
RN3.show_X()
print "\nSet R_x:"
RN3.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN3.isClosed())
print "Is SSM?:"
print "\t"+str(RN3.isSSM())
print "Is SM?:"
print "\t"+str(RN3.isSM())
print "\nList of all overproduced species in X:"
RN3.showOverproduced()

RN0 = cotlib.renet([])
RN0.add_inflow([0,6])
RN0.add_species([1])
RN0.close()

print"\n\n"
print "\nReaction Network #0 (empty):"
print "=====================\n"
print "Set X:"
RN0.show_X()
print "\nSet R_x:"
RN0.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN0.isClosed())
print "Is SSM?:"
print "\t"+str(RN0.isSSM())
print "Is SM?:"
print "\t"+str(RN0.isSM())
print "\nList of all overproduced species in X:"
RN0.showOverproduced()





#Creating reactions
cotlib.reaction("x -> y+z")
cotlib.reaction("x+y  -> x")
cotlib.species("w")

RN3 = cotlib.renet([7,8,9,10])

print"\n\n"
print "\nReaction Network #3:"
print "=====================\n"
print "Set X:"
RN3.show_X()
print "\nSet R_x:"
RN3.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN3.isClosed())
print "Is SSM?:"
print "\t"+str(RN3.isSSM())
print "Is SM?:"
print "\t"+str(RN3.isSM())
print "\nList of all overproduced species in X:"
RN3.showOverproduced()


'''