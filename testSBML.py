import sbml_to_cot as sbml2cot
import cotlib


#Parsing file

species, reactionIDs, reactions = sbml2cot.parse_file("/Disco/COT/software/COT_analyzer/sbml_to_cot/example_models/BIOMD0000000497.xml")
#species, reactionIDs, reactions = sbml2cot.parse_file("/Disco/COT/software/COT_analyzer/sbml_to_cot/example_models/curated/BIOMD0000000087.xml")
#species, reactionIDs, reactions = sbml2cot.parse_file("/Disco/COT/OrgTools/examples/ecoli.xml")

print "Number of species:"+str(len(species))
print "Number of reactions:"+str(len(reactions))

#Creating species
for s in species:
    cotlib.species(s)

#Creating reactions
for i,r in enumerate(reactions):
    #print "["+str(i)+"]: "+str(r)
    cotlib.reaction(r)


print "M := \n"
cotlib.show_M()

print "\nR := \n"
cotlib.show_R()

#Removing outliers (all reactions with more reactants than the median)
#cotlib.filter_by_std(0)

#Creating reaction network (with all the elements in M)
element_indexes_in_RN = range(len(species))
RN = cotlib.renet(element_indexes_in_RN)
RN.add_inflow(cotlib.get_inflow_indexes())
'''

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

cotlib.reaction(" a -> b")
cotlib.reaction("2b -> a")

cotlib.reaction("  -> Hello + World")

RN = cotlib.renet(range(len(cotlib.M)))
'''

print "\nReaction Network:"
print "=====================\n"
print "Set X:"
RN.show_X()
print "\nSet R_x:"
RN.show_Rx()
print "\nIs closed?:"
print "\t"+str(RN.isClosed())
print "Is SSM?:"
print "\t"+str(RN.isSSM())
print "Is SM?:"
print "\t"+str(RN.isSM())
#
# The next lines where already tested. This process needs a handful
# of minutes (it solves 294 linear programming problems), so the result
# is already written. Note: it corresponds to "lipid [intracellular]"
#
#print "\nList of all overproduced species in X:"
#RN.showOverproduced()


G,tree=RN.reactant_closure(verbose=True,add_separables=False,render=True,save_graph=True,print_log=True)