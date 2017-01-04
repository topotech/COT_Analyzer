"""
sbml_to_cot
======

Python package to turn a valid SBML input into a cotlib format

"""

import libsbml
import cotlib


#Returns lists with species IDs, reaction IDs and explicit reaction strings
def parse_file(path_to_sbml_file):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(path_to_sbml_file)
    sbml_model = document.getModel()
    if bool(sbml_model) == False:
        raise IOError("Can't extract a valid model from the path. Check if the path is correct and if it's a valid SBML file.")

    speciesIDs = [s.getId() for s in sbml_model.getListOfSpecies()] #Identifiers strings for each species.
    reactionIDs = [r.getId() for r in sbml_model.getListOfReactions()] #Identifiers strings for each reaction.
    reactionFormatted = [] #Explicit representation of each reaction in cotlib format.

    for reaction_index, reaction in enumerate(sbml_model.getListOfReactions()):

        #getting minimum

        reaction_string = ""
        curr_reaction_input = cotlib.reaction_input()

        if len(reaction.getListOfReactants()) != 0:
            last_reactant = reaction.getListOfReactants()[-1]
            for r in reaction.getListOfReactants():
                '''
                if reaction_index==282:
                    print str(r.getStoichiometry())+": "+str(r.getSpecies())
                '''
                #reaction_string += str(int(r.getStoichiometry()))+r.getSpecies()
                #print type(r.getStoichiometry())

                reaction_string += str(float(r.getStoichiometry()))+r.getSpecies()
                curr_reaction_input.reactants.add((float(r.getStoichiometry()),r.getSpecies()))
                if r != last_reactant:
                    reaction_string += "+"

        reaction_string += "->"
        #if reaction_index==282: print "A"

        if len(reaction.getListOfProducts()) != 0:
            last_product = reaction.getListOfProducts()[-1]
            for p in reaction.getListOfProducts():
                '''
                if reaction_index==282:
                    print str(p.getStoichiometry())+": "+str(p.getSpecies())
                '''
                #reaction_string += str(int(p.getStoichiometry()))+p.getSpecies()
                reaction_string += str(float(p.getStoichiometry()))+p.getSpecies()
                curr_reaction_input.products.add((float(p.getStoichiometry()),p.getSpecies()))
                if p != last_product:
                    reaction_string += "+"

        curr_reaction_input.repr = reaction_string
        reactionFormatted.append(curr_reaction_input)

    return speciesIDs, reactionIDs, reactionFormatted
