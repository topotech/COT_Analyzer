#import sbml_to_cot as sbml2cot
#import cotlib
import libsbml


#Parsing file

reader = libsbml.SBMLReader()
document = reader.readSBML("/Disco/COT/software/COT_analyzer/sbml_to_cot/example_models/curated/BIOMD0000000611.xml")
sbml_model = document.getModel()
if bool(sbml_model) == False:
    raise IOError("Can't extract a valid model from the path. Check if the path is correct and if it's a valid SBML file.")

print sbml_model.getNumSpecies()
print sbml_model.getNumReactions()
