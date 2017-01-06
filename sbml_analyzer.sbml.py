# *-* coding: utf-8 *-*

import cotlib
import sbml_to_cot as sbml2cot
import os
import datetime
import time

results_path="./results/"
input_text="./sbml_to_cot/example_models/files_to_read.txt"
input_path="./sbml_to_cot/example_models/curated/"

def final_writings(intial_time,results_path,file):
    termination_time = time.time()
    delta = (termination_time - intial_time)
    with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
        o.write("Elapsed time: "+str(datetime.timedelta(seconds=delta))+"\n")

def appendix(results_path,file):
    with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
        o.write("\nSpecies\n=======\n")
        for i, elem in enumerate(cotlib.M):
            o.write("["+str(i)+"]: "+str(elem)+"\n")
        o.write("\nReactions\n=========\n")
        for i,reac in enumerate(cotlib.R):
            o.write("["+str(i)+"]: "+str(reac)+"\n")
        o.write("\n")

with open(input_text,"r") as f:
    for line in f:
        flag_repeated_species = False
        file= line.rstrip("\n")
        print file
        intial_time = time.time()
        with open(results_path+file.rstrip(".xml")+".proc", "w") as o:
            o.write("Filename: "+file+"\n")
            o.write("File size: "+str(os.path.getsize(input_path+file))+" bytes\n")

        species, reactionIDs, reactions = sbml2cot.parse_file(input_path+file)
        with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
            o.write("Number of species:"+str(len(species))+"\n")
            o.write("Number of reactions:"+str(len(reactions))+"\n")

        if len(species) <4 or len(reactions)<3:
            with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
                o.write("Dismissed: File has less than 4 species or three reactions"+"\n")
        else:
            for s in species:
                try:
                    cotlib.species(s)
                except ValueError:
                    with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
                        o.write("Aborted because repeated species: "+str(s)+"\n")
                    flag_repeated_species= True
                    pass
            if not flag_repeated_species:
                for i,r in enumerate(reactions):
                    try:
                        cotlib.reaction(r)
                    except ValueError:
                        with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
                            o.write("Repeated reaction"+" ("+reactionIDs[i]+")"+": "+str(r)+"\n")
                        pass
                element_indexes_in_RN = range(len(species))
                RN = cotlib.renet(element_indexes_in_RN)
                RN.add_inflow(cotlib.get_inflow_indexes())
                try:
                    RN.reactant_closure(verbose=True,add_separables=False,save_graph=True,save_plot=True,print_log=True,sbml_filename=results_path+file.rstrip(".xml"),timeout=1500)
                except RuntimeError:
                    with open(results_path+file.rstrip(".xml")+".proc", "a") as o:
                        o.write("Timeout exceeded"+"\n")
                    pass
                del RN
        final_writings(intial_time,results_path,file)
        appendix(results_path,file)
        del cotlib.R[:]
        del cotlib.M[:]



'''
intial_time = time.time()
species, reactionIDs, reactions = sbml2cot.parse_file("")
print "Number of species:"+str(len(species))
print "Number of reactions:"+str(len(reactions))

os.path.getsize()



termination_time = time.time()
delta = termination_time - intial_time
string= "Elapsed time: "+str(datetime.timedelta(seconds=delta))
'''