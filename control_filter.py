#Script for filtering strings from files based on Levensthein-Damerau distance to strings in a control file -Luke Wheeler

#import dependencies
import numpy as np
import jellyfish
import glob


## define the name of a control file to load in as reference and the file extension of files to be filtered
CONTROL = ""
extension = ""
## set the distance threshold to use as a cutoff (Levenshtein-Damerau string distance metric)
cutoff = 6
## define a function to compare distance of sequences in input files to the control reference and remove those that
## are too similar

def subtract(filename):
    
    # import control file as list
    with open(CONTROL, "r") as control:
        control_list = [line.strip() for line in control.readlines()]

    # import other files to subtract from
    with open(filename, "r") as peptides:
        peptides_list = [line.strip() for line in peptides.readlines()]

    # open empty lists and then append sequences based on string distance metric
    different = []
    close = []

    # tests to see if each peptide is close to all the peptides in control file, based on the cutoff
    for i in peptides_list:
        for j in control_list:
            metric = jellyfish.damerau_levenshtein_distance(str(i), str(j))
            if metric < cutoff:
		if i not in close: #this is a modification made to reduce memory use, output file size, etc. only append unique entries
                    close.append(i)
                
    # screens out any peptides that had a hit in the control        
    for i in peptides_list:
        if i not in close:
            different.append(i)
        
    # save the filtered lists and the removed hits as text files if the lists have content       
    if len(different) > 0:
        np.savetxt(filename+"_.controlsubtracted", different, fmt="%s", delimiter="\n")
    if len(close) > 0:        
        np.savetxt(filename+"_.hitscontrol", close, fmt="%s", delimiter="\n")
        
    return different, close


for file in glob.glob(extension)[0:]:
    subtract(file)


