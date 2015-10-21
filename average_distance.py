# A script for calculating and plotting the average Levensthein-Damerau distance between all the sequences in two input files. -Luke Wheeler

# import necessary modules (need to give matplotlib a place to direct the output i.e. ".use("Pdf")")
import numpy as np
import matplotlib
matplotlib.use("Pdf")
import matplotlib.pyplot as plt
import glob
import jellyfish
import argparse


# define function tp calculate average Levenshtein-Damerau distance between two sets of peptide sequences.

def find_distances(file1, file2):
    
    #open list to dump calculated distances into
    distances = []
    # import files to compare
    with open(file1, "r") as file1:
        file1_list = [line.strip() for line in file1.readlines()]

    with open(file2, "r") as file2:
        file2_list = [line.strip() for line in file2.readlines()]

    for i in file1_list:
        for j in file2_list:
            distances.append(jellyfish.damerau_levenshtein_distance(i, j))
            
    mean = np.mean(distances)
    stdv = np.std(distances)
    
    return distances, mean, stdv

# define function that generates a histogram of the calculated distances
def plot_distances(distances, file1, file2):
    n, bins, patches = plt.hist(distances, 50, histtype='stepfilled',color='b', alpha=.75)
    plt.title('Histogram of Average Levenshtein-Damerau Distance', fontsize=10)
    plt.xlabel('Levensthein Distance', fontsize=14)
    plt.ylabel('Counts', fontsize=14)
    plt.xlim(0, 15)
    ymax = max(distances) + (5/100)*max(distances)
    plt.ylim(0, ymax)
    
    
    
#add argparse main function to take command line arguments (input file names)
def main():
    """
        Run distance calculation on input files.
    """
                             
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("file1", help="Input 1st input filename")
    parser.add_argument("file2", help="Input 2nd input filename")
    args = parser.parse_args()


    #do the distance calculation
    distances, mean, stdv = find_distances(args.file1, args.file2)
    np.savetxt(str(args.file1)+"_"+str(args.file2)+"_comparison"+".txt", np.c_[mean, stdv], fmt="%s", delimiter=",", header="mean, stdv")
    plot = plot_distances(distances,args.file1, args.file2)
    plt.savefig((str(args.file1)+"_"+str(args.file2)+"distance_histrogram.pdf")) 
                       
if __name__ == "__main__":
    main()
    
