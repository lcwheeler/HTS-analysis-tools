#!/usr/bin/python3

import pickle
import numpy as np
import jellyfish
import itertools as it
import argparse 

## calculate the medoid of sequences in a list using Damerau-Levensthein Distance as metric
def find_medoids(seqs):

    with open(seqs, "r") as seqs:
        seqs = [line.strip() for line in seqs.readlines()]

    internal_distances = {}
    
    #build a dictionary of all pairswise distances between sequences in a cluster
    for j, k in it.combinations(seqs, 2):
        internal_distances[j,k] = jellyfish.damerau_levenshtein_distance(j,k)
        
    #find mean and stdv of internal node distances
    mean = np.mean(list(internal_distances.values()))
    stdv = np.std(list(internal_distances.values()))    
    
    ## find the sum of distances for each sequence in the cluster
    dist_sums = {}
    for seq in seqs:     
        dists = []
        for key, value in internal_distances.items():
            if seq in key:
                dists.append(value)
            else:
                pass
        dist_sums[seq] = sum(dists)/len(dists)

    # find the sequence/s in the cluster with minimum sum of distances and save only one if degenerate
    medval = min(dist_sums.values())
    medoids = {}

    for key, value in dist_sums.items():
        if value == medval:
            medoids[key] = value
    return medoids, ("mean", mean), ("stdv", stdv) 

#add argparse main function to take command line arguments (input file names)
def main():
    """
        Run medoid calculation on input file.
    """
                             
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("seqs", help="Input input filename")
    args = parser.parse_args()


    #find the medoids, as well as mean+stdv of distance for the strings in the input file and 
    #outputs as a pickled dictionary and text file respectively
    medoids, mean, stdv = find_medoids(args.seqs)
    pickle.dump(medoids, open( "medoids.p", "wb"))
    np.savetxt("distance-mean_stdv.txt", np.c_[mean, stdv], fmt="%s", delimiter="\t", header="mean & stdv of distance")                  
if __name__ == "__main__":
    main()
