#!/usr/bin/python3.4

#Parses a fastq file and translates DNA sequences to protein sequences, saves output as a text file. Indicates 
#ambiguous residues with "X".
import numpy as np
import glob
import re


def translate(filename):
    sequence = np.loadtxt(filename, dtype = bytes, comments = "@", delimiter = "\n").astype(str)
    
    GENCODE = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
    
    translated = []
    regex = re.compile('[^ACGTN]')

    for s in np.arange(0,len(sequence)):
        if len(regex.findall(sequence[s])) != 0:
            
            try:
                return ''.join([GENCODE[sequence[s][3*i:3*i+3]]
                                for i in range(len(sequence[s])//3)])
            except KeyError:
                out = []
            
                for i in range(len(sequence[s])//3):
                    try:
                        out.append(GENCODE[sequence[s][3*i:3*i+3]])
                    except KeyError:
                        out.append("X")
            translated.append("".join(out))

    print(translated)
    np.savetxt(filename+".trans", translated, fmt = "%s", delimiter="/n")
        
for file in glob.glob("*.fastq")[0:]:
    translate(file)
