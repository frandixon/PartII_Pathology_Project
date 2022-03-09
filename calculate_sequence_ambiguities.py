#Counts the number of ambiguities in each sequence in an alignment
#Outputs a csv file containing the ambiguity content of each sequence and an alignment filtered to keep only
#sequences with less than a given proportion of ambiguities, default 5%
#Gaps are treated as real deletions so are not treated as ambiguities
#Expects an alignment that has been created with pangolin so is padded with Ns at the 5' and 3' end
#Ambiguities are only counted in the coding region between positions 266 and 29674
#Output csv file contains: sequence name, number ambiguities, proportion ambiguities

import argparse
from Bio import SeqIO
from collections import Counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "fasta alignment file")
    parser.add_argument("-o", help = "Name of output csv file containing ambiguity stats")
    args = parser.parse_args()

    #Import the alignment
    #alignment = AlignIO.read(args.a, "fasta")

    outFile = open(args.o, "w")
    outFile.write("sequence,number_ambiguities,proportion_ambiguities\n")

    #Used to check if the sequence is the first to be analysed
    fr = True

    for record in SeqIO.parse(args.a, "fasta"):
        if fr:
            #Identify alignment length from first sequence
            sL = float(len(record.seq))
            fr = False
        
        nC = Counter(record.seq)

        #Number of ambiguities
        nA = float(sL - (nC["A"] + nC["C"] + nC["G"] + nC["T"]))

        #Proportion of ambiguities
        pA = (nA/sL) * float(100)

        outFile.write(record.id + "," + str(nA) + "," + str(pA) + "\n")

    outFile.close()
