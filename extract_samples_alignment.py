#Extracts samples from an alignment if their sample name is in a given list

import argparse
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "Fasta alignment from which samples will be extracted")
    parser.add_argument("-s", help = "Samples to be extracted, one per line with no header")
    parser.add_argument("-o", help = "Name of output alignment")
    args = parser.parse_args()

    #Extract the samples to keep to a list
    samples = list()
    with open(args.s) as fileobject:
        for line in fileobject:
            samples.append(line.strip())
    
    outFile = open(args.o, "w")

    #Iterate through the samples in the alignment, check if they are in the list to keep and write if so
    for record in SeqIO.parse(args.a, "fasta"):
        if record.id in samples:
            outFile.write(">" + record.id + "\n" + str(record.seq) + "\n")

    outFile.close()
