#!/bin/python3
from Bio import SeqIO
import argparse
import subprocess, os

#Part 1: Building the command line interface

def commandline():
    print("Welcome to the commandline I've been told to make. I hope you are welcomed and find this interface intuitive to use.")
    parser = argparse.ArgumentParser(prog = "BA_ICA_B188292", formatter_class=argparse.RawDescriptionHelpFormatter, description = """\
Welcome to the BLAST101 Search Interface
Analyze your sequence using:
- BLAST
- Smith-Waterman
- Statistical analysis""")

    #Database Input
    parser.add_argument('--database', required = True, help = "The database for this project")

    #Query Sequence Input
    parser.add_argument('--sequence', required = True, help = "The sequence you would like to BLAST")

    #Algorithm/Analysis Input
    parser.add_argument('--algorithm', default = 'blast', choices=['blast', 'Smith-Waterman', 'statistical_analysis'], help = "If you would like to run blast or smith-waterman, select the relevant option. If you would like to conduct statistical analysis, select that option.")

    args = parser.parse_args()

    #Checking we have the input
    if not os.path.exists(args.sequence):
        raise FileNotFoundError(f"Could not find {args.sequence}")

    if not os.path.exists(args.database):
        raise FileNotFoundError(f"Could not identify {args.database}")


    if args.sequence.lower().endswith(".fasta", ".fa"):
        print("Input okay, proceeding...")

        prot_alph = set("ACDEFGHIKLMNPQRSTVWY")
        dna_alph = set("ACGT")

        record = next(SeqIO.parse(args.sequence, "fasta"))
        sequence = str(record.seq).upper()

        protein_count = sum(1 for char in sequence if char in prot_alph)
        dna_count = sum(1 for char in sequence if char in dna_alph)

        #checking the residues
        if protein_count > dna_count:
            print("This is a protein sequence.\nAccepted.")
        else:
            raise ValueError("Only protein sequences accepted.")

    else:
        raise ValueError(f"Only .fasta files are supported.")


    #Message to the user
    print(f"Running {args.sequence} against {args.database} using {args.algorithm}")

    #Running the appropriate analysis
    if args.algorithm == "blast":
        output_file = f"{database}_{query}_output.out"
        subprocess.run(f"blastx -db {database} -query {query} > {output_file}", shell = True)
        print(f"Output has been saved to {output_file}")

    elif args.algorithm == "smith-waterman":
        print("SW code goes here")

    elif args.algorithm == "statistical_analysis":
        print("SW code goes here")

    return

commandline()

##############
# References #
##############
#The following handbook was used in creating the commandline:
#https://docs.python.org/3/library/argparse.html
#https://docs.python.org/3/tutorial/datastructures.html#sets

