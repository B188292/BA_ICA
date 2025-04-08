#!/bin/python3
from Bio import SeqIO
import argparse, subprocess, os
import blosum #<- used for creating the Blosum62 matrix
import smith_waterman_p as sw_p
import blast_101_search as blast #<- runs the main search for plast
import programme_settings
import process_fasta_file as pff #<- used for generating a SW search for comparison
import build_expect_scores as bes #<- builds random sequences that are used to run SW and the scores are saved for fitting the distribution
import calc_bit_and_evalues as cbae
import biopython_e as bpe

programme_settings.read()

#Part 1: Building the command line interface

def commandline():
    print("Welcome to the commandline I've been told to make. I hope you are welcomed and find this interface intuitive to use.")
    parser = argparse.ArgumentParser(prog = "BA_ICA_B188292", formatter_class=argparse.RawDescriptionHelpFormatter, description = """\
Welcome to the BLAST101 Search Interface
Analyze your sequence using:
- BLAST101
- Smith-Waterman
- Statistical analysis""")

    #Database Input
    parser.add_argument('--database', required = True, help = "The database sequence")

    #Query Sequence Input
    parser.add_argument('--sequence', required = True, help = "The input sequence")

    #Allowing the user to pick between algorithms and statistical analysis
    parser.add_argument('--algorithm', default = 'blast101', choices=['blast101', 'smith-waterman', 'statistical_analysis'], help = "If you would like to run blast or smith-waterman, select the relevant option. If you would like to conduct statistical analysis, select that option.")

    parser.add_argument('--matrix', default = 'BLOSUM62', help = "BLAST uses the BLOSUM62 table to remove words that are low scoring when perfectly aligned.")
    
    #abbreviating to args
    args = parser.parse_args()

    #Checking we have the input
    if not os.path.exists(args.sequence):
        raise FileNotFoundError(f"Could not find {args.sequence}")

    if not os.path.exists(args.database):
        raise FileNotFoundError(f"Could not identify {args.database}")

    #Checking the file is a .fasta file
    if args.sequence.lower().endswith((".fasta", ".fa")):
        print("Input okay, proceeding...")
        
        #checking the file contains a protein sequence
        prot_alph = set("ACDEFGHIKLMNPQRSTVWY")
        dna_alph = set("ACGT")

        record = next(SeqIO.parse(args.sequence, "fasta"))
        sequence = str(record.seq).upper()

        protein_count = sum(1 for char in args.sequence if char in prot_alph)
        dna_count = sum(1 for char in args.sequence if char in dna_alph)

        if protein_count > dna_count:
            print("This is a protein sequence.\nAccepted.")
        else:
            raise ValueError("Only protein sequences accepted.")

    else:
        raise ValueError(f"Only .fasta files are supported.")

    #Processing the fasta file
    pff.process_fasta_file(fastafile, processfunction, max_scores,ats =None) #<- from the pff library


    #Message to the user
    print(f"Running {args.sequence} against {args.database} using {args.algorithm}")

    #Running the appropriate analysis
    if args.algorithm == "blast101":
       # output_file = f"{args.database}_{args.sequence}_output.out"
        print("Running BLAST101...")
        
        blast.blast101_run() # <- this is from the supplied files
       # print(f"Output has been saved to {output_file}")
        print("Executed BLAST101")

    elif args.algorithm == "smith-waterman":
        print("Running Smith-Waterman...")
        dist = blosum.BLOSUM(int(programme_settings.settings["DEFAULT"]["blosum"])
        dist = blosum.BLOSUM(62)

    elif args.algorithm == "statistical_analysis":
        print("Conducting statistical analysis...")    

    return

commandline()

##############
# References #
##############
#The following handbook was used in creating the commandline:
#https://docs.python.org/3/library/argparse.html
#https://docs.python.org/3/tutorial/datastructures.html#sets

