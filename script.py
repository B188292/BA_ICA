#!/bin/python3
from Bio import SeqIO
import argparse, subprocess, os, datetime
import smith_waterman_p as sw_p
import blast_101_search as blast #<- runs the main search for blast
import programme_settings
import process_fasta_file as pff #<- used for generating a SW search for comparison
import build_expect_scores as bes #<- builds random sequences that are used to run SW and the scores are saved for fitting the distribution
import calc_bit_and_evalues as cbae
import biopython_e as bpe

programme_settings.read()

##############################################
#Part 1: Building the command line interface #
##############################################

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
        
        blast_101_output = f"{args.sequence}_{args.database}_blast_output.txt"
        with open(blast_101_output, 'a') as bo:
            bo.write("BLAST101 Output:")
            blast.extend_diagonal(pos_s0_s1,s0,s1)
            blast.process_blast(myline_database)
            blast.process_fasta_file()
            blast.blast101_run()
            bo.write(blast.print_final_results(res))
            bo.write(blast.print_timer())
            bo.write(f"BLAST Completed at: {time.strftime('%Y-%m-%d %H:%M')}")
)
       # print(f"Output has been saved to {output_file}")
        print("Executed BLAST101")

    elif args.algorithm == "smith-waterman":
        import blosum

        print("Running Smith-Waterman...")
        sw_output = f"{args.sequence}_{args.database}_sw.txt"

        with open(sw_output, 'a') as swo:
            #The following functions are from the smith_waterman_p script
            swo.write("Smith-Waterman Output:")
            sw_p.create_matrix(rows, cols)
            swo.write(sw_p.calc_score(matrix, x, y))
            sw_p.traceback(mymatrix, maxv)
            sw_p.build_matrix(mymatrix)
            swo.write(sw_p.get_max(mymatrix))
            swo.write(sw_p.print_matrix(mymatrix, logger = None))
            swo.write(sw_p.print_traceback(mymatrix))
            swo.write(sw_p.perform_smith_waterman(seq1,seq2,print_m=False,print_a = False))
            swo.write(sw_p.test())
            swo.write(f"Smith-Waterman Completed at: {time.strftime('%Y-%m-%d %H:%M')}")


    elif args.algorithm == "statistical_analysis":
        print("Conducting statistical analysis...")    
        #Using functions from calc_bit_and_evalues
        stat_output = f"{args.sequence}_{args.database}_stats.txt"

        with open(stat_output, 'a') as so:
            so.write("The results from your statistical analysis")
            so.write(cbae.build_fit())
            so.write(cbae.get_bit_score(raw_score, k, scale))
            so.write(cbae.get_p_value(raw_score, k, scale))
            so.write(cbae.get_expect(raw_score,k, scale))
            so.write(cbae.get_bit_score_s(raw_score, k, scale, precision =0))
            so.write(cbae.get_expect_s(raw_score,k, scale,precision = 5))
            so.write(cbae.test())
            so.write(f"Analysis completed at: {time.strftime('%Y-%m-%d %H:%M')}")

    return

commandline()

############################
# Part 2: Testing the Code #
############################

##############
# References #
##############
#The following handbook was used in creating the commandline:
#https://docs.python.org/3/library/argparse.html
#https://docs.python.org/3/tutorial/datastructures.html#sets

