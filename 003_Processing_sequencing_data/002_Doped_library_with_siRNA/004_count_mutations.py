#!/usr/bin/python
# Pablo Baeza


# import libraries

import sys
from operator import ne as not_equal
from Bio import SeqIO
import time


# want to use this script from the command line with 2 arguments: input file and original DNA sequence

if len(sys.argv) != 3:
	sys.exit("ERROR - WRONG NUMBER OF ARGUMENTS.\nUsage: count_mutations.py FILE REFERENCE_SEQUENCE")

# merged and collapsed file
input_file_name = sys.argv[1]

# original DNA sequence with which all other occurrences will be compared
ref_sequence = sys.argv[2]



# define functions

def hamming(seq1, seq2):
	"""
	Count number of nucleotide positions at which 2 input sequences differ
	"""
	assert len(seq1) == len(seq2) #seqs need to be of equal length
    
    # removing for loop because it is less efficient than map()
    # mutations = 0
    # for i in range(len(seq1)):
    #    if seq1[i] != seq2[i]:
    #        mutations += 1
    # return mutations
    
    # map iterates through seq1 and seq2 at the same time, applying
    # not_equal to each character; returns a list with each result
	is_it_different = map(not_equal, seq1, seq2)
    
    # result in boolean format, so can just sum
	no_of_mutations = sum(is_it_different)
    
	return no_of_mutations
    
def timestamp():
    """
    Return local time
    """
    return time.strftime("[%a, %d %b %Y, %H:%M:%S] ", time.localtime())



# output start time
print >> sys.stderr, "**************************"
print >> sys.stderr, timestamp() + "Starting"
print >> sys.stderr, "Input file=" + str(input_file_name)
print >> sys.stderr, "Reference sequence=" + str(ref_sequence)
print >> sys.stderr, "**************************"



# MAIN BIT - count mutations and output as a table containing variant sequence, no of occurrences
# and no of mutations relative to the reference sequence

#input_file_handle = open(input_file_name, "r")
parsed_file_iterator = SeqIO.parse(input_file_name, "fasta") #want to parse as a FastA file
#input_file_handle.close()

print("Genotype\tOccurrences\tMutationCount")
for each_record in parsed_file_iterator:
	recordID = each_record.id.split("-")[1]
	print(each_record.seq + "\t" + recordID + "\t" + str(hamming(ref_sequence, each_record.seq)))
	
	
	
# output finish time
print >> sys.stderr, "..."	
print >> sys.stderr, timestamp() + "Done!"