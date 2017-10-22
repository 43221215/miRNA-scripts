# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:17:44 2017

@author: zigoto
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 17:17:29 2017

@author: zigoto
"""
import re
from collections import defaultdict
import argparse
import numpy as np

# Define parser
parser = argparse.ArgumentParser(description='This script calculates the frequency of a determined NTA respect the total counts of each miRNA (The total number of counts from sequences corresponding to each miRNA without allowing missmatches)')

# Mandatory arguments
mand = parser.add_argument_group("Required arguments")
mand.add_argument('-i', help = "Input file. It mus be the output of the script mirna_classif_and_realcounts.pl from a sam file modified with the function samtools fillmd -e.", type = str, metavar='INFILE', required = True, nargs = '+')
mand.add_argument('-o', help = "Output file name", metavar='OUTFILE', type = str, required = True)
mand.add_argument('-n', help = "Non template addition to analyze. Can be either A-C-G-T. T will be displayed as U in the output file.", metavar='NTA', type = str, required = True)


# Optinal arguments
opt = parser.add_argument_group("Optional arguments")
opt.add_argument('--counts', help = 'Output as total counts corresponding to the selected NTA instead of as percentage reads', action="store_true")

args = parser.parse_args()

files = args.i
outfile = args.o
nta = args.n
output_counts = args.counts

with open(outfile, 'w') as out:
    
    if output_counts == 1:
        out.write('miRNA\tCounts of ' + nta + ' NTA' + '\n') # write outfile header
        
    else:
        out.write('miRNA\tPerctentage of ' + nta + ' NTA' + '\n') # write outfile header
    
    for i in files:
        
        sample = re.sub("-.*", "", i)
        
        with open(i, 'r') as f:
                        
            d = defaultdict() # Create a dict to store the counts of each miRNA and the counts of its NTA. The values of this dict will be a numpy array of 2 values. The first one will represent the number of counts of the canonical miRNA (without NTA) and the second the number of counts for that miRNA with the NTA selected)..
            
            for line in f:
                if line.startswith("#"): continue
                    
                fields = line.split('\t')
                
                counts = fields[1] # Since the file is a collapsed sam, counts must be substracted. I.e: 23525-60, where 60 is the number of counts and 23525 the ID of the read.
                counts = float(re.sub("^.*-", "", counts))                
                mirna = fields[10]
                sequence = fields[6]
                                
                find_nta = '^=+' + re.escape(nta) + '+$' # Build the regex pattern as a string to be able to include the NTA input in the script as a variable inside the regex.
                
                if counts < 5: continue # If the sequence has less than 5 counts skip it.
              
                elif re.match('^=+=$', sequence): # If the sequence is a canonical sequence (= doesn't has missmatches, even if the sequence is an isomiR) sum its counts to its dict entry.
                    
                    if mirna in d: # If the miRNA is already in the hash sum the counts to he canonical miRNA (first position of the list)                       
                        d[mirna] = np.add([counts, 0], d[mirna], dtype='float')
                        
                    else: # If an entry for a mIRNA does not exists create it with a numpy array.           
                        d[mirna] = np.array([counts, 0], dtype='float')
                    
                    
                elif re.match(find_nta, sequence): 
                    
                    if mirna in d: # If the miRNA is already in the hash sum the counts to he canonical miRNA (first position of the list)                       
                        d[mirna] = np.add([0, counts], d[mirna], dtype='float')
                        
                    else: # If an entry for a mIRNA does not exists create it with a numpy array.           
                        d[mirna] = np.array([0, counts], dtype='float')
                    
            for key in d:
                
                if np.sum(d[key]) < 100: continue
                            
                if output_counts == 1:   # If the output must be as total counts corresponding to NTA instead as % of reads...
                    out.write(key + '\t' + str(d[key][1]) + '\n')
                
                else:
                    percentage_NTA = (d[key][1]/(d[key][0]+d[key][1]))*100 # Calculate the % of the NTA that represents respect the number of counts of the canonical miRNA
                    out.write(key + '\t' + str(percentage_NTA) + '\n')
                    
            
# Print the new file in the terminal
with open(outfile, 'r') as f:
    print f.read()