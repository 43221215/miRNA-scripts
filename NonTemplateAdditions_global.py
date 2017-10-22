# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 17:17:29 2017

@author: zigoto
"""
import re
from collections import defaultdict
import argparse

# Define parser
parser = argparse.ArgumentParser(description='This script calculates the frequencies of the NTAs for 1 or more files. The proportion that represents each type of NTA is calculated with the total number of counts of all the NTA')

# Mandatory arguments
mand = parser.add_argument_group("Required arguments")
mand.add_argument('-i', help = "Input file. It mus be the output of the script mirna_classif_and_realcounts.pl from a sam file modified with the function samtools fillmd -e. The filename must be name-.whatever. All the text before the '-' will be the header of the output table.", 
                  type = str, metavar='INFILE', required = True, nargs = '+')
                  
mand.add_argument('-o', help = "Output file name", metavar='OUTFILE', type = str, required = True)

# Optinal arguments
opt = parser.add_argument_group("Optional arguments")
opt.add_argument('--mir', help = "Filter by a specific miRNA to calculate its NTAs (optional)", metavar='miRNA', type = str)

opt.add_argument('--canonical', help = "Add one more variable to the table of NTAs that corresponds to the expression of canonical sequences (without NTA, isomiRs included) (optional)", 
                 action="store_true")


args = parser.parse_args()

files = args.i
outfile = args.o
filt_mirna = args.mir

if args.canonical: # Set canonical output to true or false
    canonical = 1
else: canonical = 0

with open(outfile, 'w') as out:
    
    if canonical == 1:
                out.write('Sample\tCanonical\tA\tU\tC\tG\toligo_A\toligo_U\toligo_C\toligoG\n') # write outfile header
    
    else:  
        out.write('Sample\tA\tU\tC\tG\toligo_A\toligo_U\toligo_C\toligoG\n') # write outfile header
         
    for i in files:
        
        sample = re.sub("-.*", "", i)
        
        with open(i, 'r') as f:
            
            d = defaultdict(int) # Create a dict to store the counts of each miRNA without allowing missmatches
            d_NTA = defaultdict(int)
            for line in f:
                if line.startswith("#"): continue # Skip the header of the file produced by the script mirna_classif_and_realcounts.pl
                    
                fields = line.split('\t')
                counts = fields[1] # Since the file is a collapsed sam, counts must be substracted. I.e: 23525-60, where 60 is the number of counts and 23525 the ID of the read.
                counts = float(re.sub("^.*-", "", counts))
                
                mirna = fields[10]
                sequence = fields[6]

                
                if (mirna != filt_mirna and filt_mirna != None) or counts < 5: continue # If the sequence has less than 5 counts or is not the miRNA set at --mir, skip it.
            
              
                elif canonical == 1 and re.match('^=+=$', sequence): # If the sequence is a canonical sequence (= doesn't has missmatches, even if the sequence is an isomiR) sum its counts to its dict entry
                    d_NTA['canonical'] += counts
                    
                    
                elif re.match('^=+A$', sequence): # Just 1 A as NTA, monoA
                    d_NTA['A'] += counts
                
                elif re.match('^=+T$', sequence): # Just 1 T as NTA, monoU
                    d_NTA['U'] += counts
                    
                elif re.match('^=+C$', sequence): # Just 1 C as NTA, monoC
                    d_NTA['C'] += counts   
                
                elif re.match('^=+G$', sequence): # Just 1 G as NTA, monoG
                    d_NTA['G'] += counts
                    
                
                elif re.match('^=+AA+$', sequence): # 2 A or more as NTA, oligoA
                    d_NTA['oligo_A'] += counts
                
                elif re.match('^=+TT+$', sequence): # 2 T or more as NTA, oligoU
                    d_NTA['oligo_U'] += counts
                    
                elif re.match('^=+CC+$', sequence): # 2 T or more as NTA, oligoC
                    d_NTA['oligo_C'] += counts   
                
                elif re.match('^=+GG+$', sequence): # 2 T or more as NTA, oligoG
                    d_NTA['oligo_G'] += counts
                        
                
            total_NTA = sum(d_NTA.itervalues()) # Sum all the counts stored in the dictionary to calculate the % of these counts that represents each NTA.
     
    # Print to output

            if total_NTA == 0: # if there isn't a single read that matches the previous criteria in the file... return all the NTA with 0.
                out.write(sample + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\n')

            elif canonical == 1:
                out.write(sample + '\t' + str(d_NTA['canonical']/total_NTA) + '\t' + str(d_NTA['A']/total_NTA) + '\t' + str(d_NTA['U']/total_NTA) + '\t' + str(d_NTA['C']/total_NTA) + '\t' + str(d_NTA['G']/total_NTA) + '\t' + str(d_NTA['oligo_A']/total_NTA) + '\t' + str(d_NTA['oligo_U']/total_NTA) + '\t' + str(d_NTA['oligo_C']/total_NTA) + '\t' + str(d_NTA['oligo_G']/total_NTA) + '\n')
            
            else:  
                out.write(sample + '\t' + str(d_NTA['A']/total_NTA) + '\t' + str(d_NTA['U']/total_NTA) + '\t' + str(d_NTA['C']/total_NTA) + '\t' + str(d_NTA['G']/total_NTA) + '\t' + str(d_NTA['oligo_A']/total_NTA) + '\t' + str(d_NTA['oligo_U']/total_NTA) + '\t' + str(d_NTA['oligo_C']/total_NTA) + '\t' + str(d_NTA['oligo_G']/total_NTA) + '\n')
        
# Print the new file in the terminal
with open(outfile, 'r') as f:
    print f.read()
    