# Description of the scripts

All these scripts were used to perform most of the analysis presented at Fernandez-Perez et al. 2017.

1. **mirnaclassif_and_realcounts.pl**: This script will generate a csv file with the classification of all miRNA sequences respect to the seed region (for a more detailed description, Fernandez-Perez et al. 2017). It uses 2 files as input: a sam file, processed with the samtools command `fillmd -e` (which substitutes all the matched nucleotides by "=", which will be used by the scritps NonTemplateAdditions.py to classify the NTA), and a file that is provided by miRBase called *mmu_miRbase21.str* (you can find it in this repository). To obtain the same data that was represented in the FIG 2A/C in Fernandez-Perez et al. 2017 you will have to process the output file with the script *parse_mirnaclassif_5countsfilt.pl*. 

2. **parse_mirnaclassif_5countsfilt.pl**: Parses the output of *mirnaclassif_and_realcounts.pl* to obtain a table with the total number of reads of the sequences classified as NoChange, NonCanonicalProcessing, Inseed, Outseed and Firstnt.

3. **mmu_miRbase21.str**: This file is required by the script *mirnaclassif_and_realcounts.pl*.

4. **Get_mirnas_counts_NoChange.pl**: This script will parse the output of *mirnaclassif_and_realcounts.pl* to generate a miRNA expression table. It will just consider those sequences classified as NoChange.

5. **DEseq.R**: The R script that was used at Fernandez-Perez et al. 2017 to perform the differential expression analyses. It uses as input a table generated with the output of *Get_mirnas_counts_NoChange.pl*.

6. **calculate_canonical-isomiR_expression.pl**: This script uses the output of *mirnaclassif_and_realcounts.pl* to create a table containing the % of reads corresponding to isomiR for each miRNA. 

7. **NonTemplateAdditions_global.py**: Script to calculate the global % of different NTA. The % is respect the total number of NTA reads, not total reads.

8. **NonTemplateAdditions_by_miRNA.py**: Script to calculate the % of reads that represent a determined NTA for each individual miRNA. Has an option to give the number of counts from that NTA instead of the % of reads.


For any question, contact me at daniel.fernandezperez@ieo.it
