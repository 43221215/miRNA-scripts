#!/usr/bin/perl

# This script will classify the output of mirnasclassif_and_realcounts.pl with the following criteria: Mature and isomirs with 0 missmatches will be classified as "No change".
# Then, mature and isomirs with missmatches are classified as "Inseed" (at least 1 of the missmatches is in the seed) and "outseed" (all the mismatches outseed).
# Finally, isoforms with trimmings or elongations at 5' end (=different seed and different targets) will be classified as "Non Canonical Processing".

# usage --> parse_mirnaclassif_5countsfilt.pl file.canonical.miRNA.csv. The output file will be "file.miRNA.classification.csv".
# This version ignore all sequences with less than 5 counts.


use warnings;
use strict;


open(OUTPUT, ">isomir_classification.txt") or die;
print OUTPUT "Sample\tOutseed\tInseed\tFirst Nucleotide\tNo Change\tNon Canonical Processing\n";

foreach(@ARGV) {

my $mirna_classif = $_;
my $sample = $mirna_classif;
$sample =~ s/.R2T.fq.fasta.collapsed.mapPREmir21.v3.canonical.miRNA.csv$//;

my $TotalCounts = 0; #Initiate a variable to calculate the total number of counts of the sample.
my @headers;
my %Classification;

open(CLASIF, "<$mirna_classif") or die "Couldn't open mirna classification file";

# Headers -->  Sequence ID	Real counts	miRbase	Algmt Start	Algmt End	Seq Length	Sequence	CIGAR	No. mismatches	Canonicals	Canonical Match	Canon Start	Canon End	miRNA type	Variant	Substitution

while (<CLASIF>) {
	my $line = $_;
	if ($line =~ /^#/) {next}
	@headers = split ("\t", $line);
	if ($headers[1] < 5) {next}   #If the sequence have less than 5 counts, skip it.
	$TotalCounts += $headers[1];  #Track the total number of counts


	if ($headers[8] == 0 && ($headers[13] eq "mature" || $headers[13] eq "isomiR")) {     # If the sequence has 0 mismatches and it's classified as a mature, is the canonical miRNA without any variation (the one in mirbase).
		$Classification{"NoChange"} += $headers[1];	# Sum the number of counts of that miRNA. 
	}

	if ($headers[8] != 0 && ($headers[13] eq "mature" || $headers[13] eq "isomiR")) {    # Matures that have modifications in their sequence, SNP or ADAR modification

		if ($headers[14] =~ /firstnt/) {			# 1st nt change. If variation is first nt + inseed or outseed will be classified as first nt.

			$Classification{"FirstNucleotide"} += $headers[1];
		}		
		
		elsif ($headers[14] =~ /inseed/) { 			# If 1 of the mutations is inseed, even if there are 2 more outseed, it will be classified as inseed.
			$Classification{"Inseed"} += $headers[1];
		}

		elsif ($headers[14] =~ /outseed/){							#If is not inseed, firstnt or no change, it has to be outseed.
			$Classification{"Outseed"} += $headers[1];
		}
	 }

		
	if ($headers[13] eq "circunmiR" || $headers[13] eq "paramiR") {
		$Classification{"NonCanonicalProcessing"} += $headers[1];
		if ($headers[13] eq "circunmiR") { $Classification{"circunmiR"} += $headers[1] }
		if ($headers[13] eq "paramiR") { $Classification{"paramiR"} += $headers[1] }
	}

}
close CLASIF;

# Write results

print OUTPUT "$sample\t$Classification{Outseed}\t$Classification{Inseed}\t$Classification{FirstNucleotide}\t$Classification{NoChange}\t$Classification{NonCanonicalProcessing}\n";

}


close OUTPUT;

print `head isomir_classification.txt`;
exit

