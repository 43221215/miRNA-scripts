#!/usr/bin/perl

# This script will return the number of counts for each miRNA detected (using mirbase database) from a file *.canonical.miRNA.csv. The fasta file from mirbase is used to calculate the number of couns of all the known miRNAs. So, if a miRNA is not detected, it will appear with 0 counts.

# usage --> perl Get_mirnas2DESeq.pl sample.canonical.miRNA.csv mirbase.fasta

use warnings;
use strict;

my $inputfile = $ARGV[0];
my $database = $ARGV[1];
my $outputfile = $inputfile;
$outputfile =~ s/canonical.miRNA.csv/counts.txt/;
my @headers;
my %Counts;  	 # Declare hash where the counts of each miRNA will be sotored.

open(INPUT, "<$inputfile") or die "Couln't open $inputfile";

while (<INPUT>) {
	my $line = $_;
	chomp$line;
	if ($line =~ /^#/) {next}
	@headers = split ("\t", $line);
	if ($headers[1] < 5) {next}   #If the sequence have less than 5 counts, skip it.
	
	if ( ($headers[13] eq "mature" || $headers[13] eq "isomiR") && ($headers[14] !~ /inseed|outseed|firstnt/)){
		$Counts{$headers[10]}+=$headers[1];
	}
}

close INPUT;

my %databaseIDs;
open (DATABASE, "<$database") or die "Couldn't open $database";

while (<DATABASE>) {

	if ($_ =~ /^>.*/){
		my $ID=$_;
		$ID =~ s/ .*$//;
		$ID =~ s/>//;
		chomp$ID;
		$databaseIDs{$ID} = 1;
	}

}

close DATABASE;

open(OUTPUT,">$outputfile") or die "Couldn't open $outputfile";

foreach my $mirna (keys %databaseIDs) {
	if( exists($Counts{$mirna}) ){
		print OUTPUT "$mirna\t$Counts{$mirna}\n";
	}
	else{
		print OUTPUT "$mirna\t0\n";
	}
}

close OUTPUT;

print `head $outputfile`;

exit;
