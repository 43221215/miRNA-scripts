#!/usr/bin/perl


# usage: calculate_canonical-isomiR_expression.pl file1.csv file2.csv... --> These files are the output fo the script mirnacassif_and_real_counts.pl. You can use as many files as you want as input. This will generate a table where with the % of reads that are isomirs for each miRNA (just perfect matches) 
# This script assumes that the files have a header, it skips the first line+
# The name of the files should be: sample1-whateveryouwant.csv. Everything that is before the '-' will be used to generate the output header, with the names of the samples. It is important the '-', the script will look for it.
# You will need also R installed and the perl R statistics installed

use warnings;
use strict;
use Statistics::R;

my $numberOfFiles = @ARGV; 
my @files = @ARGV; 
for (@files) {
	s/-.*$//;
}

for (@ARGV) {	
	read_file($_);
}

my %matrix_hash;

my @cell = @files;

open (ARRAYFILE, ">isomirs-canonical_table.txt");
print ARRAYFILE "miRNA\t", join("\t", @files), "\n";


foreach my $seqID (keys %matrix_hash) {
  my $counts2print = '';


  for (my $i=0 ; $i<=($numberOfFiles-1); $i++) {
    my $counts_in_cell = $matrix_hash{$seqID};
    
    if ( $matrix_hash{$seqID} =~/\b$cell[$i]-/i ) {   
	 $counts_in_cell =~ s/^.*$cell[$i]-//i; # 859665    zyg-413589      pgc-6558469     spg-4469621
	 $counts_in_cell =~ s/\t.*$//; # 859665
    }
    else { 
	 $counts_in_cell = "NA";

    }
    $counts2print .= $counts_in_cell."\t";

  }
	$counts2print =~ s/\t$//;
    print ARRAYFILE "$seqID\t$counts2print\n"; 
}

close (ARRAYFILE); 


#~~~~~~~~~~~~~Use R to filter rows that have all NA values~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $R = Statistics::R->new(); #Create a communication bridge with R and start R

#Now filter all rows that have NA values in every sample
$R->run( q`table <- read.delim("isomirs-canonical_table.txt"); 

	filt <- apply(table[ , 2:ncol(table)], 1, function(x) { any(!is.na(x))}); 

	table.ok <- table[filt , ];

	write.table(table.ok, file="isomirs-canonical_table.txt", sep="\t", quote=FALSE, row.names = FALSE)` );

$R->stop(); # Stop R

print `head isomirs-canonical_table.txt`;

exit;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub read_file {

  my %mature_hash;
  my %isomir_hash;

  my $miRNAclasif = shift @_;
  my $fileshortname = $miRNAclasif;
  $fileshortname =~ s/(-).*$/$1/;

  open (CLASIF, "<$miRNAclasif") or die "Could not open $miRNAclasif";
  my $header = <CLASIF>;
  while (<CLASIF>) {

	my $line = $_; 
	chomp $line;

	my @headers = split ("\t", $line);
	if ($headers[1] < 5) {next}   #If the sequence have less than 5 counts, skip it.

	if ( ($headers[13] eq "mature") && ($headers[8] <= 1) ){ #Count of mature miRNAs without missmatches
		$mature_hash{$headers[10]} += $headers[1];
	}

	if ( ($headers[13] ne "mature") && ($headers[8] <= 1) ){ #Counts of isomiRs, paramiRs and circunmiRs without missmatches
		$isomir_hash{$headers[10]} += $headers[1];
	}

  }


  foreach my $key (keys %isomir_hash) { #Now it's necessary to loop over the 2 hashes to be sure that both hashes have the same keys and there isn't one hash that has more than the other (for example, if for a specific miRNA just counts from the canonical isoform have been detected).
	unless (exists $mature_hash{$key}) {
		$mature_hash{$key} = 0;
	}
  }

  foreach my $key (keys %mature_hash) {
	unless (exists $isomir_hash{$key}) {
		$isomir_hash{$key} = 0;
	}
  }

  foreach my $miRNA (keys %isomir_hash) { #Now loop over all the keys in the hash to calculate the relation isomirs/canonical miRNA

	if ( $mature_hash{$miRNA} >= 50 || $isomir_hash{$miRNA} >= 50 ) {

		my $canonical_isomir = $isomir_hash{$miRNA}/($mature_hash{$miRNA}+$isomir_hash{$miRNA})*100;

		$matrix_hash{$miRNA} .= $fileshortname.$canonical_isomir."\t"; # Add to the main hash the proportion of canonical/ismomirs for each miRNA of the sample --> sample-ratio
	}
	else { # If there isn't enough counts to calculate the relation isomirs-canonicals (at least 100 counts in total, isomirs+canonicals) value = NA.
		$matrix_hash{$miRNA} .= $fileshortname."NA"."\t";
	}
  }
  close (CLASIF);
}
