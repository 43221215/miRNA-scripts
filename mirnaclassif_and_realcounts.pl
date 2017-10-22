#!/usr/bin/perl 
use strict; 
use warnings;

$|=1;


# Usage: mirnasclassif_and_real_counts.pl oocyte_aligned3_miRNA20.sam mmu_miRNA.str


my $parsedfile = $ARGV[0];
my $miranalyzerfile = $ARGV[1];

my $canonicfile = $parsedfile;
$canonicfile =~ s/sam$/canonical.miRNA.csv/;
`touch $canonicfile`;
my $precursorfile = $parsedfile;
$precursorfile =~ s/sam$/precursor.miRNA.csv/;
`touch $precursorfile`;

open (PARSEDFILE, "<$parsedfile") or die "Could not open SAMFILE $parsedfile";
open (CANONFILE, ">$canonicfile") or die "Could not open CANONFILE $canonicfile";
open (PRECFILE, ">$precursorfile") or die "Could not open PRECFILE $precursorfile";


print CANONFILE "#Sequence ID\tReal counts\tmiRbase\tAlgmt Start\tAlgmt End\tSeq Length\tSequence\tCIGAR\tNo. mismatches\tCanonicals\tCanonical Match\tCanon Start\tCanon End\tmiRNA type\tVariant\tSubstitution\n";
print PRECFILE "#Sequence ID\tReal counts\tmiRbase\tAlgmt Start\tAlgmt End\tSeq Length\tSequence\tMatch?\tNo. mismatches\tCanonicals\tmiRNA type\n";

my $linea = '';

my ($sequenceID, $miranalyzer, $bitwiseflag, $alignstart, $length, $sequence, $rawcounts, $realcounts, $cigarmismatches, $nmismatches, $alignend, $matchinmir, $campomir, $canonico, $numcanon, $rango, $canonstart, $canonend, $canonical2check, $colon, $canonstart2check, $canonend2check, $mirnatype, $variant, $variant_cumul, $substitution_ref, $substitution_alt, $substitution, $substitution_cumul, $cigarstart)= '';
my $matchingstatus = 0;
my ($i, $j)=2;
my $n = -1;
my ($substit_pos, $stepped_positions) = 0;
my (@fields, @camposmir, @canonicos, @canonrange, @substits) = ();

my $nmatches = 1;


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
my %Num_matches;

while (<PARSEDFILE>) {
  $linea = $_;
  @fields = split ("\t",$_);
  $bitwiseflag = $fields[1];
    if ($linea =~ /^@/) { next;}
    elsif ( $bitwiseflag == '4' ) { next;}
	else {
    chomp $linea;
    $sequenceID = $fields[0];
    $Num_matches{$sequenceID} += 1;
	}
}

my %Miranalyzer_file;

open (MIRANALYZERFILE, "<$miranalyzerfile") or die "Could not open MIRANALYZERFILE $miranalyzerfile";

while (<MIRANALYZERFILE>){
	my $line = $_;
	chomp $line;
	@camposmir = split (" ",$line);
	my $mirID = $camposmir[0];
	$mirID =~ s/>//;
	$Miranalyzer_file{$mirID} = $line;
}

close MIRANALYZERFILE;


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open (PARSEDFILE, "<$parsedfile") or die "Could not open SAMFILE $parsedfile";
    
while (<PARSEDFILE>) {
  $linea = $_;
  @fields = split ("\t",$_);
  $bitwiseflag = $fields[1];
  if ($linea =~ /^@/) { next;}
  elsif ( $bitwiseflag == '4' ) { next;}
  else {
	chomp $linea;
	$sequenceID = $fields[0];
	$miranalyzer = $fields[2];
	$alignstart = $fields[3];
	$length = $fields[5];
	$length =~ s/M$//;
	$sequence = $fields[9];
	$cigarmismatches = $fields[12];
	$cigarmismatches =~ s/^MD:Z://;
	$nmismatches = $fields[13];
	$nmismatches =~ s/^NM:i://;
 	$nmismatches = int ($nmismatches);
	$alignend = $alignstart + $length -1;
	
	$nmatches = $Num_matches{$sequenceID};
	
	$rawcounts = $sequenceID;	
	$rawcounts =~ s/^.*-//;
	$rawcounts = int ($rawcounts);
	$nmatches = int ($nmatches);
	$realcounts = $rawcounts/$nmatches;	
	
	$matchinmir = $Miranalyzer_file{$miranalyzer};
	@camposmir = split (" ",$matchinmir);
	$matchingstatus = 0;
	
		for ($i=2; $i<=$#camposmir; $i++) {		
	  @canonicos = split (":", $camposmir[$i]);
	  $rango = $canonicos[1];
	  @canonrange = split ("-", $rango);
	  $canonstart = $canonrange[0];
	  $canonstart = int ($canonstart);
	  $canonend = $canonrange[1];
	  $canonend =~ s/]$//;
	  $canonend = int ($canonend);

	  if ($matchingstatus == 1) { next;};
	  elsif ( ($alignstart>=$canonstart && $alignstart<=$canonend) ||	
	          ($alignend>=$canonstart && $alignend<=$canonend) ||	
	          ($alignstart<$canonstart && $alignend>$canonend) ) {
	      $matchingstatus = 1;
	      $canonical2check = $camposmir[$i];
	      $colon = index ($canonical2check, ':');
	      $canonical2check = substr ($canonical2check, 1, $colon-1);
	      $canonstart2check = $canonstart;	
	      $canonend2check = $canonend,	
	      next;			
	  }
	  else {		
	      $matchingstatus = 2; 

	  }		
	}  

 
	  if ($matchingstatus == 1) {
	      $canonico = '';
	      for ($j=2; $j<=$#camposmir; $j++) {
		  $canonico = $canonico.' '.$camposmir[$j];
	      }
	      $substitution_cumul = '';
	      $substitution = '';
	      $variant = '';
	      $variant_cumul = '';
	      $stepped_positions = 0; 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ mature ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      if ($alignstart==$canonstart2check && $alignend==$canonend2check){ 
		  $mirnatype = "mature";	
		  if ($nmismatches == 0) {
		      $variant = "nochange";
		      $variant_cumul = $variant;
		      $substitution = "nochange";
		      $substitution_cumul = $substitution;
		  }
		  elsif ($nmismatches != 0) {
		      @substits = split (/(\d+)/, $cigarmismatches);
		      for ($n=1; $n<($#substits); $n=$n+2) {
			  $substit_pos = $substits[$n];			
			  $substitution_ref = $substits[$n+1];	
			  $stepped_positions = $stepped_positions + $substit_pos;
			  
			  if ($stepped_positions == 0) {
			  $variant = "firstnt ";	
			  $variant_cumul = $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1); 
			  $substitution = 'f:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;
			  $stepped_positions++;
			  next;
			  }
			  elsif ($stepped_positions > 0 && $stepped_positions <=7) {
			  $variant = "inseed ";
			  $variant_cumul .= $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1);
			  $substitution = 'i:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;
			  $stepped_positions++;
			  next;
			  }
			  elsif ($stepped_positions >7) {
			  $variant = "outseed ";
			  $variant_cumul .= $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1);
			  $substitution = 'o:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;	
			  $stepped_positions++;	
			  next;
			  }
			  else {
			  $substitution_cumul .= $substitution;	
			  $variant_cumul .= $variant;
			  next; }
		      }
		  }   
	      }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ isomiR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      elsif ($alignstart==$canonstart2check) {
		  $mirnatype = "isomiR";	
		  if ($nmismatches == 0) {
		  $variant = "nochange";
		  $variant_cumul = $variant;
		  $substitution = "nochange";
		  $substitution_cumul = $substitution;
		  }
		  elsif ($nmismatches != 0) {
		      @substits = split (/(\d+)/, $cigarmismatches); 
		      for ($n=1; $n<($#substits); $n=$n+2) {
			  $substit_pos = $substits[$n];			
			  $substitution_ref = $substits[$n+1];		
			  $stepped_positions = $stepped_positions + $substit_pos;	
			  
			  if ($stepped_positions == 0) {
			  $variant = "firstnt ";		
			  $variant_cumul = $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1); 
			  $substitution = 'f:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;
			  $stepped_positions++;
			  next;
			  }
			  elsif ($stepped_positions > 0 && $stepped_positions <=7) {
			  $variant = "inseed ";
			  $variant_cumul .= $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1);
			  $substitution = 'i:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;
			  $stepped_positions++;
			  next;
			  }
			  elsif ($stepped_positions >7) {
			  $variant = "outseed ";
			  $variant_cumul .= $variant;
			  $substitution_alt = substr ($sequence, $stepped_positions, 1);
			  $substitution = 'o:'.$substitution_ref.'>'.$substitution_alt."\t";
			  $substitution_cumul .= $substitution;	
			  $stepped_positions++;	
			  next;
			  }
			  else {
			  $substitution_cumul .= $substitution;	
			  $variant_cumul .= $variant;
			  next; }
		      }
		  }
	      }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ paramiR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      elsif ($alignstart<$canonstart2check) {
		  $mirnatype = "paramiR";	
		  $variant = "--";
		  $variant_cumul = $variant;
		  $substitution = "--";
		  $substitution_cumul = '';
	      }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ circunmiR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      else {
		  $mirnatype = "circunmiR";
		  $variant = "--";
		  $variant_cumul = $variant;
		  $substitution = "--";
		  $substitution_cumul = '';
	      }
	      print CANONFILE "$sequenceID\t$realcounts\t$miranalyzer\t$alignstart\t$alignend\t$length\t$sequence\t$cigarmismatches\t$nmismatches\t$canonico\t$canonical2check\t$canonstart2check\t$canonend2check\t$mirnatype\t$variant_cumul\t$substitution_cumul\n";
	  }

	  elsif ($matchingstatus == 2) {
	      $mirnatype = "precursor";
	      $canonico = '';
	      for ($j=2; $j<=$#camposmir; $j++) {
	      $canonico = $canonico.' '.$camposmir[$j];
	      }
	      print PRECFILE "$sequenceID\t$realcounts\t$miranalyzer\t$alignstart\t$alignend\t$length\t$sequence\t$cigarmismatches\t$nmismatches\t$canonico\t$mirnatype\n";
	  }
	  else {print "sequence NOT MATCHING!!!";}
  }
}


close (PARSEDFILE);
print "\n";
print `head -20 $precursorfile`;
print "\n";

close (CANONFILE);
close (PRECFILE);
