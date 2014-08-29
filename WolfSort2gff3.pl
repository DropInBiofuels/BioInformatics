#!/usr/bin/perl

# Johannes Kabisch
# kabisch@uni-greifswald.de
# 12.08.2014

# This script converts signalP predictions run with custom table format:
#, based on translated AUGUSTUS output, to a gff3 for jbrowse
# Usage: blastp2gff3.pl

# Example input:
# AUGUSTUS gff3
#NC_006067	AUG5	gene	2659	5277	0.97	+	.	ID=aug5.g1
#NC_006067	AUG5	mRNA	2659	5277	0.97	+	.	ID=aug5.g1.t1;Parent=aug5.g1;Name=aug5.g1.t1;Index=1
#NC_006067	AUG5	start_codon	2659	2661	.	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	CDS	2659	5277	0.97	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	stop_codon	5275	5277	.	+	0	Parent=aug5.g1.t1
#


# Input needs to be converted to tab separated:
# sed -i -e 's/t1\s/t1\t/g' aug5.wolf.out  ## not very good because augustus specicif...looks fot t1\s
#signalP gff
##gff-version 2
##sequence-name	source	feature	start	end	score	N/A ?
## -----------------------------------------------------------
#aug5.g3.t1	WolfPsort-4.1	SIGNAL	1	15	0.667	.	.	YES
#aug5.g4.t1	WolfPsort-4.1	SIGNAL	1	17	0.878	.	.	YES
#aug5.g5.t1	WolfPsort-4.1	SIGNAL	1	15	0.533	.	.	YES
#aug5.g12.t1	WolfPsort-4.1	SIGNAL	1	18	0.765	.	.	YES
#aug5.g19.t1	WolfPsort-4.1	SIGNAL	1	16	0.821	.	.	YES


use warnings;
use strict;

#my $usage = "
#Usage:
#      WolfPsort2gff3 <iprscan.gff> <AUGUSTUS.gff>
#
#      This script takes an InterProScan gff3 format report with labeled
#      IPR domains and builds them into GFF3 formated features that can
#      be added as a track to GBrowse, JBrowse, and Apollo.
#";


#~ if(@ARGV < 2){
    #~ print $usage;
    #~ exit(1);
#~ }

#seqID	subject seqID (accession)	Subject scien name	subj. title (Annotation)	perc. ident.	length	e-val.
	my $seqid;
	my $accession;
	my $start;
	my $stop;
	my $score = 1;
	my $parent;
	my $strand;
	my $phase;
	my $ID;
	my $counter = 0;
	my $localization;
	my $type = "Localization";
	my $source = "WolfPsortv0.2";
	
#parse AUGUSTUS gff3
open(FileHandler_AUGUSTUS, "< aug5.gff3");# or die "ERROR: Could not open the file $TARGET\n";
open(FileHandler_WolfPsort, "< aug5.wolf.out");# or die "ERROR: Could not open the file $IPRSOURCE\n";		
open FileHandler_outFILE, "> TEstWolf_outfile.gff";

while(defined(my $AUGUSTUS_line = <FileHandler_AUGUSTUS>)){
    last if($AUGUSTUS_line =~ /\#\#FASTA/); # Ends as soon as FASTA begins
    next if($AUGUSTUS_line =~ /^\#/); # AUGUSTUS file does not follow strict gff3 guidelines, thus no # is there for next segment
    next unless($AUGUSTUS_line =~ /\tCDS\t/);

    chomp $AUGUSTUS_line;
    my @AUGUSTUS = split("\t", $AUGUSTUS_line);
    if($AUGUSTUS[2] eq 'CDS'){
	$_ = $AUGUSTUS[8];
	/Parent\=/;
	$parent = $';
	$start = $AUGUSTUS[3];
	$stop = $AUGUSTUS[4];
	$seqid = $AUGUSTUS[0];
	$phase = $AUGUSTUS[7];
	$strand = $AUGUSTUS[6];
	#~ print $parent;
#~ Abbrev 	Localization Site 	GO Cellular Component
#~ chlo 	chloroplast 	0009507, 0009543
#~ cyto 	cytosol 	0005829
#~ cysk 	cytoskeleton 	0005856(2)
#~ E.R. 	endoplasmic reticulum 	0005783
#~ extr 	extracellular 	0005576, 0005618
#~ golg 	Golgi apparatus 	0005794(1)
#~ lyso 	lysosome 	0005764
#~ mito 	mitochondria 	0005739
#~ nucl 	nuclear 	0005634
#~ pero 	peroxisome 	0005777(2)
#~ plas 	plasma membrane 	0005886
#~ vacu 	vacuolar membrane 	0005774(2)
	#parse WolfPsort results
		while(defined(my $WolfPsort_line = <FileHandler_WolfPsort>)){ # Filehandler zurücksetzen !!!!! Muss von Anfagn des Files starten
			chomp $WolfPsort_line;
		    my @WolfPsort = split("\t", $WolfPsort_line);
		    #print $WolfPsort[1];
		    if($WolfPsort[0] eq $parent){
				$localization = $WolfPsort[1];
					#~ my $new_localization = $localization =~ /^(.*?),/;
					#~ print $new_localization;
				$ID = $parent.".a".$counter; # For unique ID
				print FileHandler_outFILE "$seqid\t$source\t$type\t$start\t$stop\t$score\t$strand\t$phase\tID=$ID;Name=$localization\n";
			}
		$counter = $counter+1; # For unique ID
		}
		seek FileHandler_WolfPsort,0,0;	
	}
}
close (FileHandler_AUGUSTUS);
close (FileHandler_WolfPsort);
close (FileHandler_outFILE);
