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

#signalP gff
##gff-version 2
##sequence-name	source	feature	start	end	score	N/A ?
## -----------------------------------------------------------
#aug5.g3.t1	SignalP-4.1	SIGNAL	1	15	0.667	.	.	YES
#aug5.g4.t1	SignalP-4.1	SIGNAL	1	17	0.878	.	.	YES
#aug5.g5.t1	SignalP-4.1	SIGNAL	1	15	0.533	.	.	YES
#aug5.g12.t1	SignalP-4.1	SIGNAL	1	18	0.765	.	.	YES
#aug5.g19.t1	SignalP-4.1	SIGNAL	1	16	0.821	.	.	YES


use warnings;
use strict;

#my $usage = "
#Usage:
#      SignalP2gff3 <iprscan.gff> <AUGUSTUS.gff>
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
	my $score;
	my $parent;
	my $strand;
	my $phase;
	my $ID;
	my $counter = 0;
	my $AAstop;
	my $type = "Signal Peptide";
	my $source = "SignalP-4.1";
	
#parse AUGUSTUS gff3
open(FileHandler_AUGUSTUS, "< aug5.gff3");# or die "ERROR: Could not open the file $TARGET\n";
open(FileHandler_SignalP, "< aug5.sigP.gff");# or die "ERROR: Could not open the file $IPRSOURCE\n";		
open FileHandler_outFILE, "> sigP_outfile.gff";

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

	#parse SignalP results
		while(defined(my $SignalP_line = <FileHandler_SignalP>)){ # Filehandler zurücksetzen !!!!! Muss von Anfagn des Files starten
			#next unless($SignalP_line !~ /\t$organism\t/); # skip duplicate annotations, defined by organism
		    chomp $SignalP_line;
		    my @SignalP = split("\t", $SignalP_line);
		    #print $SignalP[1];
		    if($SignalP[0] eq $parent){
				$score = $SignalP[5];
				$AAstop = $SignalP[4];
				my $stopAAcalc = ($AAstop*3)+$start;
				$ID = $parent.".a".$counter; # For unique ID
				print FileHandler_outFILE "$seqid\t$source\t$type\t$start\t$stopAAcalc\t$score\t$strand\t$phase\tID=$ID\n";
			}
		$counter = $counter+1; # For unique ID
		}
		seek FileHandler_SignalP,0,0;	
	}
}
close (FileHandler_AUGUSTUS);
close (FileHandler_SignalP);
close (FileHandler_outFILE);
