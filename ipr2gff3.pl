#!/usr/bin/perl

# Johannes Kabisch
# kabisch@uni-greifswald.de
# 12.08.2014

# This script converts InterProScan predictions, based on translated AUGUSTUS output, to a gff3 for jbrowse
# Usage: ipr2gff3.pl < output.ipr.gff > cegma.gff3

# Example input:
# AUGUSTUS gff3
#NC_006067	AUG5	gene	2659	5277	0.97	+	.	ID=aug5.g1
#NC_006067	AUG5	mRNA	2659	5277	0.97	+	.	ID=aug5.g1.t1;Parent=aug5.g1;Name=aug5.g1.t1;Index=1
#NC_006067	AUG5	start_codon	2659	2661	.	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	CDS	2659	5277	0.97	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	stop_codon	5275	5277	.	+	0	Parent=aug5.g1.t1
#
# InterProScan gff3
#aug5.g1.t1_1	.	polypeptide	1	540	.	+	.	md5=46307eef4fee927632f8d6e8eeb04eba;ID=aug5.g1.t1_1
#aug5.g1.t1_1	Pfam	protein_match	151	540	1.8E-97	+	.	Name=PF03169;signature_desc=OPT oligopeptide transporter protein;Target=aug5.g1.t1_1 151 540;status=T;ID=match$1_151_540;Ontology_term="GO:0055085";date=17-07-2014;Dbxref="InterPro:IPR004813"
#aug5.g1.t1_1	TIGRFAM	protein_match	134	540	8.2E-63	+	.	Name=TIGR00728;signature_desc=OPT_sfam: oligopeptide transporter, OPT superfamily;Target=aug5.g1.t1_1 134 540;status=T;ID=match$2_134_540;Ontology_term="GO:0055085";date=17-07-2014;Dbxref="InterPro:IPR004813"


# Example output:
#NC_006067	AUG5	gene	2659	5277	0.97	+	.	ID=aug5.g1
#NC_006067	IPR-Scan	pfam	3112	4279	1.8E-97	+	.	Name=PF03169;signature_desc=OPT oligopeptide transporter protein;Target=aug5.g1.t1_1 151 540;status=T;ID=match$1_151_540;Ontology_term="GO:0055085";date=17-07-2014;Dbxref="InterPro:IPR004813"
#NC_006067	IPR-Scan	TIGRFAM	3061	4279	8.2E-63	+	.	Name=TIGR00728;signature_desc=OPT_sfam: oligopeptide transporter, OPT superfamily;Target=aug5.g1.t1_1 134 540;status=T;ID=match$2_134_540;Ontology_term="GO:0055085";date=17-07-2014;Dbxref="InterPro:IPR004813"
#NC_006067	AUG5	mRNA	2659	5277	0.97	+	.	ID=aug5.g1.t1;Parent=aug5.g1;Name=aug5.g1.t1;Index=1
#NC_006067	AUG5	start_codon	2659	2661	.	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	CDS	2659	5277	0.97	+	0	Parent=aug5.g1.t1
#NC_006067	AUG5	stop_codon	5275	5277	.	+	0	Parent=aug5.g1.t1

use warnings;
use strict;

#my $usage = "
#Usage:
#      ipr2gff3 <iprscan.gff> <AUGUSTUS.gff>
#
#      This script takes an InterProScan gff3 format report with labeled
#      IPR domains and builds them into GFF3 formated features that can
#      be added as a track to GBrowse, JBrowse, and Apollo.
#";


#~ if(@ARGV < 2){
    #~ print $usage;
    #~ exit(1);
#~ }
	my $type;
	my $AAstart;
	my $AAstop;
	my $start;
	my $stop;
	my $score;
	my $annotation;
	my $parent;
	my $seqid;
	my $strand;
	my $phase;
	my @AUGUSTUS;
	my @IPR;


#my $IPRSOURCE = "aug5.gff3";#shift; # $file
#my $TARGET = "testset.pep.gff3";#shift;		# $gff3

#die("ERROR: The file $IPRSOURCE does not exist\n") if(! -f $IPRSOURCE);
#die("ERROR: The file $TARGET does not exist\n") if(! -f $TARGET);

#parse AUGUSTUS gff3
open(FileHandler_AUGUSTUS, "< aug5.gff3");# or die "ERROR: Could not open the file $TARGET\n";
open(FileHandler_IPR, "< aug5.pep.ipr.gff3");# or die "ERROR: Could not open the file $IPRSOURCE\n";		
open FileHandler_outFILE, "> outfile.gff";

while(defined(my $AUGUSTUS_line = <FileHandler_AUGUSTUS>)){
    last if($AUGUSTUS_line =~ /\#\#FASTA/); # Ends as soon as FASTA begins
    next if($AUGUSTUS_line =~ /^\#/); # AUGUSTUS file does not follow strict gff3 guidelines, thus no # is there for next segment
    next unless($AUGUSTUS_line =~ /\tCDS\t/);

    chomp $AUGUSTUS_line;
    my @F = split("\t", $AUGUSTUS_line);
    if($F[2] eq 'CDS'){
	$_ = $F[8];
	/Parent\=/;
	$parent = $';
	$start = $F[3];
	$stop = $F[4];
	$seqid = $F[0];
	$strand = $F[6];
	$phase = $F[7];
	#print $parent;
			#parse InterProScan gff3
		while(defined(my $IPR_line = <FileHandler_IPR>)){ # Filehandler zurücksetzen !!!!! Muss von Anfagn des Files starten
		    last if($IPR_line =~ /\#\#FASTA/); # Ends as soon as FASTA begins
		    next if($IPR_line =~ /^\#/);
		    next unless($IPR_line =~ /\tprotein_match\t/);
		    chomp $IPR_line;
		    my @F = split("\t", $IPR_line);
		    if($F[0] eq $parent){
				$type = $F[1];
				$AAstart = $F[3];
				$AAstop = $F[4];
				$score = $F[5];
				$annotation = $F[8];
					my $source = "InterProScan";
					my $startAAcalc = ($AAstart*3)+$start;
					my $stopAAcalc = ($AAstop*3)+$start;
				#print FILE "$seqid\t$source\t$type\t$startAAcalc\t$stopAAcalc\t$score\t$strand\t$phase\t$annotation\n";
				print FileHandler_outFILE "$seqid\t$source\t$type\t$startAAcalc\t$stopAAcalc\t$score\t$strand\t$phase\t$annotation\n";
			}
		}
		seek FileHandler_IPR,0,0;	
	}
}
close (FileHandler_AUGUSTUS);
close (FileHandler_outFILE);
