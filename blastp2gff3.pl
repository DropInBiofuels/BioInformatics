#!/usr/bin/perl

# Johannes Kabisch
# kabisch@uni-greifswald.de
# 12.08.2014

# This script converts blastp predictions run with custom table format:
# blastp -max_target_seqs 15 -show_gis -outfmt '6 qseqid sseqid sscinames stitle pident length evalue'  -query aug5.pep -out aug5.blast -db nr -num_threads 12
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
# blastp gff3
###query seqID	subject seqID	Subject scien name	subj. title (Annotation)	perc. ident.	length	e-val.
#~ aug5.g1.t1_1	gi|50542874|ref|XP_499603.1|	Yarrowia lipolytica CLIB122	YALI0A00110p [Yarrowia lipolytica]	100.00	872	0.0
#~ aug5.g1.t1_1	gi|50545745|ref|XP_500411.1|	Yarrowia lipolytica CLIB122	YALI0B02090p [Yarrowia lipolytica]	77.11	865	0.0
#~ aug5.g1.t1_1	gi|50545932|ref|XP_500504.1|	Yarrowia lipolytica CLIB122	YALI0B04642p [Yarrowia lipolytica]	76.26	851	0.0
#~ aug5.g1.t1_1	gi|50548489|ref|XP_501714.1|	Yarrowia lipolytica CLIB122	YALI0C11253p [Yarrowia lipolytica]	72.66	867	0.0
#~ aug5.g1.t1_1	gi|50553314|ref|XP_504068.1|	Yarrowia lipolytica CLIB122	YALI0E17589p [Yarrowia lipolytica]	68.37	860	0.0
#~ aug5.g1.t1_1	gi|50557248|ref|XP_506032.1|	Yarrowia lipolytica CLIB122	YALI0F30041p [Yarrowia lipolytica]	60.91	875	0.0
#~ aug5.g1.t1_1	gi|50555622|ref|XP_505219.1|	Yarrowia lipolytica CLIB122	YALI0F09691p [Yarrowia lipolytica]	62.91	860	0.0
#~ aug5.g1.t1_1	gi|50549187|ref|XP_502064.1|	Yarrowia lipolytica CLIB122	YALI0C20823p [Yarrowia lipolytica]	61.49	875	0.0
#~ aug5.g1.t1_1	gi|50546200|ref|XP_500621.1|	Yarrowia lipolytica CLIB122	YALI0B07898p [Yarrowia lipolytica]	58.77	861	0.0
#~ aug5.g1.t1_1	gi|50556388|ref|XP_505602.1|	Yarrowia lipolytica CLIB122	YALI0F18964p [Yarrowia lipolytica]	51.80	888	0.0
#~ aug5.g1.t1_1	gi|50553458|ref|XP_504140.1|	Yarrowia lipolytica CLIB122	YALI0E19294p [Yarrowia lipolytica]	51.10	861	0.0
#~ aug5.g1.t1_1	gi|50555666|ref|XP_505241.1|	Yarrowia lipolytica CLIB122	YALI0F10241p [Yarrowia lipolytica]	52.03	765	0.0
#~ aug5.g1.t1_1	gi|50549017|ref|XP_501979.1|	Yarrowia lipolytica CLIB122	YALI0C18491p [Yarrowia lipolytica]	45.33	911	0.0
#~ aug5.g1.t1_1	gi|126139203|ref|XP_001386124.1|	Scheffersomyces stipitis CBS 6054	hypothetical protein PICST_50085 [Scheffersomyces stipitis CBS 6054]	45.22	847	0.0
#~ aug5.g1.t1_1	gi|87045969|gb|ABD17826.1|	Candida albicans	oligopeptide transporter 4 [Candida albicans]	44.25	879	0.0
#~ aug5.g2.t1_1	gi|50542876|ref|XP_499604.1|	Yarrowia lipolytica CLIB122	YALI0A00132p [Yarrowia lipolytica]	100.00	611	0.0
#~ aug5.g2.t1_1	gi|255727773|ref|XP_002548812.1|	Candida tropicalis MYA-3404	heat shock protein SSB1 [Candida tropicalis MYA-3404]	86.30	613	0.0
#~ aug5.g2.t1_1	gi|149238586|ref|XP_001525169.1|	Lodderomyces elongisporus NRRL YB-4239	heat shock protein SSB1 [Lodderomyces elongisporus NRRL YB-4239]	85.32	613	0.0
#~ aug5.g2.t1_1	gi|584390961|emb|CDK29068.1|	Kuraishia capsulata CBS 1993	unnamed protein product [Kuraishia capsulata CBS 1993]	85.97	613	0.0

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

#seqID	subject seqID (accession)	Subject scien name	subj. title (Annotation)	perc. ident.	length	e-val.
	my $seqid;
	my $accession;
	my $organism = '';
	my $identity;
	my $SUBJlength;
	my $start;
	my $stop;
	my $score;
	my $annotation;
	my $parent;
	my $strand;
	my $phase;
	my $ID;
	my $counter = 0;
	my $real_parent;
	
	my $type = "CDS";
	my $source = "blastp";
	
#parse AUGUSTUS gff3
open(FileHandler_AUGUSTUS, "< aug5.gff3");# or die "ERROR: Could not open the file $TARGET\n";
open(FileHandler_blastP, "< aug5.blast");# or die "ERROR: Could not open the file $IPRSOURCE\n";		
open FileHandler_outFILE, "> blastP_outfile.gff";

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

	#parse blastp results
		while(defined(my $blastP_line = <FileHandler_blastP>)){ # Filehandler zurücksetzen !!!!! Muss von Anfagn des Files starten
		    last if($blastP_line =~ /\#\#FASTA/); # Ends as soon as FASTA begins
			next unless($blastP_line !~ /\t$organism\t/); # skip duplicate annotations, defined by organism
		    chomp $blastP_line;
		    my @blastP = split("\t", $blastP_line);
		    #print $blastP[1];
		    if($blastP[0] eq $parent){
				$accession = $blastP[1];
				$organism= $blastP[2];
				$score = $blastP[6];
				$identity = $blastP[4];
				$annotation = $blastP[3];
				$SUBJlength = $blastP[5];
				$ID = $parent.".a".$counter; # For unique ID
				print FileHandler_outFILE "$seqid\t$source\t$type\t$start\t$stop\t$score\t$strand\t$phase\tName=$annotation;ID=$ID;Length=$SUBJlength;Perc. identy=$identity;dbxref=$accession\n";
			}
		$counter = $counter+1; # For unique ID
		}
		seek FileHandler_blastP,0,0;	
	}
}
close (FileHandler_AUGUSTUS);
close (FileHandler_blastP);
close (FileHandler_outFILE);
