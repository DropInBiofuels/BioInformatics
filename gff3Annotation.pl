#!/usr/bin/perl

# Johannes Kabisch
# kabisch@uni-greifswald.de
# 23.03.2015

#      This script takes a MAKER gff and annotates the protein_match with real annotations based on gi numbers
# Example input:
###
#
# gff3out

use warnings;
use strict;

#my $usage = "
#Usage:
#      ipr2gff3 <iprscan.gff> <AUGUSTUS.gff>
#

#~ if(@ARGV < 2){
    #~ print $usage;
    #~ exit(1);
#~ }

	my $line;
	my @front_gff;
	my $Annotation_front;
	my $GI_Number;
	my %annotation_hash = ();
	my $result_type;
	my $Annotation;
	my $DBXref;
#parse AUGUSTUS gff3
open(FileHandler_MAKERGFF, "< testset_nr.gff");# or die "ERROR: Could not open the file $TARGET\n";
open(FileHandler_nrAnno, "< fungi.annotations.unique");# or die "ERROR: Could not open the file $IPRSOURCE\n";		
open(FileHandler_outFILE, "> Ccurvatus_anno_nr.gff");

# read annotation set into hash, key = gi, value = annotation
while(defined(my $nrAnno_line = <FileHandler_nrAnno>)){ # Filehandler zurücksetzen !!!!! Muss von Anfagn des Files starten
		chomp $nrAnno_line;
		my @nrAnno = split(/\|/,$nrAnno_line);
		$annotation_hash{$nrAnno[1]} = $nrAnno[4];
	}

while(defined(my $MAKERGFF_line = <FileHandler_MAKERGFF>)){
    #next if($MAKERGFF_line !~ /.+\tblastx\t/); # next if there is no blastx in the line not sure if the regex .+ for anything makes sense here
    if($MAKERGFF_line !~ /\tprotein_match\t/){
		print FileHandler_outFILE $MAKERGFF_line;
		} # next if there is no blastx in the line not sure if the regex .+ for anything makes sense here
    #next unless($MAKERGFF_line =~ /\tprotein_match\t/); # 
    else{
    chomp $MAKERGFF_line;
    my @MAKERGFF = split("\t", $MAKERGFF_line);
   	$line = $MAKERGFF[8];
	$line =~ /Name=gi\|(.*?)\|/ ; # REGEX das nur gi hinter Namen ausgelesen wird. 
	$GI_Number = $1;
	@front_gff = @MAKERGFF[0..7];
	$Annotation_front = $MAKERGFF[8];
	$Annotation_front =~ /^(.*?)\;Name=gi/;
	$Annotation_front = $1;
				print FileHandler_outFILE join("\t",@front_gff,"$Annotation_front;Name=$annotation_hash{$GI_Number};Dbxref=$GI_Number\n");
			}
		}
close (FileHandler_MAKERGFF);
close (FileHandler_nrAnno);
close (FileHandler_outFILE);
