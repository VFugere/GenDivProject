#!/usr/bin/perl -w
### This script will get FASTA sequences from each accession number from file
### perl FastaFromAccession.pl mammals_coi_dates_km.txt.taxo

##Input file with accession numbers in particular format (4th and 5th columns)
my $ACCFile=shift;
if(!($ACCFile)){die "You forgot to indicate an input file!\n"};
if(!(-f "$ACCFile")){die "Could not find $ACCFile\n"};

##Read file with species names
my $acc="";
my $fasta="";
my $newheader="";
my $previousSp="";
my $previousYr="";
#my @List=();### accessions
#my @Full=();### lines to replace in the fasta header
open(IN,$ACCFile) or die "can't open $ACCFile\n";
open(OUT,">$ACCFile.fasta") or die "can't write to $ACCFile output file\n";

while(<IN>){
	next if(/^Species/);
	my $line=$_;chomp($line);
	my @L=split(/\t/);
	
	if($L[0] ne "$previousSp"){
		close OUTfasta if($previousSp ne "");
		open(OUTfasta,">FASTA/$L[0].fasta") or die "can't write to FASTA/$L[0].fasta output file\n";
		$previousYr="";
	}
	$previousSp=$L[0];
	
	### Separate by year too
	if($L[5] ne "$previousYr"){
		close OUTfastaYr;
		open(OUTfastaYr,">FASTA/$L[0].$L[5].fasta") or die "can't write to FASTA/$L[0].$L[5].fasta output file\n";
	}
	$previousYr=$L[5];

		
	if($L[4] eq "NA"){
		$acc=$L[3];
		$newheader=join("_",">$L[0]","$acc",$L[1],$L[2],$L[5]);
		#$fasta=system "curl -s \"http://www.boldsystems.org/index.php/API_Public/sequence?ids=${acc}|${acc}\" | grep -A 1 COI"
		$fasta=`curl -s \"http://www.boldsystems.org/index.php/API_Public/sequence?ids=${acc}|${acc}\" | grep -A 1 COI`;
	}else{
		$acc=$L[4];
		$newheader=join("_",">$L[0]","$acc",$L[1],$L[2],$L[5]);
		#$fasta=system "curl -s \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text\"";
		$fasta=`curl -s \"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text\"`;
	}
	
	#### replace header with own header, and one line fasta
	my @oneline=split(/\n/,$fasta);
	my $head=shift(@oneline);
	my $seq=join("",@oneline);
	print OUT "$newheader\n$seq\n";
	print OUTfasta "$newheader\n$seq\n";
	print OUTfastaYr "$newheader\n$seq\n";
	
	#push(@List,$acc);
	#push(@Full,$newheader);
}
close IN;
