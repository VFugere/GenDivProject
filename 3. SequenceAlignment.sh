#!/bin/sh
### This script will align files

cd FASTA
for i in *.fasta;do
if [[ ! -f ../FASTA_Aligned/$i ]];then
	mafft --quiet $i |  perl -e '{while(<>){if(/^\n/){next}if(/^>/){if($l){print "\n$_"}else{print $_}}else{$l=$_;chomp($l);print $l}}}' > ../FASTA_Aligned/$i
fi
done
