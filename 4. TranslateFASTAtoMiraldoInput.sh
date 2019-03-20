#!/bin/sh
### This script will align files

###ACTG- = 12340
cd FASTA_Aligned
#for i in Vampyressa_bidens.*.fasta;do
#for i in Artibeus_cinereus.*.fasta;do
for i in *.fasta;do
	cat $i | grep -v "^>" | tr 'actgnN-' '1234000' | sed 's/\(.\{1\}\)/\1 /g' > ../InputAlignments/$i
	cat $i | grep "^>" | cut -f 4-5 -d"_" | sed 's/_/	/' > ../InputAlignments/${i%.fasta}.coords
	## if file is empty, delete it
	if [[ ! -s ../InputAlignments/$i ]];then
		rm -f ../InputAlignments/$i ../InputAlignments/${i%.fasta}.coords
	fi

done
