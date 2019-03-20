#!/bin/sh
### This script will get collection date using accession number from file
### sh TaxIDFromAccession.sh mammals_coi.csv


###Loop through BOLD number
ACCs=`cat mammals_coi.txt | awk -F"," '{{print $4}}' | uniq`
for ACC in $ACCs;do
date=""
if  [[ $ACC != [0-9]* ]];then
	date=`curl -s "http://www.boldsystems.org/index.php/API_Public/specimen?ids={$ACC}|{$ACC}" | grep collectiondate | cut -f 2 -d">" | cut -f 1 -d"<"`
fi
echo $ACC"	"$date >> mammals_coi.txt.bold.out
done


###Loop through NCBI number
ACCs=`cat mammals_coi.txt | awk -F"," '{{print $5}}' | uniq`
for ACC in $ACCs;do
	date=""
	date=`curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={$ACC}&rettype=default&retmode=text" | grep -A 1 collection | grep -v collection | cut -f 2 -d"\""`
	echo $ACC"	"$date >> mammals_coi.txt.ncbi.out
done
