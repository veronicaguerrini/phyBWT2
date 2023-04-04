#!/bin/bash

####################################################################
####################################################################
#Script for building the three data structures (ebwt, cda and lcp)
#from a collection of fasta files
####################################################################
####################################################################

SHORT=1

#Fasta files path
PathDataset=$1

#Output files name
fastaName=$2

############################
InfoFile=$fastaName".txt"
FastaDataset=$fastaName".fasta"

############################
#Paths for installed tools
############################
pathseqtk="./Preprocessing/seqtk"
pathBCR="./Preprocessing/BCR_LCP_GSA"
pathgsufsort="./Preprocessing/gsufsort"
############################
############################

> $InfoFile
> $FastaDataset

echo -e "\nComputing a single fasta file..."

for file in $PathDataset/*
do
	NAME=$(basename "$file" .fasta)
	echo $NAME
	#Count string number
	nReads=$(grep ">" $file | wc -l)
	#Append F+RC
	$pathseqtk/seqtk seq -U $file  >> $FastaDataset
	$pathseqtk/seqtk seq -r -U $file >> $FastaDataset	
	#Write fileInfo
	nReads=$(($nReads*2))
	echo -e $NAME"\t"$nReads >> $InfoFile
done

echo -e "\nComputing eBWT/DA/LCP..."

if [ $SHORT -eq 1 ]
then
	#BCR for eBWT/LCP/DA
	/usr/bin/time -v $pathBCR/BCR_LCP_GSA $FastaDataset $FastaDataset > "BCR_"$(basename "$FastaDataset" .fasta)".stdout" 2> "BCR_"$(basename "$FastaDataset" .fasta)".stderr"
	rm *.len
	rm *.info
else
	#gsufsort for eBWT/LCP/DA
	/usr/bin/time -v $pathgsufsort/"gsufsort-64" $FastaDataset --bwt --lcp 4 --da 4 > "gsufsort_"$(basename "$FastaDataset" .fasta)".stdout" 2> "gsufsort_"$(basename "$FastaDataset" .fasta)".stderr"
	mv $FastaDataset".bwt" $FastaDataset".ebwt"
	mv $FastaDataset".4.da" $FastaDataset".da"
	mv $FastaDataset".4.lcp" $FastaDataset".lcp"
fi
echo "Done."

echo -e "\nComputing CDA..."
./create_cda $FastaDataset $InfoFile
echo "Done."
####
####

