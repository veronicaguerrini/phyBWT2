#!/bin/bash

INSTALL_PREPROCESSING=1
SHORT=1

#To install the following tools:
#seqtk
#gsufsort
#BCR

if [ $INSTALL_PREPROCESSING -eq 1 ]
then
	mkdir Preprocessing
	cd ./Preprocessing

	echo -e "\n****Dowloading seqtk...\n"
	git clone https://github.com/lh3/seqtk.git;
	cd seqtk
	echo -e "\n****Compiling seqtk...\n"
	make
	
	cd ..
	echo -e "\n****Dowloading gsufsort...\n"
	git clone https://github.com/felipelouza/gsufsort
	cd gsufsort
	echo -e "\n****Compiling gsufsort...\n"
	make TERMINATOR=0 DNA=1 
	
	cd ..
	echo -e "\n****Dowloading BCR...\n"
	git clone https://github.com/giovannarosone/BCR\_LCP\_GSA
	cd BCR\_LCP\_GSA
	cp ./../../BCR_Parameters.h Parameters.h
	echo -e "\n****Compiling BCR...\n"
	make
	
	cd ../..
fi


#Compile phyBWT
make clean 
if [ $SHORT -eq 1 ]
then
	echo -e "\n****Compiling phyBWT...\n"
	make
else
	echo -e "\n****Compiling phyBWT SHORT=0...\n"
	make SHORT=0
fi

echo -e "\nDone."

