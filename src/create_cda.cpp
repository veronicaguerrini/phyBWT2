#include "Tools.h"

#define BUFFERSIZE 1048576

int main(int argc, char **argv) {
	 
	if( argc < 3 )
    {
      std::cerr << "Error usage: " << argv[0] << " fastaFile fileInfo" << std::endl;
      exit(1);
    }
	
	string fastaFile = argv[1];
	string fileInfo = argv[2];
	
	std::ifstream in_list;
	in_list.open(fileInfo.c_str(), std::ifstream::in);
	if (!in_list.is_open()){
		std::cerr << "Error opening file " << fileInfo << "." << std::endl;
		exit (EXIT_FAILURE);
	}
	
    dataTypeNSeq numData=0;
	std::vector<uint> range;
	range.push_back(0);
	
	while(getline(in_list,fileInfo,'\t')){
		getline(in_list,fileInfo,'\n');
		range.push_back(std::stoi(fileInfo));
		numData++;
		range[numData]+=range[numData-1];
	}
    
	cout << "Intervals:" << endl;
	for(dataTypeNSeq c=1; c<=numData; c++)
		cout << range[c-1] << "," << range[c] << endl;
	
	//Create CDA
    string fnda, fncda;
    fnda = fastaFile + ".da\0";
    fncda = fastaFile + ".cda\0";
	
	FILE *InFileDA = fopen(fnda.c_str(), "rb");
    if (InFileDA==NULL) {
        std::cerr << "Error opening "  << fnda <<  "." << std::endl;
        exit (EXIT_FAILURE);
    }
	
    FILE *InFileCDA = fopen(fncda.c_str(), "wb");
    if (InFileCDA==NULL) {
        std::cerr << "Error opening "  << fncda <<  "." << std::endl;
        exit (EXIT_FAILURE);
    }
	
	uint* buffer= new uint[BUFFERSIZE];
	ulong numread = fread(buffer,sizeof(uint),BUFFERSIZE,InFileDA);
	
	while(numread>0){
		for(uint i=0; i<numread; i++){
			dataTypeNSeq j=0;
			bool written=false;
			while( j<numData && !written){
				if(buffer[i]>=range[j] && buffer[i]<range[j+1]){
					fwrite(&j, sizeof(dataTypeNSeq), 1, InFileCDA);
					written=true;
				}
				j++;
			}
			assert(written);
		}
		numread = fread(buffer,sizeof(uint),BUFFERSIZE,InFileDA);
    }
	
    fclose(InFileDA);
    fclose(InFileCDA);
	std::cout << "file cda: " << fncda <<  "." << std::endl;
    
	return 0;
}

