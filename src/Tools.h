/*
 * Basic tools and defines
 */

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <numeric>
#include <time.h>
#include <bitset>         // std::bitset
#include <deque>	  //std::deque
//Tree.hh library
#include "../external/tree.hh"
#include "../external/tree_util.hh"

#define ALF 6
#define SIZE_BITSET 64
#define BUFFERSIZE 1048576

#if SHORT
	#define MAX_LCP 255
	#define TERM '$'
#else
	#define MAX_LCP 4294967295
	#define TERM '#'
#endif


using namespace std;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#if SHORT
	#define dataTypelenSeq uchar    //---> LCP (length of the sequences)
#else
	#define dataTypelenSeq uint    //---> LCP (length of the sequences)
#endif

#define dataTypedimAlpha uchar  //---> eBWT (size of the alphabet) -- biologic case 6 symbols ($,A,C,G,N,T)

#define dataTypeNSeq uchar

#define dataTypeNumChar 1       //---> Length of the list of sorted suffixes, i.e. number of characters in the input fasta file. (If =0 -> uint)

#if dataTypeNumChar == 1
#   define dataTypeNChar ulong
#else
#   define dataTypeNChar uint
#endif

typedef std::bitset<SIZE_BITSET> type_node;
typedef std::vector<type_node> nodes;
typedef std::pair<nodes,tree<type_node>::iterator> job;
typedef std::pair<bitset<SIZE_BITSET>,dataTypeNChar> type_entry;

bool cmp_function (type_entry &i,type_entry &j) 
{ 
	return (i.second<j.second); 
}

void time_start(time_t *t_time, clock_t *c_clock){
    
	*t_time = time(NULL);
	*c_clock =  clock();
}

double time_stop(time_t t_time, clock_t c_clock){
    
	double aux1 = (clock() - c_clock) / (double)(CLOCKS_PER_SEC);
	//double aux2 = difftime (time(NULL),t_time);
	
	return aux1;
}
#endif
