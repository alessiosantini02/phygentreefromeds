#include "EDS-BWT/EDSBWTsearch.hpp"
#include "EDS-BWT/Sorting.h"
#include "EDS-BWT/Parameters.h"

#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <string.h>     // std::string, std::to_string
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>      // std::stringstream

#include <math.h>

#include <omp.h>

#include <sys/resource.h>

////////2024///////
#include <map>
#include <sdsl/bit_vectors.hpp>
#include "EDS-BWT/malloc_count/malloc_count.h"

// This function takes the position i of a dollar and the associated rank and select support (for 1) rb_1 and bsel_1
// of a compressed bitvector, computes the indexes of the dollars belonging to the starting and ending words
// of the segment on its left, and lastly returns the POSITIONS in the BWT of the characters
// which precede those dollars in the words. The positions are given as the beginning and end of the interval.
rangeElement preceding_dollars_finder(dataTypeNSeq i, rank_support_v<1> &rb_1, bit_vector::select_1_type &bsel_1){
	dataTypeNChar a;


	a = rb_1(i+1); //a is the index of the segment containing the i-th dollar

	dataTypeNChar start = bsel_1(a-1); //bsel_1(k) gives the index of the k-th 1. That is, the index of the first word of the k-th segment.
	
	
	dataTypeNChar end = bsel_1(a)-1; //same as before. We subtract 1 because we want to find the index of the last word of the (a-1)-th segment.

	#if DEBUG == 1
	std::cerr << "preceding_dollars_finder - start " << start << " end " << end << "\n";
	#endif

	rangeElement output;
	output.startPosN = start+1;
	output.endPosN = end+1;
										
	return output;
}