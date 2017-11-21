//Uses Knuth–Morris–Pratt algorithm to search global string for substring
//Cpp format: Webkit: 80 character max line

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>
#include "findORFs.h"

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;


// Input is fasta sequences and fiveUTRs as grangesList
// Output is GRanges object with ranges mapped to genomic coordinates
// Instead of passing back and forth from r, do it all in C

// [[Rcpp::export]]
List get_all_ORFs_as_List(
    CharacterVector fastaSeqs,
    std::string startCodon,
    std::string stopCodon,
    bool longestORF,
    int minimumLength)
{
  vi result_index(10000);         // index of input that output belongs to
  int index = 0;
  std::vector<vi> result_value(2,vi(10000)); // output values
  
  // create result in C++
  for(int i =0 ; i < fastaSeqs.size();i++){
    
    String f = fastaSeqs[i];
    std::string fastaSeq = std::string(f);
    
    std::vector<int> ORFdef = get_all_orfs_as_vector(fastaSeq, startCodon, stopCodon,
                                        longestORF, minimumLength);
    //check if resize is needed
    
    for(int j = 0; j < ORFdef.size()/2;j++){
      result_value[0][index+j] = ORFdef[j*2] ; // set starts
      result_value[1][index+j] = ORFdef[(j*2)+1] ; // set ends
      
      result_index[index+j] = i+1;
    }
    index+=ORFdef.size()/2; // Increase counter
  }
  
  result_index.resize(index);
  result_value[0].resize(index);
  result_value[1].resize(index);
  // then return to R, e.g., allowing many orf per input
  return List::create(
    Named("index") = wrap(result_index),
    Named("orf") = wrap(result_value));
}



