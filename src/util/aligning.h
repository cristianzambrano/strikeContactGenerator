

#include <map>
#include <string>
#include <climits>


#include "../classes/Matrices.h"

typedef
struct{
	std::string seq;
	std::map<unsigned int, unsigned int> mapping;
} Merged_seq;

Merged_seq*
_merge_seqs(const string &seq1, const string &seq2);
