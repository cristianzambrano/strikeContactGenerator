/**
*   This file is part of STRIKE.
*
*   STRIKE is free software: you can redistribute it and/or modify
*   it under the terms of the GNU Lesser General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   STRIKE is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU Lesser General Public License for more details.
*
*   You should have received a copy of the GNU Lesser General Public License
*   along with STRIKE.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#include "aligning.h"


using namespace std;


Merged_seq*
_merge_seqs(const string &seq1, const string &seq2)
{

	unsigned int len1 = seq1.length() + 1;
	unsigned int len2 = seq2.length() + 1;
	unsigned int i, j;

	int **matrix_m = new int*[len1];
	int **matrix_h = new int*[len1];
	int **matrix_v = new int*[len1];
	for (i = 0; i < len1; ++i)
	{
		matrix_m[i] = new int[len2];
		matrix_h[i] = new int[len2];
		matrix_v[i] = new int[len2];
	}


	//  set scoring values
	int **blosum62 = Matrices::blosum62();
	int gop = -11;
	int gep = -1;
	int tmp_value;


	//  calculate Matrix
	matrix_m[0][0] = 0;
	matrix_v[0][0] = INT_MIN;
	matrix_h[0][0] = INT_MIN;

	for (i=1; i < len1; i++)
	{
		matrix_v[i][0] = matrix_m[0][0] + gop + i * gep;
		matrix_m[i][0] = matrix_v[i][0];
		matrix_h[i][0] = INT_MIN;
	}
	for (j=1; j < len2; j++)
	{
		matrix_h[0][j] = matrix_m[0][0] + gop + j * gep;
		matrix_m[0][j] = matrix_h[0][j];
		matrix_v[0][j] = INT_MIN;
		for (i=1; i < len1; i++)
		{
			matrix_h[i][j] = max(matrix_m[i][j-1] + gop, matrix_h[i][j-1]) + gep;
			matrix_v[i][j] = max(matrix_m[i-1][j] + gop, matrix_v[i-1][j]) + gep;
			tmp_value      = max(matrix_h[i][j], matrix_v[i][j]);
// 			printf("%c %c\n", seq1[i-1], seq2[i-1]);
			matrix_m[i][j] = max(tmp_value, matrix_m[i-1][j-1]+ blosum62[seq1[i-1]-65][seq2[j-1]-65]);
		}
	}
	for (unsigned int z=0; z <26; ++z)
		delete[] blosum62[z];
	delete[] blosum62;

	//  backtracking
	--i;
	--j;
	Merged_seq *merged_seq = new Merged_seq();
	merged_seq->seq.reserve(len1);


	char mode = 'm';
	if (matrix_h[i][j] > matrix_m[i][j])
		mode = 'h';
	if (matrix_v[i][j] > matrix_m[i][j])
		mode = 'v';

	while ((i > 0) && (j > 0))
	{
		if (mode == 'v')
		{
			if (matrix_v[i][j] != matrix_v[i-1][j] + gep)
			{
				mode = 'm';
			}
			--i;
			merged_seq->seq.push_back(seq1[i]);
		}
		else if (mode == 'h')
		{
			if (matrix_h[i][j] != matrix_h[i][j-1] + gep)
			{
				mode = 'm';
			}
			--j;
			merged_seq->seq.push_back(seq2[j]);
		}
		else if (mode == 'm')
		{
			if (matrix_h[i][j] == matrix_m[i][j] )
			{
				mode = 'h';
				continue;
			}
			if (matrix_v[i][j] == matrix_m[i][j] )
			{
				mode = 'v';
				continue;
			}
			--j;
			--i;
			merged_seq->seq.push_back(seq1[i]);
			merged_seq->mapping[i] = j;
		}
	}

	while (j > 0)
		merged_seq->seq.push_back(seq2[--j]);
	while (i > 0)
		merged_seq->seq.push_back(seq2[--i]);

	//  free memory
	for (i = 0; i < len1; ++i)
	{
		delete[] matrix_m[i];
		delete[] matrix_h[i];
		delete[] matrix_v[i];
	}
	delete[] matrix_m;
	delete[] matrix_h;
	delete[] matrix_v;

	unsigned int half = merged_seq->seq.size()/2;
	j = merged_seq->seq.size()-1;
	char tmp;

	for (i = 0; i < half; ++i)
	{
		tmp = merged_seq->seq[i];
		merged_seq->seq[i] = merged_seq->seq[j];
		merged_seq->seq[j] = tmp;
		--j;
	}
	return merged_seq;
}
