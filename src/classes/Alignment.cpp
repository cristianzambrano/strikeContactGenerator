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

#include "Alignment.h"


Alignment::Alignment()
{
}



Alignment::~Alignment()
{

}




map<unsigned int, unsigned int> *
Alignment::_gotoh_align(const string &seq1, const string &seq2)
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
	map<unsigned int, unsigned int> *mapping = new map<unsigned int, unsigned int>();
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
		}
		else if (mode == 'h')
		{
			if (matrix_h[i][j] != matrix_h[i][j-1] + gep)
			{
				mode = 'm';
			}
			--j;
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
			(*mapping)[i] = j;
		}
	}


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


	return mapping;
}


void
Alignment::print()
{
	for (unsigned int i = 0; i < _names.size(); ++i)
	{
		printf("%s\n%s\n", _names[i].c_str(), _sequences[i].c_str());
	}
}


void
Alignment::_read_fasta_aln(FILE *aln_F)
{
	const unsigned int READ_LENGTH = 201;
	char line[READ_LENGTH];

	_num_sequences = 0;
	while (fgets(line, READ_LENGTH, aln_F) != NULL)
	{

		if (line[0] == '>')
		{
			line[strlen(line)-1] = '\0';
			_names.push_back(&line[1]);
			_sequences.push_back("");
			++_num_sequences;
		}
		else
		{
			if (line[strlen(line)-1] == '\n')
			line[strlen(line)-1] = '\0';
			_sequences[_num_sequences-1].append(line);
		}
	}


	unsigned int j, length;
	length = _sequences[0].size();
	for (unsigned int i = 0; i < _num_sequences; ++i)
	{
		for (j = 0; j < length; ++j)
			_sequences[i][j] = toupper(_sequences[i][j]);
	}

}



void
Alignment::_read_clustalw_aln(FILE *aln_F)
{
	const unsigned int READ_LENGTH = 201;
	char line[READ_LENGTH];
	fgets(line, READ_LENGTH, aln_F);
	fgets(line, READ_LENGTH, aln_F);
	char *seq_name, *seq;
	unsigned int counter = 0;
	while (fgets(line, READ_LENGTH, aln_F) != NULL)
	{

		seq_name = strtok(line, " \n");
		seq = strtok(NULL, " \n");

		if ((line[0] != '\n') && (line[0] != ' '))
		{
			if (counter == _sequences.size())
			{
				_sequences.push_back(seq);
				_names.push_back(seq_name);
			}
			else
			{
				_sequences[counter].append(seq);
			}
			++counter;
		}
		else
			counter = 0;
	}
	_num_sequences = _sequences.size();


	unsigned int j;
	for (unsigned int i = 0; i < _num_sequences; ++i)
	{
		for (j = 0; j < _sequences[i].size(); ++j)
			if (_sequences[i][j] == '.')
				_sequences[i][j] = '-';
			else
				_sequences[i][j] = toupper(_sequences[i][j]);
	}
}



void
Alignment::_read_msf_aln(FILE *aln_F)
{
	_num_sequences = 0;
	const unsigned int READ_LENGTH = 201;
	char line[READ_LENGTH];

	char *seq_name, *seq;
	while (fgets(line, READ_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '/')
			break;
		if (line[0] == '\n')
			continue;

		seq_name = strtok(line, " \n");
		if (!strcmp(seq_name, "Name:"))
		{
			seq_name = strtok(NULL, " \n");
			++_num_sequences;
			_names.push_back(seq_name);
		}
	}

	_sequences.resize(_names.size());
	unsigned int counter = -1;
	while (fgets(line, READ_LENGTH, aln_F) != NULL)
	{

		if ((line[0] != '\n') && (line[0] != ' '))
		{
			++counter;
			seq_name = strtok(line, " \n");
			while ((seq = strtok(NULL, " \n")) != NULL)
				_sequences[counter].append(seq);
		}
		else
			counter = -1;
	}

	unsigned int j;
	for (unsigned int i = 0; i < _num_sequences; ++i)
	{
		for (j = 0; j < _sequences[i].size(); ++j)
			if (_sequences[i][j] == '.')
				_sequences[i][j] = '-';
			else
				_sequences[i][j] = toupper(_sequences[i][j]);
	}
}


int
Alignment::_detect_aln_format(FILE *aln_F)
{
	const unsigned int READ_LENGTH = 201;
	char line[READ_LENGTH];
	fgets(line, READ_LENGTH, aln_F);

	if (line[0] == '>')
		return 2; //Fasta format
	if (strstr(line, "CLUSTAL"))
		return 1;
	while (fgets(line, READ_LENGTH, aln_F) != NULL)
	{
		if (strstr(line, "MSF:"))
			return 0;
	}
	return -1;
}

void
Alignment::read_alignment(const char *aln_f)
{
	FILE *aln_F = fopen(aln_f, "r");
	if (aln_F == NULL)
	{
		fprintf(stderr, "Error: File %s could not be opened!\n", aln_f);
		exit(1);
	}
	else
	{
		int format = _detect_aln_format(aln_F);
		fseek ( aln_F, 0, SEEK_SET );
		sprintf(_aln_f, "%s", aln_f);
		switch( format )
		{
			case 0 :
				_read_msf_aln(aln_F);
				break;
			case 1 :
				_read_clustalw_aln(aln_F);
				break;
			case 2 :
				_read_fasta_aln(aln_F);
				break;
			default:
				fclose(aln_F);
				fprintf(stderr, "Error: File format has not been recognized! Please check if it is in fasta, clustal or msf format.\n");
				exit(1);
		}
	}
	fclose(aln_F);
}



pair<unsigned int, unsigned int> *
Alignment::_map_contacts(const string &aln_string, const string &aln_seq_name, const Contacts &contact, unsigned int &num_contacts)
{

	// mapping from seq to aligned seq
	unsigned int aln_length = aln_string.size();
	string gapless_temp;
	unsigned int *seq2aln = new unsigned int[aln_length];
	unsigned int i, j = -1;

	for (i = 0; i < aln_length; ++i)
	{
		if (aln_string[i] != '-')
		{
			gapless_temp.push_back(aln_string[i]);
			seq2aln[++j] = i;
		}
	}


	// mapping vom pdb 2 sequence
	map<unsigned int, unsigned int> *mapping = _gotoh_align(contact.get_sequence(), gapless_temp);
        
	num_contacts = contact.num_contacts();
        printf("Num COntacts %d\n" , num_contacts);

	pair<unsigned int, unsigned int> *contact_pairs = new pair<unsigned int, unsigned int>[num_contacts];
	const vector< vector< unsigned int > > contacts = contact.get_contacts();

	unsigned int tmp_length, seq_length = contact.get_seq_length();
	unsigned int mapped_i;

	num_contacts = 0;
	for (i = 0; i < seq_length; ++i)
	{
		if (!mapping->count(i))
			continue;
		mapped_i =seq2aln[mapping->find(i)->second];
		tmp_length = contacts[i].size();
		for (j = 0; j < tmp_length; ++j)
		{
			if (!mapping->count(contacts[i][j]))
				continue;
			contact_pairs[num_contacts].first  = mapped_i;
			contact_pairs[num_contacts].second = seq2aln[mapping->find(contacts[i][j])->second];
			++num_contacts;
		}
	}

	delete[] seq2aln;
	delete mapping;


        printf("Num COntacts Final %d \n" , num_contacts);
	return contact_pairs;
}


double
Alignment::score_cs(const Contacts &contact, unsigned int min_dist, bool normalize)
{

	// identify the sequence belonging to the contacts
	string seq_name = contact.get_seq_name();
	unsigned template_seq = 0;
	while ((_names[template_seq].compare(seq_name)) && (template_seq < _num_sequences))
	{
		++template_seq;
	}
	if (template_seq == _num_sequences)
	{
		printf("ERROR: Sequence '%s' not found in alignment.\n", seq_name.c_str());
		exit(2);
	}


	unsigned int num_contacts = 0;
        printf("template_seq %d _sequences %s\nName %s\n",template_seq, _sequences[template_seq].c_str(), _names[template_seq].c_str());
         
	pair<unsigned int, unsigned int> *contact_pairs = _map_contacts(_sequences[template_seq], _names[template_seq], contact, num_contacts);


        printf("Size %d \n",num_contacts);
        for(int k=0;k<num_contacts;k++){
            
            printf("%d %d - ",contact_pairs[k].first, contact_pairs[k].second);
        }
        
        printf("\n");
        //exit(-1);
       
        
	if ( num_contacts < 50)
	{
		fprintf(stderr, "**************   !WARNING!   **************\n");
		fprintf(stderr, "      Only few contacts avilable!\n");
		fprintf(stderr, "    Scoring might be very inaccurate\n");
		fprintf(stderr, "   # contacts: %i\n", num_contacts);
		fprintf(stderr, "   Alignment: %s\n", _aln_f);
		fprintf(stderr, "   PDB: %s   -   Chain: %c\n", contact.pdb_name(), contact.chain_name());
		fprintf(stderr, "**************   !WARNING!   **************\n");
	}


	//calculate normalize value if necessary
	char char1, char2;
	unsigned int i;
	double **cs_mat = Matrices::cs();
	double normalize_value = 0;
	unsigned int contact_counter = 0;
	if (normalize)
	{
		for (i = 0; i < num_contacts; ++i)
		{
			if (contact_pairs[i].second -contact_pairs[i].first >= min_dist)
			{
				char1 = _sequences[template_seq][contact_pairs[i].first];
				char2 = _sequences[template_seq][contact_pairs[i].second];
				if ((char1 != 'X') && (char2 != 'X') && (char1 != 'B') && (char2 != 'B') && (char1 != 'Z') && (char2 != 'Z') && (char1 != 'J') && (char2 != 'J'))
				{
					++contact_counter;
					normalize_value += cs_mat[char1-65][char2-65];
				}
			}
		}
		normalize_value /= contact_counter;
	}


        // printf("normalize_value %f\n",normalize_value);
         
	//calculate score
	unsigned int j;
	double cs_score = 0;

	contact_counter = 0;
	for (i = 0; i < _num_sequences; ++i)
	{
		if (i != template_seq)
		{

			for (j = 0; j < num_contacts; ++j)
			{
				if (contact_pairs[j].second -contact_pairs[j].first >= min_dist)
				{
					char1 = _sequences[i][contact_pairs[j].first];
					char2 = _sequences[i][contact_pairs[j].second];
					if ((char1 != 'X') && (char2 != 'X') && (char1 != 'B') && (char2 != 'B') && (char1 != 'Z') && (char2 != 'Z') && (char1 != 'J') && (char2 != 'J'))
					{
						++contact_counter;
						if ((char1 != '-') && (char2 != '-'))
							cs_score += cs_mat[char1-65][char2-65];
					}
				}
			}
		}
	}


	//  free memory
	for (unsigned int tmp =0; tmp <26; ++tmp)
		delete[] cs_mat[tmp];
	delete[] cs_mat;
	delete[] contact_pairs;


        double sc;
	// return score
	if (normalize)
		sc=(cs_score/contact_counter)/normalize_value;
	else
		sc=cs_score/contact_counter;
        
        //printf("Score %f\n",sc);
        //exit(-1);
        
        return sc;
}



/*

double
Alignment::score_cao(const Contacts &contact, unsigned int min_dist)
{

	// identify the sequence belonging to the contacts
	string seq_name = contact.get_seq_name();
	unsigned template_seq = 0;
	while ((_names[template_seq].compare(seq_name)) && (template_seq < _num_sequences))
	{
		++template_seq;
	}
	if (template_seq == _num_sequences)
	{
		printf("ERROR: Sequence '%s' not found in alignment.\n", seq_name.c_str());
		exit(2);
	}


	unsigned int num_contacts = 0;
	pair<unsigned int, unsigned int> *contact_pairs = _map_contacts(_sequences[template_seq], _names[template_seq], contact, num_contacts);


	if ( num_contacts < 50)
	{
		fprintf(stderr, "**************   !WARNING!   **************\n");
		fprintf(stderr, "      Only few contacts avilable!\n");
		fprintf(stderr, "    Scoring might be very inaccurate\n");
		fprintf(stderr, "   # contacts: %i\n", num_contacts);
		fprintf(stderr, "   Alignment: %s\n", _aln_f);
		fprintf(stderr, "   PDB: %s   -   Chain: %c\n", contact.pdb_name(), contact.chain_name());
		fprintf(stderr, "**************   !WARNING!   **************\n");
	}


	//calculate score
	unsigned int i,j;
	double cao_score = 0;
	unsigned int contact_counter = 0;
	char temp_char1, temp_char2, char1, char2;
	double ****cao_mat = Matrices::cao();
	for (i = 0; i < _num_sequences; ++i)
	{
		if (i != template_seq)
		{

			for (j = 0; j < num_contacts; ++j)
			{
				if (contact_pairs[j].second -contact_pairs[j].first >= min_dist)
				{
					temp_char1 = _sequences[template_seq][contact_pairs[j].first];
					temp_char2 = _sequences[template_seq][contact_pairs[j].second];
					char1 = _sequences[i][contact_pairs[j].first];
					char2 = _sequences[i][contact_pairs[j].second];
					if ((char1 != 'X') && (char2 != 'X') && (char1 != 'B') && (char2 != 'B') && (char1 != 'Z') && (char2 != 'Z') && (char1 != 'J') && (char2 != 'J') & (temp_char1 != 'X') && (temp_char2 != 'X') && (temp_char1 != 'B') && (temp_char2 != 'B') && (temp_char1 != 'Z') && (temp_char2 != 'Z') && (temp_char1 != 'J') && (temp_char2 != 'J'))
					{
						++contact_counter;
						if ((char1 != '-') && (char2 != '-'))
						{
							cao_score += (cao_mat[temp_char1-65][temp_char2-65][char1-65][char2-65] + cao_mat[char1-65][char2-65][temp_char1-65][temp_char2-65])/2;
						}
					}
				}
			}
		}
	}


	//  free memory
	for (unsigned int tmp1 =0; tmp1 <25; ++tmp1)
	{
		for (unsigned int tmp2 =0; tmp2 <25; ++tmp2)
		{
			for (unsigned int tmp3 =0; tmp3 <25; ++tmp3)
			{
					delete[] cao_mat[tmp1][tmp2][tmp3];
			}
			delete[] cao_mat[tmp1][tmp2];
		}
		delete[] cao_mat[tmp1];
	}
	delete[] cao_mat;
	delete[] contact_pairs;

	return cao_score/contact_counter;
}




*/

