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

/*
 * Contacts.cpp
 *
 *  Created on: Nov 13, 2010
 *      Author: Carsten Kemena
 */



#include "Contacts.h"
#include <sstream>

Contacts::Contacts()
{
	_num_contacts=0;
}


Contacts::Contacts(const char *pdb_f, char chain):_pdb_f(pdb_f),_chain(chain)
{
	_num_contacts = 0;
}


Contacts::Contacts(unsigned int seq_length, const char *pdb_f, char chain):_pdb_f(pdb_f),_chain(chain)
{
	_num_contacts = 0;
	_pairs.resize(seq_length);
}



// Contacts::Contacts(char *seq_name, char *pdb_f, char *chain, unsigned int min_distance):_seq_name(seq_name), _pdb_f(pdb_f), _chain(chain)
// {
// 	//cal perl
// 	char command[1000];
//
// 	sprintf(command, "perl /users/cn/ckemena/projects/PDB/src/PDB_TOOL.pl %s %s > tmp_cont", _pdb_f.c_str(), _chain.c_str());
// 	printf("%s\n", command);
// 	system(command);
//
// 	//read contacts
// 	read_contact_file("tmp_cont", min_distance);
// }


// Contacts::Contacts(string sequence):_sequence(sequence)
// {
// 	_num_contacts = 0;
// 	_pairs.resize(sequence.length());
// }




Contacts::~Contacts()
{
}


void
Contacts::print(char *contact_f)
{
	FILE *contact_F;
	if (contact_f != NULL)
		contact_F = fopen(contact_f, "w");
	else
		contact_F = stdout;
	unsigned int i, end = _sequence.length(), j, end2;
	fprintf(contact_F, "%s\n", _seq_name.c_str());
	fprintf(contact_F, "%s\n", _sequence.c_str());
	for (i = 0; i < end; ++i)
	{
		end2 = _pairs[i].size();
		fprintf(contact_F, "%i", i);
		for (j = 0; j < end2; ++j)
		{
			fprintf(contact_F, " %i", _pairs[i][j]);
		}
		fprintf(contact_F, "\n");
	}
	if (contact_f != NULL)
		fclose(contact_F);
}



void
Contacts::read_contact_file(char *contact_f, int min_distance)
{
	FILE *contact_F = fopen(contact_f, "r");
	if (contact_F == NULL)
	{
		fprintf(stderr, "Contact file %s could not be opened !\n", contact_f );
		exit(1);
	}
	
	const unsigned int READ_LENGTH = 10000;
	char line[READ_LENGTH];


	_seq_name.clear();
	fscanf(contact_F,"%s",line);
       _seq_name.append(line);
	//printf("Sequence Name %s\n",_seq_name.c_str());

	//Read Sequence
	_sequence.clear();
	fscanf(contact_F,"%s",line);
       _sequence.append(line);

	//printf("Sequence %s L=%d\n",_sequence.c_str(),_sequence.length());

	unsigned int seq_length = _sequence.length();
	for (unsigned int i = 0; i < seq_length; ++i)
	{
		_sequence[i] = toupper(_sequence[i]);
	}


	//Read contacts
	_pairs.clear();
	_pairs.resize(seq_length);
	string tmp;
        _num_contacts = 0;

	string item,linea;
	int a,b;
	while (fgets (line , READ_LENGTH , contact_F) != NULL )
	{
		 
		 linea = string(line); 
                 stringstream ss(linea);
		 a = -1;
	         while (getline(ss, item, ' ')) {
			if(a>=0){
			    b= atoi(item.c_str());
			    if (b-a >= min_distance){
				    add_pair(a,b);
			     }
			 }else{ //Primer item de la fila
				a=atoi(item.c_str());
			 }
			
		 }		
		 
	}
	fclose(contact_F);
}

