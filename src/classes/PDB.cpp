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
 * PDB.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: Carsten Kemena
 */




#include "PDB.h"



using namespace std;

Atom::Atom()
{

}


Atom::Atom(char chain_ID, unsigned int res_ID, char amino_accid, char *atom_name, unsigned int atom_serial_number, double x, double y, double z):_chain_ID(chain_ID),_res_ID(res_ID),_amino_accid(amino_accid),_atom_serial_number(atom_serial_number),_x(x),_y(y),_z(z)
{
	sprintf(_atom_name, "%s", atom_name);
}


Atom::~Atom()
{

}



double
Atom::radius()
{
	if (_atom_name[0] == 'O')
		return 1.40;
	else if (_atom_name[0] == 'N')
	{
		if (_atom_name[1] == 'Z')
			return 1.5;
		else
			return 1.65;
	}
	else if (_atom_name[0] == 'C')
	{
		if (_atom_name[1] == '\0')
			return 1.76;
		else
		{
			if ((_amino_accid == 'R') || (_amino_accid == 'N') || (_amino_accid == 'D') || (_amino_accid == 'E') || (_amino_accid == 'Q') || (_amino_accid == 'H') || (_amino_accid == 'F') || (_amino_accid == 'W') || (_amino_accid == 'Y') || (_amino_accid == 'B') || (_amino_accid == 'Z'))
			{
				if ( (_atom_name[1] == 'A') || (_atom_name[1] == 'B') || ( (_atom_name[1] == 'G') && ((_amino_accid == 'R') || (_amino_accid == 'E') || (_amino_accid =='Q'))) || ((_atom_name[1] == 'D') && (_amino_accid == 'R' )))
					return 1.87;
				else
					return 1.76;
			}
			else
				return 1.87;
		}
	}
	else if (_atom_name[0] == 'S')
		return 1.85;
	else if (_atom_name[0] == 'A')
		return 1.5;
	else if (_atom_name[0] == 'E')
		return 1.9;
	else if (_atom_name[0] == 'H')
		return 1.0;
	else
		return 0;

}


string
aa1_3(char amino_accid_name)
{
	switch (amino_accid_name)
	{
		case 'A': return "ALA";
		case 'R': return "ARG";
		case 'N': return "ASN";
		case 'D': return "ASP";
		case 'B': return "ASX";
		case 'C': return "CYS";
		case 'E': return "GLU";
		case 'Q': return "GLN";
		case 'Z': return "GLX";
		case 'G': return "GLY";
		case 'H': return "HIS";
		case 'I': return "ILE";
		case 'L': return "LEU";
		case 'K': return "LYS";
		case 'M': return "MET";
		case 'F': return "PHE";
		case 'P': return "PRO";
		case 'S': return "SER";
		case 'T': return "THR";
		case 'W': return "TRP";
		case 'Y': return "TYR";
		case 'V': return "VAL";
		case 'X': return "XAA";
		case 'U': return "CSE";
		default: return "NNN";
	}
}


PDB_chain::PDB_chain(char id, unsigned int reserve):_id(id)
{
	_seq_res_string.reserve(reserve);
	_atom_res_string.reserve(reserve);
}


PDB_chain::~PDB_chain()
{
}



PDB::PDB() {
	_aa["ALA"] = 'A';
	_aa["ARG"] = 'R';
	_aa["ASN"] = 'N';
	_aa["ASP"] = 'D';
	_aa["ASX"] = 'B';
	_aa["CYS"] = 'C';
	_aa["GLU"] = 'E';
	_aa["GLN"] = 'Q';
	_aa["GLX"] = 'Z';
	_aa["GLY"] = 'G';
	_aa["HIS"] = 'H';
	_aa["ILE"] = 'I';
	_aa["LEU"] = 'L';
	_aa["LYS"] = 'K';
	_aa["MET"] = 'M';
	_aa["PHE"] = 'F';
	_aa["PRO"] = 'P';
	_aa["SER"] = 'S';
	_aa["THR"] = 'T';
	_aa["TRP"] = 'W';
	_aa["TYR"] = 'Y';
	_aa["VAL"] = 'V';
	_aa["XAA"] = 'X';
	_aa["CSE"] = 'U';
}

PDB::~PDB() {
	delete[] _pdb_f;
}


void
PDB::read_pdb(const char *pdb_f, bool ignore_H)
{
	_num_chains = 0;
	char last_chain = ' ';
	FILE *pdb_F = fopen(pdb_f, "r");
	if (pdb_F == NULL)
	{
		fprintf(stderr, "ERROR: Could not open PDB file: %s\n", pdb_f);
		exit(1);
	}
	const char *tmp_name = strrchr(pdb_f,  '/');
	if (tmp_name == NULL)
		tmp_name = pdb_f;
	else
		++tmp_name;
	_pdb_f = new char[strlen(tmp_name)+1];
	sprintf(_pdb_f, tmp_name);
	int const LINE_LENGTH = 201;
	char line[LINE_LENGTH];
	int file_res_ID, old_res_ID = -9999;
	char iCode, old_iCode='%';
	unsigned int res_ID= 0, atom_serial_number = 0;
	char amino_accid;
	char chain_ID;
	char *atom_name = new char[5];

	double x, y, z;
	unsigned int i, j;
	char *tmp_aa;

	pair<map<char, PDB_chain>::iterator, bool> inserted;
	map<char, PDB_chain>::iterator chain_it;
	while (fgets(line, LINE_LENGTH, pdb_F) != NULL)
	{
		if (!strncmp( "ENDMDL", line, 6))
		{
			while (fgets(line, LINE_LENGTH, pdb_F) != NULL)
			{
				if (!strncmp( "ATOM", line, 6))
					if (last_chain != line[11])
					{
						break;
					}
			}
		}

		if (!strncmp( "SEQRES", line, 6))
		{
			if (last_chain != line[11])
			{
				last_chain = line[11];
				++_num_chains;

				inserted = _pdb_chains.insert(pair<char, PDB_chain>(last_chain, PDB_chain(last_chain, atoi(&line[14]))));
				chain_it = inserted.first;
				if (_num_chains == 1)
					_first_chain = chain_it;
			}

			strtok(&line[16], " \n");
			while ((tmp_aa = strtok(NULL, " \n")) != NULL)
			{
				if ((amino_accid = _aa3_1(tmp_aa)) != '#')
					chain_it->second.push_back_seq_res(amino_accid);
				else
					continue;
			}
		}
		else
		if (!strncmp( "ATOM", line, 4))
		{

			if (line[21] != last_chain)
			{
				last_chain = line[21];
				chain_it = _pdb_chains.find(last_chain);
				//Cases with ATOM-FIELD but no SEQRES FIELD are ignored!
				if  (chain_it == _pdb_chains.end())
				{
					fprintf(stderr, "Warning: Chain %c does not have an SEQRES field and is ignored!\n", last_chain);
					while (fgets(line, LINE_LENGTH, pdb_F) != NULL)
					{
						if (!strncmp( "TER", line, 3))
							break;
					}
					continue;
				}
				old_iCode='%';
				old_res_ID = -999;
				res_ID = 0;
			}

			j = -1;
			for (i = 12; i < 16; ++i)
			{
				if (line[i] != ' ')
					atom_name[++j] = line[i];
			}
			atom_name[++j] = '\0';

			if ((ignore_H) && (atom_name[0] == 'H'))
				continue;

			line[11] = '\0';
			++atom_serial_number;

			line[20] = '\0';
			amino_accid = _aa3_1(&line[17]);
			if (amino_accid == '#')
				continue;
			chain_ID = line[21];

			iCode = line[26];

			line[26] = '\0';
			file_res_ID = atoi(&line[22]);

			if ((old_res_ID != file_res_ID) || (iCode != old_iCode))
			{
				++res_ID;
				old_res_ID = file_res_ID;
				old_iCode = iCode;
// 				printf("last_chain: %c %c %s\n", last_chain, amino_accid, line);
				chain_it->second._atom_list.push_back(map<string, Atom>());
				chain_it->second.push_back_atom_res(amino_accid);
			}

			char tmp;
			tmp = line[38];
			line[38] = '\0';
			x = atof(&line[30]);
			line[38] = tmp;

			tmp = line[46];
			line[46] = '\0';
			y = atof(&line[38]);
			line[46] = tmp;

			tmp = line[54];
			line[54] = '\0';
			z = atof(&line[46]);
			line[54] = tmp;

			chain_it->second._atom_list[res_ID-1][atom_name]=Atom(chain_ID, res_ID, amino_accid, atom_name, atom_serial_number, x, y, z);
		}
	}

	fclose(pdb_F);
	delete[] atom_name;
}


Contacts *
PDB::calculate_contacts(char chain, double min_dist, char *seq_name)
{
	if (chain == '#')
		chain = _first_chain->first;

	map<char, PDB_chain>::iterator it;
	if ((it = _pdb_chains.find(chain)) == _pdb_chains.end())
		return NULL;

	vector< map<string, Atom> > &atoms = it->second._atom_list;
	const double r_probe=1.4 * 2;
	double ca_dist, dist, lim;
	unsigned int num_residues = atoms.size();
	map<string, Atom>::iterator it1, it2, end1, end2;
	bool in_contact = false;

	Merged_seq* merged_seq  = _merge_seqs(it->second._atom_res_string, it->second._seq_res_string);
	Contacts* contacts = new Contacts(merged_seq->seq.size(), _pdb_f, chain);
	contacts->set_seq_name(seq_name);
	contacts->set_seq(merged_seq->seq);
	unsigned int res_id_A, res_id_B, shifted_A, shifted_B;

	#pragma omp parallel for shared(num_residues, min_dist, atoms, contacts) private(shifted_A, shifted_B, res_id_A, it1, it2, end1, end2, in_contact, ca_dist, dist, lim, res_id_B) schedule(dynamic, 10)
	for (res_id_A = 0; res_id_A < num_residues; ++res_id_A)
	{
		shifted_A = merged_seq->mapping[res_id_A];
		for (res_id_B = res_id_A; res_id_B < num_residues; ++res_id_B)
		{
			shifted_B = merged_seq->mapping[res_id_B];
			if ((shifted_B - shifted_A) < min_dist)
				continue;
			if ((atoms[res_id_A].count("CA")) && (atoms[res_id_B].count("CA")))
				ca_dist = atoms[res_id_A]["CA"].compute_distance(atoms[res_id_B]["CA"]);
			else
				ca_dist = 0;
			if (ca_dist < 22.5)
			{
				end1 = atoms[res_id_A].end();
				for ( it1=atoms[res_id_A].begin() ; it1 != end1; ++it1 )
				{
					if (!(it1->first.compare("C")) || (!it1->first.compare("CA")) || (!it1->first.compare("N")) || (!it1->first.compare("O")) || (!it1->first.compare("OXT")))
						continue;

					end2 = atoms[res_id_B].end();
					for ( it2=atoms[res_id_B].begin() ; it2 != end2; ++it2 )
					{
						if ((!it2->first.compare("C")) || (!it2->first.compare("CA")) || (!it2->first.compare("N")) || (!it2->first.compare("O")) || (!it2->first.compare("OXT")))
							continue;

						dist = it1->second.compute_distance(it2->second);
						lim= it1->second.radius() + it2->second.radius() + r_probe;

						// Once a contact is found no other atoms between these two residues have to be examined.
						if(dist < lim) {
							#pragma omp critical (add_pair)
							{
								contacts->add_pair(shifted_A, shifted_B);
							}
							in_contact = true;
							break;
						}
					}
					if (in_contact)
					{
						in_contact = false;
						break;
					}
				}
			}
		}
	}
	
	return contacts;
}

