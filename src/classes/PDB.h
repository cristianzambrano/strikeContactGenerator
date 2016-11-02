/*
 * PDB.h
 *
 *  Created on: Nov 12, 2010
 *      Author: Carsten Kemena
 */

#ifndef PDB_H_
#define PDB_H_

// #include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include <string>
#include <vector>
#include <map>

#ifdef _OPENMP
	#include "omp.h"
#endif

#include "Contacts.h"
#include "../util/aligning.h"
// using namespace std;



/*! \file PDB.h
 \brief This file contains two classes needed to work with PDB files. The Atom class and the PDB class.
*/




//! Saves information of an atom entry.
class Atom {


//  ATOM ENTRY
public:

	//   Constructors

	/**
	 * Standard Constructors
	 */
	Atom();

	/**
	 * Constructor to initalize all values.
	 * \param chain_ID The chain ID.
	 * \param res_ID The residue ID.
	 * \param amino_accid The name of the amino accid in one letter code.
	 * \param atom_name The name of the atom.
	 * \param atom_serial_number The serial number of the atom.
	 * \param x The x-coordinate.
	 * \param y The y-coordinate.
	 * \param z The z-coordinate.
	 */
	Atom(char chain_ID, unsigned int res_ID, char amino_accid, char *atom_name, unsigned int atom_serial_number, double x, double y, double z);
	virtual ~Atom();


	//   Methods
	/**
	 * Computes the euclidian distance between two atoms.
	 * \param other_atom The second atom.
	 */
	double compute_distance(const Atom &other_atom)
	{
		return sqrt((_x-other_atom._x)*(_x-other_atom._x) + (_y-other_atom._y)*(_y-other_atom._y) + (_z-other_atom._z)*(_z-other_atom._z));
	}


	/**
	 * Retuns the radius of the atom.
	 * \return The radius of the atom.
	 */
	double radius();


private:
	char _chain_ID;
	unsigned int  _res_ID;
	char _amino_accid;
	char _atom_name[5];
	unsigned int  _atom_serial_number;
	double _x, _y, _z;
};




class PDB_chain
{
public:

	char _id;
	std::string _seq_res_string;
	std::string _atom_res_string;
	std::vector<std::map< std::string, Atom> > _atom_list;

	PDB_chain()
	{};
	PDB_chain(char id, unsigned int reserve);
	virtual ~PDB_chain();

	void push_back_seq_res(char amino_accid)
	{
		_seq_res_string.push_back(amino_accid);
	}

	void push_back_atom_res(char amino_accid)
	{
		_atom_res_string.push_back(amino_accid);
	}

};



//PDB FILE


//! Class to save a PDB file.
class PDB {

private:

	std::map<char, PDB_chain>::iterator _first_chain;
	std::map<char, PDB_chain> _pdb_chains;
	unsigned int _num_chains;
// 	map<char, string> _sequences;
// 	map<char, vector<map< string, Atom> > > _atom_list;

	map<string, char> _aa;
	char *_pdb_f;

	char
	_aa3_1(char* amino_accid_name)
	{
		if (_aa.count(amino_accid_name))
			return _aa[amino_accid_name];
		else
			return '#';
	}




public:
	/**
	 * Standard Constructor
	 */
	PDB();

	/**
	 *  Standard Destructor
	 */
	virtual ~PDB();

	/**
	 * Reads a pdb file and saves the content in the object.
	 *
	 * \param pdb_f The name of the pdb_file.
	 * \param ignore_H If set to true hydrogen atoms will be ignored.
	 */
	void read_pdb(const char *pdb_f, bool ignore_H);

	/**
	 * Calculates the contacts inside a given chain.
	 *
	 * \param chain The chain to use.
	 * \param min_dist The minimum distance which should be allowed in the sequence between two residues.
	 * \param seq_name Name of the sequence.
	 * \return A Contacts object with all found contacts.
	 */
	Contacts* calculate_contacts(char chain, double min_dist, char *seq_name);

};

#endif /* PDB_H_ */
