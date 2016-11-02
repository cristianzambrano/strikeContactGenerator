/*
 * Contacts.h
 *
 *  Created on: Nov 13, 2010
 *      Author: Carsten Kemena
 */


/*! \file Contacts.h
	\brief Header file for a class to save contacts between protein residues.
*/


#ifndef CONTACTS_H_
#define CONTACTS_H_


#include <string>
#include <vector>
#include <map>
#include <algorithm>


#include "stdio.h"
#include "ctype.h"
#include "string.h"
#include "stdlib.h"


using namespace std;


//! A class to save protein contacts.


class Contacts
{


private:

	//---   Variables   ---
	string _seq_name; /*!< The sequence name. */
	string _pdb_f; /*!< The file of which the contacts were created. */
	char _chain; /*!< The chain which was used to create this contacts. */
	unsigned int _num_contacts; /*!< The number of contacts extracted. */
	string _sequence; /*!< The sequence belonging to the contacts. */

    //! The contacts
    /*! Fist dimension has the length of the sequence. This corresponds to the left contact. The second dimension gives the right contact.
     * The length depends on the number of contacts found for the residue in the first dimension.*/
	vector<vector <unsigned int> > _pairs;

	//---   Methods   ---


public:



	
	//--- Constructors  ---

	Contacts();
	

	/**
	* Standard Constructor
	*/
	Contacts(const char *pdb_f, char chain);

	/**
	 * Constructor giving the sequence length.
	 */
	Contacts(unsigned int seq_length, const char *pdb_f, char chain);

	/**
	 * Standard destructor
	 */
	virtual ~Contacts();



	
	//--- Methods ---


	/**
	* \brief Returns the PDB file from which the contacts where derived.
	* \return File name
	*/
	const char*
	pdb_name() const
	{
		return _pdb_f.c_str();
	}

	/**
	* \brief Returns the PDB file from which the contacts where derived.
	* \return File name
	*/
	char
	chain_name() const
	{
		return _chain;
	}

	/**
	 * \brief Sets the sequence.
	 *
	 * \param sequence The sequence belonging to the contacts.
	 */
	void set_seq(string sequence)
	{
		_sequence = sequence;
	}


	/**
	 * \brief Adds a contact.
	 *
	 * \param a frist residue.
	 * \param b second residue.
	 */
	void add_pair(unsigned int a, unsigned int b)
	{
		_pairs[a].push_back(b);
		++_num_contacts;
	}

	/**
	 * \brief The number of contacts in this object.
	 * \return The number of contacts saved.
	 */
	unsigned int num_contacts() const
	{
		return _num_contacts;
	}


	/**
	* Reads the contact file in CAO format.
	*
	* \param contact_f The file containing the contacts
	* \param min_distance The minimum distance which sould be used between contacts
	*/
	void read_contact_file(char *contact_f, int min_distance);


	/**
	 * \brief prints the contacts to standard output.
	 *
	 */
	void print()
	{
		print(NULL);
	}
	
	/**
	* \brief prints the contacts to a file.
	*
	* \param contact_f The file name.
	*/
	void print(char *contact_f);
	
	/**
	 * \brief Returns all contacts.
	 *
	 *  Fist dimension has the length of the sequence. This corresponds to the left contact. The second dimension gives the right contact. The length depends on the number of contacts found for the residue in the first dimension.
	 *
	 * \return The contacts.
	 */
	const vector<vector< unsigned int> > & get_contacts() const
	{
		return _pairs;
	}

	/**
	 * \brief Returns the sequence.
	 *
	 * \return The sequence.
	 */
	const string & get_sequence() const
	{
		return _sequence;
	}



	/**
	 * \brief Returns the sequence name.
	 *
	 * \return The sequence name.
	 */
	const string & get_seq_name() const
	{
		return _seq_name;
	}


	/**
	 * \brief Sets the sequence name.
	 *
	 * \param seq_name The sequence name.
	 */
	void set_seq_name(const char *seq_name)
	{
		_seq_name = seq_name;
	}


	/**
	 * \brief Returns the sequence length.
	 *
	 * \brief The sequence length.
	 */
	unsigned int get_seq_length() const
	{
		return _sequence.size();
	}




};




#endif /* CONTACTS_H_ */





