
/*! \file Alignment.h
	\brief Functions for reading and scoring of alignments.
*/



#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_


//STL
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <climits>

//my classes
#include "Contacts.h"
#include "Matrices.h"


//! A class to read and score alignments.
class Alignment
{

	private:
		vector<string> _names;
		vector<string> _sequences;
		unsigned int _num_sequences;
		char _aln_f[100];

		/**
		 * Returns the matches/mismatches from sequence1 to sequence 2
		 *
		 * Uses the simple needleman wunsch algorithm to produce the alignment.
		 * \param seq1 Sequence 1.
		 * \param seq2 Sequence 2.
		 * \param gap_cost The gap penalty to be used.
		 * \return The mapping of sequence 1 to sequence 2.
		*/
		map<unsigned int, unsigned int> * _gotoh_align(const string &seq1, const string &seq2);


		//Alignment reading
		int _detect_aln_format(FILE *aln_F);
		void _read_fasta_aln(FILE *aln_F);
		void _read_clustalw_aln(FILE *aln_F);
		void _read_msf_aln(FILE *aln_F);


		/**
		 * Maps the given contacts to an aligned sequence.
		 * \param aln_string The string to which the contacts should be mapped.
		 * \param contact The contacts.
		 * \param core_region_f A file name containing the core regions either BALIBASE or PREFAB format. If null the whole alignment will be evalulated.
		 * \param num_contacts will save the number of contacts.
		 */
		pair<unsigned int, unsigned int> * _map_contacts(const string &aln_string, const string &aln_seq_name, const Contacts &contact, const char *core_region_f, unsigned int &num_contacts);


		pair<unsigned int, unsigned int> *
		_map_contacts(const string &aln_string, const string &aln_seq_name, const Contacts &contact, unsigned int &num_contacts);


	public:
		/**
		* Standard constructor.
		*/
		Alignment();

		/**
		* Standard destructor.
		*/
		virtual ~Alignment();

		const string& name(unsigned int i) const
			{return _names[i];};

		const string& sequence(unsigned int i) const
		{return _sequences[i];};

		unsigned int n_seqs() const
			{return _names.size();};
		/**
		 * \brief Reads an alignment.
		 *
		 * It is able to read msf, fasta and clustalw formatas
		 * \param aln_f The alignment file.
		 */
		void read_alignment(const char *aln_f);



		/**
		 * \brief Scores alignment using
		 *
		 * \param contacts The contacts to use for the scoring.
		 * \param min_dist The minimum distance in the sequence to use.
		 * \param normalize Normalizes the score.
		 * \returns The calculated contact score value.
		 */
		double score_cs(const Contacts &contacts, unsigned int min_dist, bool normalize);


// 		double score_cao(const Contacts &contacts, unsigned int min_dist);

		/**
		 * \brief This function prints the alignment in fasta format to standard output.
		 */
		void print();

};



#endif /* ALIGNMENT_H_ */