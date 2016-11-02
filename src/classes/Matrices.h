

#ifndef MATRICES_H_
#define MATRICES_H_



using namespace std;


/*! \file Matrices.h
\brief This file contains The matrix class which provides several matrices.
*/


//! A class with static functions to provide different scoring matrices.
class Matrices {
	
	
	public:

	private:
		

		
	public:
		/**
		 * Gives the Blosum62 matrix.
		 * \returns A two dimensional matrix giving the blosum values.
		 */
		static int ** blosum62();

		/**
		 * Gives the contact score matrix.
		 * \return A two dimensional matrix giving the contact values.
		 */
		static double ** cs();

		/**
		* Gives the CAO score matrix.
		* \return A four dimensional matrix giving the contact values.
		*/
		static double **** cao();
};




#endif /* MATRICES_H_ */