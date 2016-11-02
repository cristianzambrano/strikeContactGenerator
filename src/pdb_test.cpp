/*
 * main.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: Carsten Kemena
 */


#include <cstring>
#include <cmath>

// My own
#include "classes/PDB.h"
#include "classes/Contacts.h"
#include "classes/Alignment.h"
#include "classes/RestServices.h"
#include "util/xml.h"


// #include <curl/curl.h>

#ifdef _OPENMP
	#include <omp.h>
#endif



using namespace std;





int
main(int argc, char *argv[])
{
	char *pdb_f = argv[1];
	printf("pdb_f %s\n", argv[1]);
	PDB pdb;
	pdb.read_pdb(pdb_f, true);
// 	pdb.print("j");

}
