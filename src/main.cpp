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
 * main.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: Carsten Kemena
 */


#include <cstring>
#include <cmath>

#include "classes/PDB.h"
#include "classes/Contacts.h"

using namespace std;

int main(int argc, char *argv[])
{

		PDB pdb;	
		pdb.read_pdb(argv[2], true);

		char contact_filename[1000];
		strcpy(contact_filename,argv[4]); //Path
		strcat(contact_filename,"/");
		strcat(contact_filename,argv[1]); //SequenceName
		strcat (contact_filename,".contacts"); 

		Contacts *contactsg = pdb.calculate_contacts(argv[3][0], 5.0, argv[1]);
		if (contactsg != NULL){
		   contactsg->print(contact_filename);
		   printf("1");
		}else
		   printf("ERROR Generating Contacts for %s \n",argv[1]);

		delete contactsg;
}


