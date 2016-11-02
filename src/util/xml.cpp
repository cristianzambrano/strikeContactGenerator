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

#include "xml.h"

using namespace std;


int
find_next(const string &tag, char *value, FILE *xml_F)
{
	static const unsigned int READ_LENGTH = 501;
	static char line[READ_LENGTH];
	char *tag_found;
	while (fgets(line, READ_LENGTH, xml_F) != NULL)
	{
		strtok(line, "<>");
		tag_found = strtok(NULL, "<>");
		if (!strcmp(tag_found ,tag.c_str()))
		{
			strcpy(value, strtok(NULL, "<>"));
			return 0;
		}
	}
	return 1;
}

int
find_next_int(const string &tag, FILE *xml_F)
{
	static const unsigned int READ_LENGTH = 501;
	static char line[READ_LENGTH];
	char *value;
	while (fgets(line, READ_LENGTH, xml_F) != NULL)
	{
		strtok(line, "<>");
		value = strtok(NULL, "<>");
		if (!strcmp(value ,tag.c_str()))
		{
			return atoi(strtok(NULL, "<>"));
		}
	}
	return 0;
}



string
find_best_pdb(unsigned int min_coverage, unsigned int min_identity, const string &xml_f )
{
	FILE *xml_F = fopen(xml_f.c_str(), "r");
	const unsigned int READ_LENGTH = 201;
	char line[READ_LENGTH];

	char value[100];
	char name[20];
	char *tag;
	unsigned int query_len = 0;
	unsigned int start = 0, end = 0, identity = 0, align_len = 0;
	while (fgets(line, READ_LENGTH, xml_F) != NULL)
	{
		strtok(line, "<>");
		tag = strtok(NULL, "<>");
		if (!strcmp(tag ,"BlastOutput_query-len"))
		{
			tag = strtok(NULL, "<>");
			query_len = atoi(tag);
			break;
		}
	}

	int stop = 0;
	int coverage;
	int id;
	while (fgets(line, READ_LENGTH, xml_F) != NULL)
	{
		stop = find_next("Hit_accession", value, xml_F);
		if (stop)
			break;
		else
			strcpy(&name[0], &value[4]);
		start = find_next_int("Hsp_query-from", xml_F);
		end = find_next_int("Hsp_query-to", xml_F);
		identity = find_next_int("Hsp_identity", xml_F);
		align_len = find_next_int("Hsp_align-len", xml_F);

		coverage = (end-start)*100.0/query_len;
		id = (identity*100/align_len);
		if ((coverage > 80) && (id > 80))
			return string(name);
	}
	return "";

}



