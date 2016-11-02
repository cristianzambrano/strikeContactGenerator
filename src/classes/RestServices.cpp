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

#include "RestServices.h"

using namespace std;


// EBI_Rest::EBI_Rest(const string &server,
// 				   const string &email):_server(server)
// {
//
//
// }

// EBI_Rest::~EBI_Rest()
// {
// }


// unsigned int EBI_Rest::_max_jobs = 25;
// unsigned int EBI_Rest::_running_jobs = 0;

std::string EBI_Rest::_email="";
void EBI_Rest::lib_start(const string &email)
{

	unsigned int pos = email.find('@');
	unsigned int email_len = email.size();
	_email.reserve(email_len+2);
	_email.append(email, 0, pos);
	_email.append("%40");
	_email.append(email, pos+1, email_len-(pos+1));
	curl_global_init(CURL_GLOBAL_ALL);
}



void writefunc(void *ptr, size_t size, size_t nmemb, string *s)
{
	s->assign((char*)ptr);
}


short
EBI_Rest::get_job_status(const string &job_id)
{
	CURL *easyhandle = curl_easy_init();
	string status;
	string url("http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/");
	url.append(job_id);
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  &status);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	if (status == "RUNNING")
		return 0;
	if (status == "FINISHED")
		return 1;
	else
		return -1;
}


void
EBI_Rest::get_job_result(const string &job_id,
						 const string &out_f)
{
	short stat;
	FILE *out_F = fopen(out_f.c_str(), "w");
	while (!(stat = get_job_status(job_id)))
	{
		sleep(5);
	}
	if (stat == -1)
	{
		fclose(out_F);
		return;
	}

	CURL* easyhandle = curl_easy_init();
	string status;
	string url("http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/");
	url.append(job_id);
	url.append("/out");
// 		printf("%s\n", job_id.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
// 	curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  out_F);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	fclose(out_F);
}


void
EBI_Rest::submit_blast_job(const std::string &program,
						   const std::string &type,
						   const std::string &database,
						   double ecut,
						   const std::string &name,
						   const std::string &seq,
						   std::string &job_id)
{
	char exp[10];
	sprintf(exp, "%.1f",ecut);
	string request;
	request.reserve(seq.size()+200);
	request.append("sequence=%3E");
	request.append(name);
	request.append("%0A");
	request.append(seq);
	request.append("&database=");
	request.append(database);
	request.append("&stype=");
	request.append(type);
	request.append("&exp=");
	request.append(exp);
	request.append("&align=7");
	request.append("&program=");
	request.append(program);
	request.append("&email=");
	request.append(_email);
// 	printf("%s\n", request.c_str());
	CURL* easyhandle = curl_easy_init();

	if(easyhandle) {
		curl_easy_setopt(easyhandle, CURLOPT_POSTFIELDS, request.c_str());
		curl_easy_setopt(easyhandle, CURLOPT_WRITEFUNCTION, writefunc);
		curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  &job_id);
		curl_easy_setopt(easyhandle, CURLOPT_URL, "http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/");
		curl_easy_perform(easyhandle);
		curl_easy_cleanup(easyhandle);
	}
}


void PDB_Rest::get_pdb_file(const string &pdb_id, const string &out_f)
{

	FILE *out_F = fopen(out_f.c_str(), "w");
	CURL* easyhandle = curl_easy_init();
	string url("http://www.pdb.org/pdb/files/");
	url.append(pdb_id);
	url.resize(url.size()-2);
	url.append(".pdb");
// 	printf("%s\n", url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_URL, url.c_str());
	curl_easy_setopt(easyhandle, CURLOPT_WRITEDATA,  out_F);
	curl_easy_perform(easyhandle);
	curl_easy_cleanup(easyhandle);
	fclose(out_F);
}



void
get_pdb_files(char *email, char *pdb_dir, const Alignment &aln, char *aln_f )
{
	EBI_Rest::lib_start(email);
	unsigned int n_seqs = aln.n_seqs();
	string template_f = pdb_dir;
	if (template_f == ".")
		template_f = "";
	char *pos =strrchr(aln_f, '/');

	if (pos == NULL)
		template_f.append(aln_f);
	else
		template_f.append(pos+1);
	template_f.append("_template.txt");

	FILE *template_F = fopen(template_f.c_str(), "w");
	if (template_F == NULL)
	{
		printf("Could not write template file %s\n", template_f.c_str());
		exit(0);
	}
	for (unsigned int i = 0; i < n_seqs; ++i)
	{
		string job_id;
		string tmp_out(pdb_dir);
		tmp_out.push_back('/');
		printf("Getting Structure for %s\n", aln.name(i).c_str());
		tmp_out.append(aln.name(i));
		tmp_out.append(".blast");
		const string &seq = aln.sequence(i);
		string ungapped_seq;
		ungapped_seq.reserve(seq.size());

		for (unsigned j = 0; j < seq.size(); ++j)
		{
			if (seq[j] != '-')
				ungapped_seq.push_back(seq[j]);
		}

		EBI_Rest::submit_blast_job("blastp", "protein", "pdb", 1.0, aln.name(i), ungapped_seq, job_id);
		EBI_Rest::get_job_result(job_id, tmp_out);
		string pdb_name = find_best_pdb(1, 1, tmp_out );
		if (pdb_name != "")
		{
			string pdb_file = pdb_name;
			pdb_file.append(".pdb");
			PDB_Rest::get_pdb_file(pdb_name, pdb_file);
			char chain = *(strrchr(pdb_name.c_str(),'_')+1);
			fprintf(template_F, "%s %s %c\n", aln.name(i).c_str(), pdb_file.c_str(), chain );
		}
		remove(tmp_out.c_str());
	}
	fclose(template_F);


	EBI_Rest::lib_stop();
}



// 1AAB.pdb