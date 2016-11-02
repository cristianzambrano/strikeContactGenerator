#ifndef RESTSERVICES_H_
#define RESTSERVICES_H_


#include <string>

#include <curl/curl.h>



class EBI_Rest
{

private:

	static std::string _email;

public:

	// Static Methods
	void static lib_start(const std::string &email);

	void static lib_stop()
	{
		curl_global_cleanup();
	}

	// Methods
	void static submit_blast_job(const std::string &program, const std::string &type, const std::string &database, double ecut, const std::string &name, const std::string &seq, std::string &job_id);

	short static get_job_status(const std::string &job_id);

	void static get_job_result(const std::string &job_id, const std::string &out_f);
};


class PDB_Rest
{
private:

public:
	void static get_pdb_file(const std::string &pdb_id, const std::string &out_f);
};





#endif // RESTSERVICES_H_