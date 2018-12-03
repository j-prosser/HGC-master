#ifndef _candidate_h_
#define _candidate_h_

#include <cstring>
#include <time.h>
#include <cmath>
#include <vector>


typedef struct Entry {
	double x,y,p; 
} entry; 


class Candidate {
private:
	// Data
	std::vector<entry> _Entries;
	// Truth Data
	double _XT, _YT, _ET;
	double _R;

	/*Filled by calculations*/
	double e_res;
	double x_res;
	double y_res;
	double r_res;
	double e_sum;	

public:
	Candidate(const double& XT, const double& YT, const double& ET);
	void importDetails(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& P);
	//bool inCircle(const entry& E);
	void crop(const double& R); 

	// Sets variables, 
	void calculate_resolutions();
	

	//populate/update the below doubles.
	/*Maybe implement get_ methods?*/
	double getERes();
	double getYRes();
	double getXRes();
	double getRRes();
	double getESum();

};

std::vector<double> generate_R (const double& start, const double& stop, const double& inc);

double Deviation(std::vector<double>& V, double& mean);

double Average(std::vector<double>& V);


#endif
