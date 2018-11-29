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

public:
	Candidate(const double& XT, const double& YT, const double& ET);
	void importDetails(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& P);
	//bool inCircle(const entry& E);
	void crop(const double& R); 
};

#endif
