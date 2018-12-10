#include <candidate.h>
#include <algorithm> // For remove_if

/*Constructor saves the truth value to instance variables,
 * Could be exchanged for seed values, as required.*/
Candidate::Candidate(const double& XT, const double& YT, const double& ET) : _XT(XT), _YT(YT), _ET(ET) {}

/*Method imports event 'details' from 3 vectors*/
void Candidate::importDetails( const std::vector<double>& X, const std::vector<double>& Y, 
		const std::vector<double>& P) {
	/*Import TC readouts by iterating over X,Y,P*/
	entry tmp;
	auto xi = X.begin();
	auto yi = Y.begin();
	auto pi = P.begin();
	/*Should add 'and' conditions for each vector, however currently assuming all saeme length*/	
	while (xi != X.end()) {
		tmp.x = *xi; tmp.y = *yi; tmp.p = *pi; 
		_Entries.push_back(tmp);
		/*Increment iterators*/
		xi++; yi++; pi++;
	}
}

struct inCircle {
	/* Declaration + Definition for the inCricle unary operator*/
    double _R, _XT,_YT;    	
	/* Constructor */
	inCircle(const double& R, const double& XT, const double& YT): _R(R), _XT(XT), _YT(YT) {}
	/* Unary operator */
	bool operator()(const entry& E) const {
		// Comparision here, (true if outside R!)
		return ( ((E.x-_XT)*(E.x-_XT) + (E.y-_YT)*(E.y-_YT)) > _R*_R ); 
	}
}; 

void Candidate::crop(const double& R) {
	/*Removes entries in the _Entries if outside circle!*/
	_R = R;	
	// Delete entries which are outside circle
	//std::cout << _Entries.size() <<"\t";
	_Entries.erase( 
			std::remove_if ( _Entries.begin(), _Entries.end(), inCircle(_R, _XT, _YT) ),
			_Entries.end() 
		);
	//std::cout << _Entries.size() << "\n";
} 


void Candidate::calculate_resolutions() {
	/*Calculate the e_res, x_res, y_res, e_sum*/
	e_res =0; x_res=0; y_res=0; r_res=0; e_sum=0;	
	// Loop over all candidate entries
	for(auto const& e: _Entries) {// weighted sums 
		e_sum += e.p; x_res += e.p*e.x; y_res += e.p*e.y; 
	}
	// Divide by weights to obtain weighted average
	x_res /= e_sum; y_res /= e_sum;
    // Energy resolution -> Sum(E) - E_truth	
	e_res = e_sum-_ET; 
	// Subtract the TRUTH values to obtain resolution for single event
	x_res -= _XT; y_res -= _YT; 
	// Calcualte the radial resolution
	r_res = std::sqrt (x_res*x_res + y_res*y_res) ;
}

/*GET-er methods below*/
double Candidate::getERes() { return e_res;} 
double Candidate::getXRes() { return x_res;}
double Candidate::getYRes() { return y_res;}
double Candidate::getRRes() { return r_res;}
double Candidate::getESum() { return e_sum;} 

