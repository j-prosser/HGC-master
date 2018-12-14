#include <functions.h>
//#include <string>
//#include <vector>
//#include <map> 
//#include <TGraph>
using namespace std;

void plotFB(const map<string,vector<double>> &DV, const vector<string>  &sV,string fname,  map<string, TGraph*> &outmap){
	string ffname = "forward_"	+	fname; 
	string bfname = "backward_" +	fname;

	//f
	outmap[ffname] = new TGraph(DV.at(sV.at(0)).size(), &DV.at(sV.at(0)).at(0), &DV.at(sV.at(1)).at(0) );
	outmap[ffname]->SetName(ffname.c_str());
	outmap[ffname]->SetMarkerStyle(37);
	//b
	outmap[bfname] = new TGraph(DV.at(sV.at(2)).size(), &DV.at(sV.at(2)).at(0), &DV.at(sV.at(3)).at(0));
	outmap[bfname]->SetName(bfname.c_str());
	outmap[bfname]->SetMarkerStyle(37);
}

void CalcXY() {} 
