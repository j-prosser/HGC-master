#ifndef _functions_
#define _functions_
// functions
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <cstring>
#include <time.h>
#include <TLorentzVector.h>
#include <TMultiGraph.h>
#include <TMatrixD.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEllipse.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TKey.h>
#include <TLatex.h>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>

#include <TGraph.h>
using namespace std; 

void plotFB(const map<string,vector<double>> &DV, const vector<string>  &sV,string fname,  map<string, TGraph*> &outmap);

#endif
