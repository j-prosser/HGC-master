/********************************************************************************
HCG plotting code
Samuel Webb
Imperial College
********************************************************************************/

#ifndef _HGCPlotting_h_
#define _HGCPlotting_h_
#define XXX std::cout << "I am here: " << __LINE__ << std::endl << std::endl;

#include "BuildTreeBase.h"
#include "tdrstyle.h"
#include "CmdLine.h"
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


typedef std::map< std::string, std::map<double, double> > stringdoubledouble;
typedef std::map<std::string, std::map<double,std::vector<double>>> stringdoublemap;

class HGCPlotting : public BuildTreeBase {
	public: 
		
	private :
  
  /***** General *****/
  CmdLine * _cmd ;
  std::string _in_directory;
  std::string _out_directory;
  TDirectory * _origDir ;
  
  /*Globals in HGCPlotting*/
  int _max_events;
  std::vector<std::string> _HistoSets;
  TChain * _chain    ;  

  typedef std::map< std::string, TH1D* > histmap;
  typedef std::map< std::string, histmap > clonemap;
  typedef std::map< std::string, double > doublemap;


  clonemap _cloned_hists;  // map string -> histmap -> string-> TH1D
  doublemap _event_variables; // map string -> double



	// define vector of 
	typedef std::map<std::string, TH2D*> plot_2d_map;
	typedef std::map<std::string, TGraph*> map_tgraphs;


  //doublemap _event_single;
	typedef std::map<std::string, std::vector<double>> vectormap;
  
	typedef std::map<double, std::vector<double>> doublevectormap;

  /* string : double : vector(doubles) 
   * data_id : r_curr : vector(data) */
  stringdoublemap _radial_reconstruction;
  
  vectormap _event_details; // map of vectors 

  stringdoubledouble _radial_results;

	// Datastructures for eta-separated radial reconstruction...
	typedef std::map<std::string, stringdoublemap> stringstringdoublemap;
	typedef std::map<std::string, stringdoubledouble> stringstringdoubledouble;
	stringstringdoublemap _radial_eta_reconstruction; // NEW
	stringstringdoubledouble _radial_eta_results; // NEW


  plot_2d_map _2d_plots; 
  map_tgraphs _graphs; 

 public :
  HGCPlotting( CmdLine * cmd );
  ~HGCPlotting();
  
  void DoNothing();

  void SetupRoot();

  //void SetupFillSingleEvent(); 

  void LoadHistoTemplates( std::string name );

  void MakeAllHists( std::vector<std::string> &HistoSets);

  void Loop();

  void Fill();

  void CalculateTriggerCellVariables();

  //void CalculateReducedCircle(const double& R);

  void CalculateCircleStats();

  
  void GraphReducedCircle(stringdoublemap& dataset, std::string graph_name);

  void FillAllHists( std::string name );

  //void FillPlot(std::string name); // potentially use this?

  bool FileExists( std::string file );

  void Save();
  
};


#endif
