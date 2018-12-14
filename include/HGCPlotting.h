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

#include "VPMC.h" 


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
  
	/*Command line HGCPlotting instance variables*/
	int _verbose;

	int _max_events;			// -1 default : all events.
	int _energy_resolution;		// 1 : on (default), 0 : off
	int _position_resolution;	// 1 : on, 0 : off, NOT IMPLEMENTED
	int _PU;					// pile up, default 0, NOT IMPLEMENTED
	int _view_event;			// View/debug a view of a single event 
	std::string _compareC3D;			// compare CMSSW ntuples with HGC-CMSSW-ANLYSIS

	/*Ploting Structures*/
	std::vector<std::string> _HistoSets;
	TChain * _chain    ;  

	typedef std::map< std::string, TH1D* > histmap;
	typedef std::map< std::string, histmap > clonemap;
	typedef std::map< std::string, double > doublemap;

	clonemap _cloned_hists;  // map string -> histmap -> string-> TH1D
	doublemap _event_variables; // map string -> double

	typedef std::map<std::string, std::vector<double>> vectormap;
	typedef std::map<double, std::vector<double>> doublevectormap;
	/* string : double : vector(doubles) 
	* data_id : r_curr : vector(data) */
	stringdoublemap		_radial_reconstruction;
	vectormap			_event_details; // map of vectors 
	stringdoubledouble	_radial_results;

	// Datastructures for eta-separated radial reconstruction...
	typedef std::map<std::string, stringdoublemap>		stringstringdoublemap;
	typedef std::map<std::string, stringdoubledouble>	stringstringdoubledouble;
	stringstringdoublemap		_radial_eta_reconstruction;	
	stringstringdoubledouble	_radial_eta_results;		

	// Pointer Maps {std::string : TPointer}
	typedef std::map<std::string, TH2D*> plot_2d_map;
	typedef std::map<std::string, TGraph*> map_tgraphs;	
	plot_2d_map _2d_plots; 
	map_tgraphs _graphs; 

	VPMC * _vcomp; // HCA
	typedef std::map<std::string, std::vector<std::vector<double>> > map2dvector;
	map2dvector dcl;

	public :
  HGCPlotting( CmdLine * cmd );
  ~HGCPlotting();
  
  void DoNothing();

  void SetupRoot();

  void SetupVPMC();

  //void SetupFillSingleEvent(); 

  void LoadHistoTemplates( std::string name );

  void MakeAllHists( std::vector<std::string> &HistoSets);

  void Loop();

  void Fill();

  void CalculateTriggerCellVariables();

  //void CalculateReducedCircle(const double& R);

  void PlotEnergyResolution();

  void CalculateCircleStats();
  
  void GraphReducedCircle(stringdoublemap& dataset, std::string graph_name);

  void FillAllHists( std::string name );

  //void FillPlot(std::string name); // potentially use this?

  bool FileExists( std::string file );

  void Save();
 
};

#endif
