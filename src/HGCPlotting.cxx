/******************************************************************************
HGC analysis code
Samuel Webb
Imperial College
***************************************************************************/

// Time-stamp: <2018-11-01 13:44:34 (snwebb)>

#include "HGCPlotting.h"
#include "algorithm"
#include "functions.h"

HGCPlotting::HGCPlotting( CmdLine * cmd ){

  _cmd = cmd; 
  _origDir = gDirectory ;
  _in_directory = _cmd->string_val( "--in_directory" ) ;
  _out_directory = _cmd->string_val( "--out_directory" ) ;
  _max_events =  _cmd->int_val( "--max_events"  , -1);


  // ANYTHING ADDED TO _HistoSets will create a histogram, if defined in src/Histogramming.cxx, 
  // where "NAME" is the argument checked for in Histogramming

//  _HistoSets.push_back( "PU0_General" );
  _HistoSets.push_back( "PU0_forward" );
  _HistoSets.push_back( "PU0_backward" );

  //yoyo
//  _HistoSets.push_back( "TriggerCells" );
  //_HistoSets.push_back( "" );

}


HGCPlotting::~HGCPlotting() {

  gDirectory = _origDir ;

  for(auto &it1 : _cloned_hists) {
    for(auto &it2 : it1.second) {
      it2.second->Delete();
    }
  }

}

void HGCPlotting::DoNothing(){
	std::cout << "DoNothing" <<std::endl;
}


void HGCPlotting::SetupRoot(){


  _chain = new TChain ("hgcalTriggerNtuplizer/HGCalTriggerNtuple");

  //  std::string remotedir = "root://cms-xrd-global.cern.ch//store/user/sawebb/SingleGammaPt25Eta1p6_2p8/crab_SingleGammaPt25_PU0-stc/181025_100629/0000/";

  for ( int i = 1; i < 2; i++ ){
        if (FileExists( (_in_directory + "/ntuples/ntuple_" + std::to_string(i) + ".root"   ).c_str()  )  )   
          _chain  ->Add ( (_in_directory + "/ntuples/ntuple_" + std::to_string(i) + ".root"   ).c_str() );
    //  _chain  ->Add ( (_in_directory + "/ntuple_" + std::to_string(i) + ".root"   ).c_str() );
  }
  // Finally make them
  MakeAllHists( _HistoSets ); 
}


void HGCPlotting::Fill(){

  Init( _chain );
  Loop( );

}


void HGCPlotting::Loop( ){

  //Loop over events

  std::cout << "Beginning Event Loop" << std::endl;

  if (fChain == 0) return;

  std::cout << "Getting Number of Entries" << std::endl;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Entries: " << nentries << std::endl;
  
  Long64_t nbytes = 0, nb = 0;

  /////////////////////////////////
  // SET BRANCHES TO BE 'PROCESSED'
  // 	unset all branches
  fChain->SetBranchStatus("*",0);
  //	set specificed/needed branches
  fChain->SetBranchStatus("tc_n",1);
  fChain->SetBranchStatus("tc_pt",1);
  fChain->SetBranchStatus("tc_phi",1);
  fChain->SetBranchStatus("tc_eta",1);
  //TRUTH DATA
  fChain->SetBranchStatus("gen_phi",1);
  fChain->SetBranchStatus("gen_eta",1);  
  fChain->SetBranchStatus("gen_energy",1); 
  fChain->SetBranchStatus("gen_pt",1);

  //Detector Geometry- DOESNT WORK :(
  
  fChain->SetBranchStatus("tc_x",1);
  fChain->SetBranchStatus("tc_y",1);
  fChain->SetBranchStatus("tc_z",1);
  fChain->SetBranchStatus("tc_layer",1);
  fChain->SetBranchStatus("tc_zside",1);
         
  //fChain->SetBranchStatus("",1);
  // MAIN for  LOOP, FILLS HISTOGRAMS
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // yoyo change 10000 devisor to 100
    if ( jentry % 10000 == 0 ) {
      std::cout << jentry << "/" << nentries  << "\n" ;
    }
    // DISPLAY TRUTH INFOMATION <- Move this somewhere
    if (1 == 0) {
    // Every jentry has different truth values!    
    std::cout << "***TRUTH DATA***"  << std::endl;
    std::cout << "phi: "<< gen_phi->at(0)<<"/"<<gen_phi->at(1) <<std::endl;
    std::cout << "eta: "<< gen_eta->at(0)<<"/"<<gen_eta->at(1) <<std::endl;
    std::cout << "energy: "<< gen_energy->at(0)<<"/"<<gen_energy->at(1)<<std::endl;
    // NOT TESTED std::cout << "pt:" << gen_pt->at(0)<<"/"<<gen_pt->at(1)<<std::endl;
    //std::cout << "SIZE TRUTH: "<< gen_eta->size()<< std::endl;
    }
    if ( _max_events != -1 ){
      if ( jentry > _max_events ) break;
    }
    for (auto& names : _HistoSets ){
	  //std::cout << "TESTING1 "<< names << std::endl;
      FillAllHists( names );
    }
	if (jentry==0) {
		/*Create plot of single event*/
	    // Create 2D histogram of event 0	
		
		
		//int xs = _event_details["xnf"].size();
		//int ys = _event_details["ynf"].size();
		
		//auto result = std::minmax_element (foo.begin(),foo.end());
	   
		//std::vector<double> test = {1,2,3,2,-1};
	    //auto lim = std::minmax_element(test.begin(), test.end());
		//std::cout << *lim.first << " here" << std::endl;	
	
		std::cout << "PLOTTING EVENT 0" << std::endl;		

//		auto x_lims = std::minmax_element (_event_details["xnf"].begin(), _event_details["xnf"].end());	
//		auto y_lims = std::minmax_element (_event_details["ynf"].begin(), _event_details["ynf"].end());
        
        std::cout << "size of entire event " << _event_details["xnf"].size() <<std::endl; 		
		
        //std::cout << "x lims: " << *x_lims.first << " "<< *x_lims.second << std::endl; 

//	    _2d_plots["test_plot"] = new TH2D ("normalised_coord", "",500,*x_lims.first,*x_lims.second,100,*y_lims.first,*y_lims.second);
		//_2d_plots["test_plot"]->Fill(_event_details["xnf"],_event_details["ynf"],_event_details["ptf"]);
		//_2d_plots["test_plot"]->Draw("COLZ");	
//		for (unsigned i =0; i < _event_details["ptf"].size(); ++i) {
//			_2d_plots["test_plot"]->Fill(_event_details["xnf"][i], _event_details["ynf"][i], _event_details["ptf"][i]);	
//		}
		//_2d_plots["test_plot"]->Draw("COLZ");
		// other plot
//		auto xlc = std::minmax_element(_event_details["xnfc"].begin(),_event_details["xnfc"].end());
//        auto ylc = std::minmax_element(_event_details["ynfc"].begin(),_event_details["ynfc"].end());
//		_2d_plots["cand_plot"] = new TH2D ("normalised_coord_cand_cirlce_0.05", "", 100,*xlc.first,*xlc.second,100,*ylc.first,*ylc.second);
		
		std::cout <<"size of event with circle "<< _event_details["xnfc"].size()<< std::endl;
		
//		_2d_plots["cand_plot"] = new TH2D ("normalised_coord_cand_cirlce_0.05", "", 100,-1,1,100,-1,1);
//		for (unsigned i=0; i<_event_details["ptfc"].size(); ++i) {
//			_2d_plots["cand_plot"]->Fill(_event_details["xnfc"][i],_event_details["ynfc"][i]);
//		}	
	}
  }

}


bool HGCPlotting::FileExists( std::string file ){

  struct stat buf;
  if (  stat (  ( file ).c_str(), &buf ) == 0)
    return true;
  else{
    return false;
  }
}

void HGCPlotting::Save(){

  // This needs rework!!!
  std::system( ("mkdir -p output/" + _out_directory  )   .c_str() );  
  TFile * f_hists = new TFile(  ("output/" + _out_directory + "/output.root").c_str(), "RECREATE" );
 
  for(auto &it1 : _cloned_hists) {
    for(auto &it2 : it1.second) {
      it2.second->Write();
    }
  }
  for(auto &it1 : _2d_plots){
	it1.second->Write();
  } 

  f_hists->Close();
}


//////////





