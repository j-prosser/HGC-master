
#include "HGCPlotting.h"
#include <cmath> 

#include <candidate.h>
#include <TVectorD.h> 


void HGCPlotting::MakeAllHists( std::vector<std::string> &HistoSets){
	/*Runs LoadHistoTemplates for name in Histosets*/
	std::cout << "Creating All Histograms: " << std::endl;
	for (auto& names : HistoSets ){
            LoadHistoTemplates (names);
			std::cout << "\t" << names << "\n";
    }  
}
// TH1D Syntax: new TH1D(name, title, nbinsx, xlow, xhigh) 
//		xlow-> edge of lowest bin, xhigh -> edge of highest bin
// or: new TH1D(name, title, nbinsx, xbins) ; xbins -> low edge of of everybin, contains nbinsx+1 entries. 

void HGCPlotting::LoadHistoTemplates( std::string name ) { 
	/* Initialise all empty histograms, for a given name */	
	if ( name == "TriggerCells" ) {
		_cloned_hists[ name ] [ "tc_eta" ] = new TH1D ( (name + "_tc_eta").c_str(), "", 100,-5,5 );  
		_cloned_hists[ name ] [ "tc_phi" ] = new TH1D ( (name + "_tc_phi").c_str(), "", 100,-M_PI,M_PI );  // M_PI == Pi 
	} else if ( name == "PU0_General" ){
		  _cloned_hists[ name ] [ "tc_n" ] = new TH1D ( (name + "_tc_n").c_str(), "", 100,0,1000 );  
	} else if ( name == "PU200" ){		
		_cloned_hists[ name ] [ "tc_n" ] = new TH1D ( "default_tc_n", "", 100,0,100000 );  
		_cloned_hists[ name ] [ "ex_sum" ] = new TH1D ( "default_ex_sum", "", 150,-150,150);
		_cloned_hists[ name ] [ "ey_sum" ] = new TH1D ( "default_ey_sum", "", 150,-150,150);
		_cloned_hists[ name ] [ "er_sum" ] = new TH1D ( "default_er_sum", "", 150,0,150);
		_cloned_hists[ name ] [ "ephi_sum" ] = new TH1D ( "default_ephi_sum", "", 150,-M_PI,M_PI);
	} else if ( name == "PU0_forward" || name == "PU0_backward" ){
		_cloned_hists[ name ] [ "tc_n" ] = new TH1D ( (name + "_tc_n").c_str(), "", 100,0,1000 );  
		_cloned_hists[ name ] [ "ex_sum" ] = new TH1D ( (name + "_ex_sum").c_str(), "", 150,-50,50);
		_cloned_hists[ name ] [ "ey_sum" ] = new TH1D ( (name + "_ey_sum").c_str(), "", 150,-50,50);
		_cloned_hists[ name ] [ "er_sum" ] = new TH1D ( (name + "_er_sum").c_str(), "", 150,0,30); // CHANGED 2->30, edge of last bin...
		_cloned_hists[ name ] [ "ephi_sum" ] = new TH1D ( (name + "_ephi_sum").c_str(), "", 150,-M_PI,M_PI);// energy sum over all phi (-pi to pi)
		_cloned_hists[ name ] [ "dphi_met" ] = new TH1D ( (name + "_dphi_met").c_str(), "", 150,-0.05,0.05);// delta Phi, difference between phi
		_cloned_hists[ name ] [ "denergy" ] = new TH1D ( (name + "_denergy").c_str(), "", 150, -5.0,5.0); // difference in energy
		//_cloned_hists[ name ] [ "dpos" ] = new TH1D ( (name + "_dpos").c_str(), "", 150, 0,.1);
		//_cloned_hists[ name ] [ "dpos_2" ] = new TH1D ( (name + "_dpos_E").c_str(), "", 150, 0,.1);
		//_cloned_hists[ name ] [ "dpos_X" ] = new TH1D ( (name + "_dpos_x").c_str(), "", 150, -.04,.04);
		//_cloned_hists[ name ] [ "dpos_Y" ] = new TH1D ( (name + "_dpos_y").c_str(), "", 150, -.04,.04);
		//_cloned_hists[ name ] [ "dpos_X_E" ] = new TH1D ( (name + "_dpos_x_E").c_str(), "", 150, -.04,.04);
		//_cloned_hists[ name ] [ "dpos_Y_E" ] = new TH1D ( (name + "_dpos_y_E").c_str(), "", 150, -.04,.04);
		_cloned_hists[ name] [ "denergy_R" ] = new TH1D ( (name+"_denergy_R").c_str(), "", 150, -25.,25.);											   
	} else if ( name == "Radial_Reconstruction" )  { 
	//_graphs["sig_e_e_r"] = new TGraph(); // n,x,y	
	// NO NEED TO INITALISE GRAPHS HERE!
	}
}
void HGCPlotting::CalculateTriggerCellVariables() {
	/*Calculates the Basic trigger readouts, without doing intensive calculations!*/
	
	//std::cout << tc_zside->size() << " " << tc_layer->size() <<" "<<tc_z->size()<<std::endl;
    /* ENERGY RESOLUTION */
	//Ex and Ey sums
	double exsum_forward = 0;
	double eysum_forward = 0;
	double exsum_backward = 0;
	double eysum_backward = 0;
  
	// intialise _Event_details by deleting, for every loop!
	_event_details.clear();
	
	// LOOP OVER ALL ENTRIES
	for (unsigned int i = 0; i < tc_pt->size(); i++){
		// IF-statement to filter between forward and backward calorimeters...
		if  ( tc_eta->at(i)>0){
		//FORWARD
		exsum_forward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_forward += tc_pt->at(i)*std::sin(tc_phi->at(i));
		//
		_event_details["fX"].push_back( tc_x->at(i) / tc_z->at(i) ); 
		_event_details["fY"].push_back( tc_y->at(i) / tc_z->at(i) );
		_event_details["fP"].push_back( tc_pt->at(i) ); 
		} else {
		//BACKWARD
		exsum_backward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_backward += tc_pt->at(i)*std::sin(tc_phi->at(i));
		//
		_event_details["bX"].push_back(- tc_x->at(i) / tc_z->at(i) );
		_event_details["bY"].push_back(- tc_y->at(i) / tc_z->at(i) );
		_event_details["bP"].push_back( tc_pt->at(i) );
		}
	}
	/* Total energy sums in r, forward & backward*/
	double ersum_forward = std::sqrt( exsum_forward*exsum_forward + eysum_forward*eysum_forward );
	double ersum_backward = std::sqrt( exsum_backward*exsum_backward + eysum_backward*eysum_backward );
	  
	/* Total energy sum in phi, forwar & backward*/ 
	// std::atan2 -> arctan(y/x), uses signs of y,x to find the correct quadrant!
	double ephisum_forward = std::atan2( eysum_forward, exsum_forward ); 
    double ephisum_backward = std::atan2( eysum_backward, eysum_backward );  

	/* Save calculated variables to be used in histograms */
	_event_variables[  "ex_sum_forward"  ] = exsum_forward;
	_event_variables[  "ey_sum_forward"  ] = eysum_forward;
	_event_variables[  "er_sum_forward"  ] = ersum_forward;
	_event_variables[  "ephi_sum_forward"  ] = ephisum_forward;

	_event_variables[  "ex_sum_backward"  ] = exsum_backward;
	_event_variables[  "ey_sum_backward"  ] = eysum_backward;
	_event_variables[  "er_sum_backward"  ] = ersum_backward;
	_event_variables[  "ephi_sum_backward"  ] = ephisum_backward;

	/* Naive Calculations between the TRUTH and the measured values to find the difference in measured phi
	 * and truth phi */
	// Declare temp variables
	TVector2 met, truth; // init measured, truth 2-d (x-y) vectors
	double dphi_met; 
	// Forward Case
	met.Set( exsum_forward, eysum_forward ); 
	truth.SetMagPhi( 1, gen_phi->at(0) ); // gen_phi is the TRUTH VALUE for phi.
	dphi_met = met.DeltaPhi( truth );
	_event_variables[  "dphi_met_forward"  ] =  dphi_met;
	// Backward Case
	met.Set( exsum_backward, eysum_backward );
	truth.SetMagPhi( 1, gen_phi->at(1) );
	dphi_met = met.DeltaPhi( truth );
	_event_variables[  "dphi_met_backward"  ] =  dphi_met;
	/* Naive calculation for difference in measured energy and truth energy */
	_event_variables["denergy_forward"] = ersum_forward - gen_pt->at(0);
	_event_variables["denergy_backward"] = ersum_backward - gen_pt->at(1);

	/*FIND TRUTH VALUES*/ 
	// in our normalised co-ordinates, Y = y/z, X = x/z
	_event_variables["xnft"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(0));
	_event_variables["ynft"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(0));
	_event_variables["xnbt"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(1)); 
	_event_variables["ynbt"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(1));  
}

/* Mean and Standard Deviation calculations for a vector of double,
 * (could implement a template such that any data_type coulde be used) */
double Mean(std::vector<double>& V) {
	/*Finds the mean (average) of a vector of doubles*/	
	double sum = std::accumulate(V.begin(), V.end(), 0.0);
	return sum / static_cast<double>(V.size());
}
double Deviation(std::vector<double>& V, double& mean){
	/*Finds the standard deviation of a vector given its mean.*/
	double sq_sum = 0;
	for(unsigned n = 0; n < V.size() ; n++ )
	{
		sq_sum += (V[n] - mean) * (V[n] - mean);
	}
	sq_sum /= static_cast<double>(V.size());
	return std::sqrt(sq_sum);
} 

std::vector<double> generate_R (const double& start, const double& stop, const double& inc) {
	/*Generate a vector of decreasing R, returns as a vector of edges*/
	std::vector<double> tmp;
	for (double r_tmp=start; r_tmp > stop; r_tmp -= inc ) { tmp.push_back(r_tmp); }
	return tmp; 
}

void HGCPlotting::FillAllHists( std::string name ){
	/* Method to Fill histograms and perform other functions for specified 'name' in _HistoSets */
	/* RUN FOR EVERY __name__ IN EVERY EVENT */	
	// syntax for filling an ininialised histogram below.:
		//_cloned_hist [ name ] [ HIST NAME ] ->Fill ([_event_variables["var_name"]]); 
	if ( name == "TriggerCells" ){
		// Visualise TCs...
		for (unsigned int i = 0; i < tc_eta->size(); i++){
			_cloned_hists[ name ] [ "tc_eta" ] ->Fill ( tc_eta->at(i) );
			_cloned_hists[ name ] [ "tc_phi" ] ->Fill ( tc_phi->at(i) );
		}
	} else if ( name == "PU0_General" ){
    	_cloned_hists[ name ] [ "tc_n" ] ->Fill ( tc_n );
	} else if ( name == "PU0_forward" ){
		// Histograms for entire event, sums over all TCs
		_cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_forward"  ] );
	    _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_forward"  ] );
		_cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_forward"  ] );
	    _cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_forward"  ] );        
		_cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_forward"  ] );    
		_cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_forward"] );        
	} else if ( name == "PU0_backward" ){
		_cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_backward"  ] );
	    _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_backward"  ] );
	    _cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_backward"  ] );
		_cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_backward"  ] );   
		_cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_backward"  ] );                
		_cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_backward"] );
	} else if ( name == "Radial_Reconstruction") {
		/*Routine to calculate each event*/		
	
		// USE: _radial_reconstruction INSTEAD, <std::string, vector<double> >
		/*Generate a vector of decreasing R (for each event)*/
		// Define range of R's to be calculated

		//Move this to a different scope?
		std::vector<double> Rs = generate_R(0.1, 0.005, 0.002); 
		//std::vector<double> ETAs = generate_R(3.0, 1.47, 0.51);

		// Initialise with TRUTH values, could use seeds here? 
		Candidate fCand( _event_variables["xnft"], _event_variables["ynft"], gen_pt->at(0)); 
		Candidate bCand( _event_variables["xnbt"], _event_variables["ynbt"], gen_pt->at(1));

		// Import event details (readout)
		fCand.importDetails(_event_details["fX"], _event_details["fY"], _event_details["fP"]);
		bCand.importDetails(_event_details["bX"], _event_details["bY"], _event_details["bP"]);
		
		unsigned r_idx = 0;
		for (auto& r_curr : Rs) {	
			// Crop everything outside r
			fCand.crop(r_curr); 
			bCand.crop(r_curr);
			// Do calculations
			fCand.calculate_resolutions();
			bCand.calculate_resolutions();
			// Now we can get the data that we want 
			if (0) {
				//DEBUG	
				std::cout << r_curr<<"\n";	  
				std::cout << "fERes/ESum" <<"\t\t" <<fCand.getERes() << "\\" << fCand.getESum() <<"\n"; 
				std::cout << "fPosition Res.\t" << fCand.getXRes() << "\n";
				std::cout << "bERes/ESum" <<"\t\t" <<bCand.getERes() << "\\" << bCand.getESum() <<"\n"; 
				std::cout << "bPosition Res.\t" << bCand.getXRes() << "\n";
			}
			// Store results for each R, in _radial_reconstruction
			_radial_reconstruction["e_res"][r_curr].push_back( fCand.getERes() );
			_radial_reconstruction["e_res"][r_curr].push_back( bCand.getERes() );
			_radial_reconstruction["e_sum"][r_curr].push_back( fCand.getESum() );
			_radial_reconstruction["e_sum"][r_curr].push_back( bCand.getESum() );
             

			//Separate by ETA, assume; gen_eta_forward == -gen_eta_backward 
			if (1.47<= gen_eta->at(0) &&  gen_eta->at(0) < 1.98) {
				_radial_eta_reconstruction["1.47_1.98"]["e_res"][r_curr].push_back(fCand.getERes());
				_radial_eta_reconstruction["1.47_1.98"]["e_res"][r_curr].push_back(bCand.getERes());		
				_radial_eta_reconstruction["1.47_1.98"]["e_sum"][r_curr].push_back(fCand.getESum());
				_radial_eta_reconstruction["1.47_1.98"]["e_sum"][r_curr].push_back(bCand.getESum());
			} else if (1.98 <= gen_eta->at(0) && gen_eta->at(0) < 2.49) {
			
				_radial_eta_reconstruction["1.98_2.49"]["e_res"][r_curr].push_back(fCand.getERes());
				_radial_eta_reconstruction["1.98_2.49"]["e_res"][r_curr].push_back(bCand.getERes());		
				_radial_eta_reconstruction["1.98_2.49"]["e_sum"][r_curr].push_back(fCand.getESum());
				_radial_eta_reconstruction["1.98_2.49"]["e_sum"][r_curr].push_back(bCand.getESum());
			} else if (2.49 <= gen_eta->at(0) && gen_eta->at(0) < 3.00) {
			
				_radial_eta_reconstruction["1.98_3.00"]["e_res"][r_curr].push_back(fCand.getERes());
				_radial_eta_reconstruction["1.98_3.00"]["e_res"][r_curr].push_back(bCand.getERes());		
				_radial_eta_reconstruction["1.98_3.00"]["e_sum"][r_curr].push_back(fCand.getESum());
				_radial_eta_reconstruction["1.98_3.00"]["e_sum"][r_curr].push_back(bCand.getESum());
			}

			// Increment loop (important!)
			r_idx +=1;
		}
		
	}
}

void HGCPlotting::CalculateCircleStats(  ) {
	/*Does stuff on _radial_reconstruction dataset*/
	/*Output onto _radial_results*/
	
	for (auto& r_pair : _radial_reconstruction["e_sum"]) {
		_radial_results["mean_e_sum"][r_pair.first] = Mean(r_pair.second);
	}

	/*TEMP vectors for plotting*/
	std::vector<double> tmp_r;
	std::vector<double> tmp_sigee; 
	for (auto & r_pair : _radial_reconstruction["e_res"]) {
		/*Find the mean, */ 
		_radial_results["mean_e_res"][r_pair.first] = Mean(r_pair.second); 
		_radial_results["stdev_e_res"][r_pair.first] = Deviation(r_pair.second, _radial_results["mean_e_res"][r_pair.first] );
		//std::cout << r_pair.first << "\t" << "Mean e_res " << Mean(r_pair.second) << "\n";
		_radial_results["sig_e_e"][r_pair.first] = _radial_results["stdev_e_res"][r_pair.first] / _radial_results["mean_e_sum"][r_pair.first]; 
		/*Print results?*/
		//std::cout << r_pair.first << "\t" << "sigma_E/E\t" << _radial_results["sig_e_e"][r_pair.first] << "\n"; 	
		tmp_r.push_back(r_pair.first);
		tmp_sigee.push_back(_radial_results["sig_e_e"][r_pair.first]); 
	}

	/*Plottign Stuff -> Move somewhere else!?*/ 
	//Make some plots! 
	// TVector<float> tvf(svf.size(), &svf[0]);
	/*Convert from std::vector<double> to TVectorD*/
	TVectorD tv_r(tmp_r.size(), &tmp_r[0]);	
	TVectorD tv_sigee(tmp_sigee.size(), &tmp_sigee[0]);

	/*Add graph to map of graphs created*/
	_graphs["sigma_EE_r"] = new TGraph( tv_r,tv_sigee );  
	_graphs["sigma_EE_r"]->SetName("Sigma_EE_r"); //set name of plot!
	/*Perhaps add other stuff like a title*/
}


void HGCPlotting::GraphReducedCircle(stringdoublemap& dataset, std::string graph_name){
	/*Does stuff on dataset dataset*/
	/*Output onto tmp_results*/
	stringdoubledouble tmp_results;

	for (auto& r_pair : dataset["e_sum"]) {
		tmp_results["mean_e_sum"][r_pair.first] = Mean(r_pair.second);
	}

	/*TEMP vectors for plotting*/
	std::vector<double> tmp_r;
	std::vector<double> tmp_sigee; 
	for (auto & r_pair : dataset["e_res"]) {
		/*Find the mean, */ 
		tmp_results["mean_e_res"][r_pair.first] = Mean(r_pair.second); 
		tmp_results["stdev_e_res"][r_pair.first] = Deviation(r_pair.second, tmp_results["mean_e_res"][r_pair.first] );
		//std::cout << r_pair.first << "\t" << "Mean e_res " << Mean(r_pair.second) << "\n";
		tmp_results["sig_e_e"][r_pair.first] = tmp_results["stdev_e_res"][r_pair.first] / tmp_results["mean_e_sum"][r_pair.first]; 
		/*Print results?*/
		//std::cout << r_pair.first << "\t" << "sigma_E/E\t" << tmp_results["sig_e_e"][r_pair.first] << "\n"; 	
		tmp_r.push_back(r_pair.first);
		tmp_sigee.push_back(tmp_results["sig_e_e"][r_pair.first]); 
	}

	/*Plottign Stuff -> Move somewhere else!?*/ 
	//Make some plots! 
	// TVector<float> tvf(svf.size(), &svf[0]);
	/*Convert from std::vector<double> to TVectorD*/
	TVectorD tv_r(tmp_r.size(), &tmp_r[0]);	
	TVectorD tv_sigee(tmp_sigee.size(), &tmp_sigee[0]);

	/*Add graph to map of graphs created*/
	graph_name = "sigma_EE_r_"+graph_name; 
	std::cout << "Graphing:\t" << graph_name <<"\n";
	_graphs[graph_name] = new TGraph( tv_r,tv_sigee );  
	_graphs[graph_name]->SetName(graph_name.c_str()); //set name of plot!
	/*Perhaps add other stuff like a title*/
}


/* TODO: 
 *	- Implement a way / variable to seed
 * */
