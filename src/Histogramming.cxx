
#include "HGCPlotting.h"
#include <cmath> 

#include <candidate.h>
#include <algorithm> // For remove_if


//Creates List of radii over which the program can run
const int NumberOfEntries = 40;
double LowLimit = 0;
double HighLimit = 0.08;
double RadiiList [NumberOfEntries] = {};


void HGCPlotting::MakeAllHists( std::vector<std::string> &HistoSets){

	std::cout << "Creating All Histograms" << std::endl;

	for (auto& names : HistoSets ){
            LoadHistoTemplates ( names );
    }  

}

	

// TH1D Syntax: new TH1D(name, title, nbinsx, xlow, xhigh) ; xlow-> edge of lowest bin, xhigh -> edge of highest bin
// or: new TH1D(name, title, nbinsx, xbins) ; xbins -> low edge of of everybin, contains nbinsx+1 entries. 

void HGCPlotting::LoadHistoTemplates( std::string name ) { 
	/* Initialise all empty plots */
	
	// Code to run multiple Rs (stephans)
	/*
	for (unsigned j=0; j < NumberOfEntries; j++) {
		double Value = round((LowLimit + j*(HighLimit-LowLimit)/NumberOfEntries)*1000)/1000 ;
		RadiiList[j] = Value;
	}
	*/		
	
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
		_cloned_hists[ name ] [ "dpos" ] = new TH1D ( (name + "_dpos").c_str(), "", 150, 0,.1);
		_cloned_hists[ name ] [ "dpos_2" ] = new TH1D ( (name + "_dpos_E").c_str(), "", 150, 0,.1);
		_cloned_hists[ name ] [ "dpos_X" ] = new TH1D ( (name + "_dpos_x").c_str(), "", 150, -.04,.04);
		_cloned_hists[ name ] [ "dpos_Y" ] = new TH1D ( (name + "_dpos_y").c_str(), "", 150, -.04,.04);
		_cloned_hists[ name ] [ "dpos_X_E" ] = new TH1D ( (name + "_dpos_x_E").c_str(), "", 150, -.04,.04);
		_cloned_hists[ name ] [ "dpos_Y_E" ] = new TH1D ( (name + "_dpos_y_E").c_str(), "", 150, -.04,.04);
		_cloned_hists[ name] [ "denergy_R" ] = new TH1D ( (name+"_denergy_R").c_str(), "", 150, -25.,25.);											   
	}

	//} else if (name="single_event") {
	//scatter in normalised co-ordinates
	//	_clone_2d_map[ name ] [ "scatter_norm" ] = new TH2D ((name+"_scatter_norm").c_str(), )	
	//}
	// Multiple R code
	/*
	else if ( name == "PU0_forward" || name == "PU0_backward" ){
		for (int j=0; j < NumberOfEntries; j++) {
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] = new TH1D ( (name+HistName).c_str(), "", 150, -10,10);
		}
  }
  */	

}




void HGCPlotting::CalculateTriggerCellVariables() {
	//std::cout << tc_zside->size() << " " << tc_layer->size() <<" "<<tc_z->size()<<std::endl;
    /* ENERGY RESOLUTION */
    //Ex and Ey sums
	double exsum_forward = 0;
	double eysum_forward = 0;
	double exsum_backward = 0;
	double eysum_backward = 0;
  
	// intialise _Event_details by deleting
	_event_details.clear();
	
	// LOOP OVER ALL ENTRIES
	for (unsigned int i = 0; i < tc_pt->size(); i++){
		// IF-statement to filter between forward and backward calorimeters...
		if  ( tc_eta->at(i)>0){
		//FORWARD
		exsum_forward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_forward += tc_pt->at(i)*std::sin(tc_phi->at(i));
		
		_event_details["fX"].push_back( tc_x->at(i) / tc_z->at(i) ); 
		_event_details["fY"].push_back( tc_y->at(i) / tc_z->at(i) );
		_event_details["fP"].push_back( tc_pt->at(i) ); 
		
		} else {
		//BACKWARD
		exsum_backward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_backward += tc_pt->at(i)*std::sin(tc_phi->at(i));
		
		//
		_event_details["bX"].push_back( tc_x->at(i) / tc_z->at(i) );
		_event_details["bY"].push_back( tc_y->at(i) / tc_z->at(i) );
		_event_details["bP"].push_back( tc_pt->at(i) );

		}
	}

	double ersum_forward = std::sqrt( exsum_forward*exsum_forward + eysum_forward*eysum_forward ); // forward energy sum in r

	double ersum_backward = std::sqrt( exsum_backward*exsum_backward + eysum_backward*eysum_backward ); // backward energy sum in r
	//std::cout << "er_sum_forward " << ersum_forward << std::endl; //yoyo-db  
  
	// std::atan2 -> arctan(y/x), uses signs of y,x to find the correct quadrant!
	// sum of energy in phi
	double ephisum_forward = std::atan2( eysum_forward, exsum_forward ); 
    double ephisum_backward = std::atan2( eysum_backward, eysum_backward );  

	_event_variables[  "ex_sum_forward"  ] = exsum_forward;
	_event_variables[  "ey_sum_forward"  ] = eysum_forward;
	_event_variables[  "er_sum_forward"  ] = ersum_forward;
	_event_variables[  "ephi_sum_forward"  ] = ephisum_forward;

	_event_variables[  "ex_sum_backward"  ] = exsum_backward;
	_event_variables[  "ey_sum_backward"  ] = eysum_backward;
	_event_variables[  "er_sum_backward"  ] = ersum_backward;
	_event_variables[  "ephi_sum_backward"  ] = ephisum_backward;

	//Difference betwen photon truth phi and MET
	TVector2 met, truth; // init measured, truth 2-d (x-y) vectors
	met.Set( exsum_forward, eysum_forward ); 
	truth.SetMagPhi( 1, gen_phi->at(0) ); // gen_phi is the TRUTH VALUE for phi.
	// Find angle (in phi) between met & truth; Named dphi_met (delta-phi)
	double dphi_met = met.DeltaPhi( truth );
	_event_variables[  "dphi_met_forward"  ] =  dphi_met;
	// SAME for backward!!
	met.Set( exsum_backward, eysum_backward );
	truth.SetMagPhi( 1, gen_phi->at(1) );
	dphi_met = met.DeltaPhi( truth );
	 _event_variables[  "dphi_met_backward"  ] =  dphi_met;
 
	//delta-energy (denergy), need forward and back! 
	// CHANGE to er_sum MINUS pen_pt ... 
	_event_variables["denergy_forward"] = ersum_forward - gen_pt->at(0);
	_event_variables["denergy_backward"] = ersum_backward - gen_pt->at(1);

	/*TRUTH VALUES*/ 
	// in our normalised co-ordinates  
	_event_variables["xnft"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(0));
	_event_variables["ynft"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(0));
	_event_variables["xnbt"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(1)); 
	_event_variables["ynbt"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(1));  

	/* Debug for position Resolution */
        //std::cout <<  "T: "<<_event_variables["xnbt" ]<< " | "<< _event_variables["ynbt"] << std::endl; 
	//std::cout << "C: "<< _event_variables["bX_weighted_pt"]<< " | " <<_event_variables["bY_weighted_pt"] << std::endl;
	//std::cout << _event_variables["bX_sum"]<<std::endl;
	//std::cout << std::endl;    


}



void HGCPlotting::CalculateReducedCircle(const double& R) {
	/* Calculates values for energy/ position resolution in the circle of radius R = r/z_0 */
	
	// Loop to create circle and _event_details
	//double av_pt_xy=0; // expectation value of pt (?)
	// init for the weighted sums
	_event_variables["fX_weighted_pt"] = 0;
	_event_variables["fY_weighted_pt"] = 0;
    _event_variables["fY_sum"] = 0 ;
    _event_variables["fX_sum"] = 0;	
    _event_variables["bX_weighted_pt"] = 0;
	_event_variables["bY_weighted_pt"] = 0;
	_event_variables["fY_weighted_Et"] = 0;
	_event_variables["fX_weighted_Et"] = 0;
	_event_variables["bX_sum"] = 0;
	_event_variables["bY_sum"] = 0; 
	_event_variables["fE_sum"] = 0;
    _event_variables["bE_sum"] = 0; 

	// Loop over all TCs in a single event
	for (unsigned j=0; j < tc_pt->size(); j++){
		
		//double tmpphi = tc_eta->at(j);   
		double tmpx,tmpy,tmppt; //,tmpptx,tmppty;

		tmpx = tc_x->at(j) / tc_z->at(j);
		tmpy = tc_y->at(j) / tc_z->at(j);
		tmppt = tc_pt->at(j);

		/* pt in X and Y directions*/
		//tmpptx = tmppt*std::cos(tmpphi); 
        //tmppty = tmppt*std::sin(tmpphi);

		// Choose Between Forward and backward
		if (tc_eta->at(j) > 0) { /* FORWARD: Implemented sum(p_x.x)/sum(p_x)*/
 
			//std::cout << tmpx << "|" << tmpy << "|" << tmppt << std::endl;
			//std::cout <<"tc_x"<< tc_x->at(j) <<std::endl;
	 
			_event_details["xnf"].push_back(tmpx); 
			_event_details["ynf"].push_back(tmpy);
			_event_details["ptf"].push_back(tmppt);
		
			/* CHECK IF TC IS IN RADIUS R
			 * (X-X_t)^2 + (Y-Y_t)^2 < R^2 */
			if ( (tmpx-_event_variables["xnft"])*(tmpx-_event_variables["xnft"]) 
					+ (tmpy - _event_variables["ynft"])*(tmpy - _event_variables["ynft"]) < R*R ) {
				
				/*Add TC variables to vectors */
				/* Fill Candidate Arrays */
				_event_details["xnfc"].push_back( tmpx );
				_event_details["ynfc"].push_back( tmpy );
				_event_details["ptfc"].push_back( tmppt );
		  		
				_event_variables["fE_sum"] +=tmppt;
				_event_variables["fX_weighted_Et"] += tmpx*tmppt;
				_event_variables["fY_weighted_Et"] += tmpy*tmppt;
			}

		} else { /* Backward: IMPLEMENTED sum(pt.x)/sum(pt)*/
			/*Invert co-ordinates for no good reason*/
			tmpx = -tmpx; tmpy = -tmpy;
             
			_event_details["xnb"].push_back(tmpx);
			_event_details["ynb"].push_back(tmpy);
			_event_details["ptb"].push_back(tmppt);
        	
            //std::cout << "tmp x,y: " << tmpx<<","<<tmpy<<std::endl;
			//std::cout << "e.v. TRUTH: x,y: "<<_event_variables["xnbt"]<<","<<_event_variables["ynbt"]<<std::endl;
			
			/* Check if in radius */
			if ( (tmpx - _event_variables["xnbt"])*(tmpx - _event_variables["xnbt"]) 
					+(tmpy - _event_variables["ynbt"])*(tmpy - _event_variables["ynbt"]) < R*R ) {
				/* Fill Candidate Arrays */
				_event_details["xnbc"].push_back(tmpx);
				_event_details["ynbc"].push_back(tmpy);
				_event_details["ptbc"].push_back(tmppt);

				/*Add position resolution calculations here to form sum of pt.x */
				_event_variables["bX_weighted_Et"] += tmppt*tmpx;
				_event_variables["bY_weighted_Et"] += tmppt*tmpy;	
				_event_variables["bE_sum"] += tmppt;
			} 
		}
	} /* end of loop*/
	
	/* Forward Calorimeter XY-position resolution weighted by energy*/
	_event_variables["fX_weighted_Et"] /= _event_variables["fE_sum"];
	_event_variables["fY_weighted_Et"] /= _event_variables["fE_sum"]; 	
	/* Forward radial resolution */	
	_event_variables["fd_pos_E"] = std::sqrt( 
			(_event_variables["fX_weighted_Et"]-_event_variables["xnft"]) 
			*(_event_variables["fX_weighted_Et"] - _event_variables["xnft"]) 
			+ (_event_variables["fY_weighted_Et"] - _event_variables["ynft"])
			*(_event_variables["fY_weighted_Et"] - _event_variables["ynft"]) ); 
    	
	/* Backward Calorimeter XY-position resolution weighted by energy*/
	_event_variables["bX_weighted_Et"] /= _event_variables["bE_sum"];
	_event_variables["bY_weighted_Et"] /= _event_variables["bE_sum"];
	/* Backward radial resolution */	
	_event_variables["bd_pos_E"] = std::sqrt(
			(_event_variables["bX_weighted_Et"] - _event_variables["xnbt"])
			*(_event_variables["bX_weighted_Et"] - _event_variables["xnbt"])
			+(_event_variables["bY_weighted_Et"] - _event_variables["ynbt"])
			*(_event_variables["bY_weighted_Et"] - _event_variables["ynbt"]) ); 

    /*difference in energy for circle*/
	_event_variables["fd_energy_R"] = _event_variables["fE_sum"] - gen_pt->at(0);
	_event_variables["bd_energy_R"] = _event_variables["bE_sum"] - gen_pt->at(1);

}


Candidate::Candidate(const double& XT, const double& YT, const double& ET) : 
	_XT(XT), _YT(YT), _ET(ET) 
{	
}

void Candidate::importDetails(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& P) {
	//
	entry tmp;
	auto xi = X.begin();
	auto yi = Y.begin();
	auto pi = P.begin(); 
	while (xi != X.end()) {
		tmp.x = *xi; tmp.y = *yi; tmp.p = *pi; 
		_Entries.push_back(tmp);
		xi++;
		yi++;
		pi++;
	}
}



struct inCircle {
    double _R, _XT,_YT;
    	
	// Params
	inCircle(const double& R, const double& XT, const double& YT): _R(R), _XT(XT), _YT(YT) {}
	bool operator()(const entry& E) const {
		// Comparision here
		return ( ((E.x-_XT)*(E.x-_XT) + (E.y-_YT)*(E.y-_YT)) > _R*_R ); 
	}
}; 


void Candidate::crop(const double& R) {
	// Crop _Entries  such that they are constrained in the circle defined by R
	// Have vector _Entries
	_R = R;
	
	//MarkEntries for deletion
	
	
	//std::cout << "No. of entries: " << _Entries.size() << "\t"; 
	// Delete entries which are outside circle
	_Entries.erase( 
			std::remove_if ( 
				_Entries.begin(), _Entries.end(), inCircle(_R, _XT, _YT) 
				),
		   _Entries.end() );	
	//std::cout << "Croped size: " << _Entries.size() << "\n"; 
} 

void Candidate::calculate_resolutions() {
	/*Calculate the e_res, x_res, y_res, e_sum*/

	// Initialise result doubles to zero.
	e_res =0; x_res=0; y_res=0; r_res=0; e_sum=0;

	for(auto const& e: _Entries) {
		// Loop over all candidate entries
		
		// sum energy
		e_sum += e.p;
		
		// sum weighted position
		x_res += e.p*e.x;
		y_res += e.p*e.y;

	}
	// Divide by weights to obtain weighted average
	x_res /= e_sum; y_res /= e_sum;


	// Subtract the TRUTH values to obtain resolution for single event
	x_res -= _XT;
	y_res -= _YT;
	e_res -= _ET;

	r_res = std::sqrt(x_res*x_res + y_res*y_res);
}

double Candidate::getERes() { return e_res;} 
double Candidate::getXRes() { return x_res;}
double Candidate::getYRes() { return y_res;}
double Candidate::getRRes() { return r_res;}
double Candidate::getESum() { return e_sum;} 

void HGCPlotting::FillAllHists( std::string name ){
  /* RUN FOR EVERY __name__ IN EVERY EVENT*/

  CalculateTriggerCellVariables();

  //double r = 0.04;
  // Run calcs for radius r about truth values, remember, SINGLE EVENT at a time, saved to _event_variables.
  // If run multiple times, will currently replace itself for each R
  //CalculateReducedCircle(r); 
 
  /*
  std::cout << _event_variables["fX_weighted_Et"]-_event_variables["xnft"] << "\t";
  // Get the data into the required form

  // Initialise Candidate instance 
  Candidate fCand(_event_variables["xnft"], _event_variables["ynft"], gen_pt->at(0));
  fCand.importDetails(_event_details["fX"], _event_details["fY"], _event_details["fP"]);
  fCand.crop(r);  
  
  fCand.calculate_resolutions();	

  std::cout << fCand.getXRes() << "\n";
  std::cout << _event_variables["fd_energy_R"] << "\t";
  std::cout << fCand.getESum() << "\n \n";
 
  Candidate bCand(_event_variables["xnbt"], _event_variables["ynbt"], gen_pt->at(0));
  bCand.importDetails(_event_details["bX"], _event_details["bY"], _event_details["bP"]);
  bCand.crop(r);
  */
 

  
  std::pair<double,double> r_lims(0.1, 0.01); 
  unsigned r_num = 10;
  double r_curr;
  double r_inc = (r_lims.first/r_lims.second)/r_num; 

  std::map<unsigned, std::vector<double>> data_r;  // r : vector<E_sum>/<>

  for (unsigned i = 0; i < r_num ;++i ) {
	  // Initialise with vars
	  Candidate fCand(_event_variables["xnft"], _event_variables["ynft"], gen_pt->at(0)); 
	  Candidate bCand(_event_variables["xnbt"], _event_variables["ynbt"], gen_pt->at(0));
	  // Import event details (readout)
	  fCand.importDetails(_event_details["fX"], _event_details["fY"], _event_details["fP"]);
	  fCand.importDetails(_event_details["bX"], _event_details["bY"], _event_details["bP"]);

	  // 
	  fCand.crop(r_curr); fCand.crop(r_curr);
	  // Now we can get the data that we want 	  
       

	  data_r[i].push_back(fCand.getERes());
	  //bCand.getERes(); 

	  // Increment loop (important!)
	  r_curr += r_inc;
  }


  /*TODO 
   *	Implement methods for:
   *		- d-position (radial,x,y)
   *		- d-energy (E_sum - E_truth)
   * */


  // Implement Scheme to loop over varying R, 
  /*
	unsigned r_num = 10;
	pair<double,double> r_range;
	double r_curr;
	r_range.first = 0.1;
	r_range.second = 0.01;	
	r_curr = r_range.first;
	for (unsigned i=0; i < r_num;i++) {

		r_curr -= (r_range.first - r_range.second)/r_num ;
		
		// create candidates for f/b; 
		// calculate res
		// store res for specified r_check 

		}
   */

    
  //  if ( name == "PU0" ||  name == "PU200" )

  if ( name == "TriggerCells" ){
    // Visualise TCs...
	// 
    for (unsigned int i = 0; i < tc_eta->size(); i++){
      _cloned_hists[ name ] [ "tc_eta" ] ->Fill ( tc_eta->at(i) );
      _cloned_hists[ name ] [ "tc_phi" ] ->Fill ( tc_phi->at(i) );
    }
  } else if ( name == "PU0_General" ){
    	_cloned_hists[ name ] [ "tc_n" ] ->Fill ( tc_n );
  } else if ( name == "PU0_forward" ){

	    /*
		for (int j=0; j<NumberOfEntries; j++) {
			CalculateReducedCircle(RadiiList[j]); 
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] ->Fill( _event_variables["fd_energy_R"]);
		}
		*/

			// Histograms for entire event, sums over all TCs
		     _cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_forward"  ] );
		     _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_forward"  ] );
		     _cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_forward"  ] );
		     //std::cout << _event_variables["er_sum_forward"] << std::endl;//yoyo-db  
		     _cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_forward"  ] );        
		     _cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_forward"  ] );    
		     _cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_forward"] );
		        
			 // Radius specific histograms, position, energy resolutions
			 _cloned_hists[ name ] [ "dpos_2" ] ->Fill(_event_variables["fd_pos_E"]);
			 //_cloned_hists[ name ] [ "dpos" ] ->Fill(_event_variables["fd_pos"]);
			 // Redacted momentum weight
			 //_cloned_hists[ name ] [ "dpos_X" ] ->Fill(_event_variables["fX_weighted_pt"]- _event_variables["xnft"]);
		     //_cloned_hists[ name ] [ "dpos_Y" ] ->Fill(_event_variables["fY_weighted_pt"] - _event_variables["ynft"]);
			 _cloned_hists[ name ] [ "dpos_X_E" ] ->Fill(_event_variables["fX_weighted_Et"]- _event_variables["xnft"]);
		     _cloned_hists[ name ] [ "dpos_Y_E" ] ->Fill(_event_variables["fY_weighted_Et"] - _event_variables["ynft"]);

  } else if ( name == "PU0_backward" ){
	    /*
		for (int j=0; j<NumberOfEntries; j++) {
			CalculateReducedCircle(RadiiList[j]); 
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] ->Fill( _event_variables["bd_energy_R"]);
		}
		*/
		    _cloned_hists[ name ] [ "tc_n" ] ->Fill ( tc_n );
		    _cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_backward"  ] );
		    _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_backward"  ] );
		    _cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_backward"  ] );
		    _cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_backward"  ] );   
		    _cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_backward"  ] );                
		    _cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_backward"] );
		    //_cloned_hists[ name ] [ "dpos" ] ->Fill( _event_variables["bd_pos"]);
			
			// Radius specific
			//_cloned_hist [ name ] [ HIST NAME ] ->Fill ([_event_variables["bd_pos_E"]]); 
  
  }

}


