
#include "HGCPlotting.h"
#include <cmath> 


//Creates List of radii over which the program can run
const int NumberOfEntries = 100;
double LowLimit = 0;
double HighLimit = 0.06;
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
for (unsigned j=0; j < NumberOfEntries; j++){
	double Value = round((LowLimit + j*(HighLimit-LowLimit)/NumberOfEntries)*1000)/1000 ;
	RadiiList[j] = Value;
}
  if ( name == "TriggerCells" ){

    _cloned_hists[ name ] [ "tc_eta" ] = new TH1D ( (name + "_tc_eta").c_str(), "", 100,-5,5 );  
    _cloned_hists[ name ] [ "tc_phi" ] = new TH1D ( (name + "_tc_phi").c_str(), "", 100,-M_PI,M_PI );  // M_PI == Pi 

  }/*else if ( name == "PU0_General" ){
		  
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
											
													
													       //
													        
													   } //else if (name="single_event") {
													          //_clone_2d_map[ name ] [ "scatter_norm" ] = new TH2D ((name+"_scatter_norm").c_str(), ) //scatter in normalised co-ordinates
													 }*/
	else if ( name == "PU0_forward" || name == "PU0_backward" ){
		for (int j=0; j < NumberOfEntries; j++) {
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] = new TH1D ( (name+HistName).c_str(), "", 300, 0,88);
		}
  } }


/* POSITION RESOLUTION */
// find x/z , y/z for each TC, then...
// Truth: 
  // using: gen_phi, gen_eta
  // find x/z y/z for truth means projecting it onto the first layer
  // 	eta == - ln [ tan(theta/2) ]
  // where theta is the angle with the z-axis
  // Also
  // 	tan(theta) = rho/z (for all rho,z) where rho is in the x-y plane
  // Hence using tan(2A) == 2tan(A)/(1-tan^2(A)),
  // 	exp{-eta} = tan(theta/2) 
  // 	tan(theta) = 2exp{-eta}/(1-exp(-2*eta)) = rho/z
  // Also
  // 	x/z = rho/z cos(phi)
  // 	y/z = rho/z sin(phi)
  // Therefore: (TRUTH)
  // 	x/z = 2exp{-eta}/(1-exp(-2eta)) cos(phi)
  // 	y/z = 2exp{-eta}/(1-exp(-2eta)) sin(phi)
  // Let X = x/z, Y = y/z, // USING normalised co-ordinates..
  //
  // Draw a circle of radius R=(r/z)~0.3 around truth co-ordinates. 
  // Obtain from this a list of TCs inside the circle
  // 	CALCULATE
  // 		X_cluster = sum_i(E_i X_i) / sum_i(E_i)
 // 		E_cluster = sum_i(E_i)
  // 		   .: E_c X_c = sum_i(E_i X_i) //check
  // 		   	
  // 	 	Also E_X = E_C sin(theta_c) cos(phi_c) // Projection onto C coordinate
  


void HGCPlotting::CalculateTriggerCellVariables() {
	//std::cout << tc_zside->size() << " " << tc_layer->size() <<" "<<tc_z->size()<<std::endl;
    /* ENERGY RESOLUTION */
    //Ex and Ey sums
	double exsum_forward = 0;
	double eysum_forward = 0;
	double exsum_backward = 0;
	double eysum_backward = 0;
  
	for (unsigned int i = 0; i < tc_pt->size(); i++){
		// IF-statement to filter between forward and backward calorimeters...
		if  ( tc_eta->at(i)>0){
		//FORWARD
		exsum_forward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_forward += tc_pt->at(i)*std::sin(tc_phi->at(i));

		} else {
		//BACKWARD
		exsum_backward += tc_pt->at(i)*std::cos(tc_phi->at(i));
		eysum_backward += tc_pt->at(i)*std::sin(tc_phi->at(i));
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
	_event_variables["xnft"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(0));  //std::sin(2 * std::atan2(std::exp( -gen_eta->at(0) ),1. ) ) * std::cos(gen_phi->at(0));
	_event_variables["ynft"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(0));  //std::sin(2 * std::atan2(std::exp( -gen_eta->at(0) ),1. ) ) * std::sin(gen_phi->at(0));
	_event_variables["xnbt"] = std::cos(gen_phi->at(0)) / sinh(gen_eta->at(1));  //std::sin(2 * std::atan2(std::exp( -gen_eta->at(1) ),1. ) ) * std::cos(gen_phi->at(1));
	_event_variables["ynbt"] = std::sin(gen_phi->at(0)) / sinh(gen_eta->at(1));  //std::sin(2 * std::atan2(std::exp( -gen_eta->at(1) ),1. ) ) * std::sin(gen_phi->at(1)); 


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
			//_event_details["ptfx"].push_back(tmpptx);
			//_event_details["ptfy"].push_back(tmppty);
		
			/* CHECK IF TC IS IN RADIUS R
			 * (X-X_t)^2 + (Y-Y_t)^2 < R^2 */
			if ( (tmpx-_event_variables["xnft"])*(tmpx-_event_variables["xnft"]) 
					+ (tmpy - _event_variables["ynft"])*(tmpy - _event_variables["ynft"]) < R*R ) {
				
				/*Add TC variables to vectors */
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
				//_event_details["xnbc_pt"].push_back( tmppt*std::cos(tmpphi));
				//_event_details["ynbc_pt"].push_back( tmppt*std::sin(tmpphi));
                //_event_details["ptfxc"].push_back( tmpptx );
	            //_event_details["ptfyc"].push_back( tmppty );
	         
                //_event_variables["fX_weighted_pt"] += tmpx*tmpptx;  //tmppt*std::cos(tmpphi);
                //_event_variables["fY_weighted_pt"] += tmpy*tmppty;//tmppt*std::sin(tmpphi);
                //_event_variables["fX_sum"] +=tmpptx;
                //_event_variables["fY_sum"] +=tmppty;
                _event_variables["fE_sum"] +=tmppt;

				/*Add position resolution calculations here to form sum of pt.x */
				_event_variables["bX_weighted_Et"] += tmppt*tmpx;
				_event_variables["bY_weighted_Et"] += tmppt*tmpy;
				
				//_event_variables["bX_sum"] += tmppt;
				//_event_variables["bY_sum"] += tmppt;
				
				_event_variables["bE_sum"] += tmppt;
				//std::cout << "pt: " << tmppt << std::endl;	
			} 
		}
	} /* end of loop*/
	/*This section computes the weighted averages for both the one weighted by ptx(transverse energy/momentum) in the forward calorimeter*/
	/*
	_event_variables["fX_weighted_pt"] /= _event_variables["fX_sum"];
    _event_variables["fY_weighted_pt"] /= _event_variables["fY_sum"]; 
    _event_variables["fd_pos"] = std::sqrt((_event_variables["fX_weighted_pt"] - _event_variables["xnft"])*(_event_variables["fX_weighted_pt"] - _event_va
	*/

	//This section does the same for the backwards calorimeter. The pt calculations are not yet included in the backwards calorimeter, and NEED TO BE ADDE
	/*
	_event_variables["bX_weighted_pt"] /= _event_variables["bX_sum"];
	_event_variables["bY_weighted_pt"] /= _event_variables["bY_sum"]; 
	_event_variables["bd_pos"] = std::sqrt((_event_variables["bX_weighted_pt"] - _event_variables["xnbt"])*(_event_variables["bX_weighted_pt"] - _event_va
	*/
   
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
	_event_variables["fd_energy_R"] = _event_variables["fE_sum"]; //- gen_pt->at(0);
	_event_variables["bd_energy_R"] = _event_variables["bE_sum"]; //- gen_pt->at(1);
}


void HGCPlotting::FillAllHists( std::string name ){

  /* RUN FOR EVERY __name__ IN EVERY EVENT*/

  CalculateTriggerCellVariables();

  double r = 0.04;
  // Run calcs for radius r about truth values, remember, SINGLE EVENT at a time, saved to _event_variables.
  // If run multiple times, will currently replace itself for each R
  CalculateReducedCircle(r); 
  //CalculateReducedCircle(r-0.02);  
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
		for (int j=0; j<NumberOfEntries; j++) {
			CalculateReducedCircle(RadiiList[j]); 
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] ->Fill( _event_variables["fd_energy_R"]);

		}
		//     _cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_forward"  ] );
		//     _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_forward"  ] );
		//     _cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_forward"  ] );
		//     //std::cout << _event_variables["er_sum_forward"] << std::endl;//yoyo-db  
		//     _cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_forward"  ] );        
		//     _cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_forward"  ] );    
		//     _cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_forward"] );
		//        // POSITION RES. FILLED HERE
		//     _cloned_hists[ name ] [ "dpos_2" ] ->Fill(_event_variables["fd_pos_E"]);
		//     _cloned_hists[ name ] [ "dpos" ] ->Fill(_event_variables["fd_pos"]);
		//     _cloned_hists[ name ] [ "dpos_X" ] ->Fill(_event_variables["fX_weighted_pt"]- _event_variables["xnft"]);
		//     _cloned_hists[ name ] [ "dpos_Y" ] ->Fill(_event_variables["fY_weighted_pt"] - _event_variables["ynft"]);
		//     _cloned_hists[ name ] [ "dpos_X_E" ] ->Fill(_event_variables["fX_weighted_Et"]- _event_variables["xnft"]);
		//     _cloned_hists[ name ] [ "dpos_Y_E" ] ->Fill(_event_variables["fY_weighted_Et"] - _event_variables["yn



  } else if ( name == "PU0_backward" ){
		for (int j=0; j<NumberOfEntries; j++) {
			CalculateReducedCircle(RadiiList[j]); 
			char o[50];
			sprintf(o, "%f", RadiiList[j]);
			o[5] = '\0';
			std::string p = o;
			std::string HistName  = "_denergy_R_" + p;
			_cloned_hists[ name] [ HistName  ] ->Fill( _event_variables["bd_energy_R"]);

		}
		//    _cloned_hists[ name ] [ "tc_n" ] ->Fill ( tc_n );
		//    _cloned_hists[ name ] [ "ex_sum" ] ->Fill (  _event_variables[  "ex_sum_backward"  ] );
		//    _cloned_hists[ name ] [ "ey_sum" ] ->Fill (  _event_variables[  "ey_sum_backward"  ] );
		//    _cloned_hists[ name ] [ "er_sum" ] ->Fill (  _event_variables[  "er_sum_backward"  ] );
		//    _cloned_hists[ name ] [ "ephi_sum" ] ->Fill (  _event_variables[  "ephi_sum_backward"  ] );   
		//    _cloned_hists[ name ] [ "dphi_met" ] ->Fill (  _event_variables[  "dphi_met_backward"  ] );                
		//    _cloned_hists[ name ] [ "denergy" ] ->Fill ( _event_variables[ "denergy_backward"] );
		//    _cloned_hists[ name ] [ "dpos" ] ->Fill( _event_variables["bd_pos"])


  }
}


