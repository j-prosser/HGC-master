#include "TROOT.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "stdio.h"



TGraph* getGraph(TFile *file, std::string name = "Graph") {
		TGraph *graph = (TGraph*) file->Get(name.c_str());
		return graph;
}

void setAxis(TMultiGraph *mg) {

    // implement latex #frac{2s}{#pi#alpha^{2}} 
		mg->GetXaxis()->SetTitle("Cone Size #left(R = #frac{r}{z} #right) ");
		mg->GetYaxis()->SetTitle("Energy Width #left( #frac{#sigma_{E}}{<E>} #right)");
        mg->GetHistogram()->SetTitle("Energy Widths for Varying Radius n=10000");
}


void plot(int firstbin, int secondbin, int thirdbin, int total, int plotzero) {
		TCanvas *c1 = new TCanvas("c1", "multigraph", 700, 500);
		
        c1->SetGrid(); // display the grid
		
        TMultiGraph *mg_1 = new TMultiGraph();

		std::string filepath_0  = "output/PU0/output.root";
		std::string filepath_200  = "output/PU200/output.root";

		TFile *file_0 = new TFile(filepath_0.c_str(), "READ");
		
        std::string graphname_1 = "sigma_EE_r_1.47_1.98";
		std::string graphname_2 = "sigma_EE_r_1.98_2.49";
		std::string graphname_3 = "sigma_EE_r_2.49_3.00";
		std::string graphname_4 = "Sigma_EE_r";
        
		// Graph names are the legend labels!
        std::string leg_graphname_1 = "PU200 1.47<eta<1.98";
		std::string leg_graphname_2 = "PU200 1.98<eta<2.49";
		std::string leg_graphname_3 = "PU200 2.49<eta<3.00";
		std::string leg_graphname_4 = "PU200";
		
        std::string leg0_graphname_1 = "PU0 1.47<eta<1.98";
		std::string leg0_graphname_2 = "PU0 1.98<eta<2.49";
		std::string leg0_graphname_3 = "PU0 2.49<eta<3.00";
		std::string leg0_graphname_4 = "PU0";
        
        TGraph *Graph_0_1 = getGraph(file_0, graphname_1.c_str());
		Graph_0_1->SetMarkerStyle(2);
        //Graph_0_1->SetLineColor(kBlue);
        Graph_0_1->SetLineColor(11);
		Graph_0_1->SetTitle(leg0_graphname_1.c_str());
		Graph_0_1->SetMarkerColor(11);
		
        TGraph *Graph_0_2 = getGraph(file_0, graphname_2.c_str());
		Graph_0_2->SetMarkerStyle(2);
        Graph_0_2->SetLineColor(21);
		Graph_0_2->SetTitle(leg0_graphname_2.c_str());
		Graph_0_2->SetMarkerColor(21);
		TGraph *Graph_0_3 = getGraph(file_0, graphname_3.c_str());
		Graph_0_3->SetMarkerStyle(2);
        Graph_0_3->SetLineColor(31);
		Graph_0_3->SetTitle(leg0_graphname_3.c_str());
		Graph_0_3->SetMarkerColor(31);
		TGraph *Graph_0_4 = getGraph(file_0, graphname_4.c_str());
		Graph_0_4->SetMarkerStyle(2);
        Graph_0_4->SetLineColor(41);
		Graph_0_4->SetTitle(leg0_graphname_4.c_str());
		Graph_0_4->SetMarkerColor(41);

		TFile *file_200 = new TFile(filepath_200.c_str(), "READ");

		TGraph *Graph_200_1 = getGraph(file_200, graphname_1.c_str());
		Graph_200_1->SetLineColor(2);
		Graph_200_1->SetTitle(leg_graphname_1.c_str());
		Graph_200_1->SetMarkerColor(2);
		Graph_200_1->SetMarkerStyle(5);

		TGraph *Graph_200_2 = getGraph(file_200, graphname_2.c_str());
		Graph_200_2->SetTitle(leg_graphname_2.c_str());
		Graph_200_2->SetLineColor(3);
		Graph_200_2->SetMarkerColor(3);
		
		Graph_200_2->SetMarkerStyle(5);
        TGraph *Graph_200_3 = getGraph(file_200, graphname_3.c_str());
		Graph_200_3->SetTitle(leg_graphname_3.c_str());
		Graph_200_3->SetLineColor(4);
		Graph_200_3->SetMarkerColor(4);
		Graph_200_3->SetMarkerStyle(5);
		
        TGraph *Graph_200_4 = getGraph(file_200, graphname_4.c_str());
		Graph_200_4->SetTitle(leg_graphname_4.c_str());
		Graph_200_4->SetLineColor(5);
		Graph_200_4->SetMarkerColor(5);
		Graph_200_4->SetMarkerStyle(5);


		if (firstbin) {
			if (plotzero) {
			mg_1->Add(Graph_0_1);
			}
			mg_1->Add(Graph_200_1);
		}
		if (secondbin) {
			if (plotzero) {
				mg_1->Add(Graph_0_2);
			}
			mg_1->Add(Graph_200_2);
		}
		if (thirdbin) {
			if (plotzero) {
				mg_1->Add(Graph_0_3);
			}
			mg_1->Add(Graph_200_3);
		}
		if (total) {
			if (plotzero) {
				mg_1->Add(Graph_0_4);
			}
			mg_1->Add(Graph_200_4);
		}
		mg_1->Draw("AL*"); /// prev: ac*
		setAxis(mg_1);
		c1->BuildLegend(.9, .21, .9, .21); 
}
void ploteta() {
		plot(1,1,1,1,1);
}
