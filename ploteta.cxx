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
		mg->GetXaxis()->SetTitle("Radius (Reduced coordinates) ");
		mg->GetYaxis()->SetTitle("Sigma_E/E ");
}


void plot(int firstbin, int secondbin, int thirdbin, int total, int plotzero) {
		TCanvas *c1 = new TCanvas("c1", "multigraph", 700, 500);
		c1->SetGrid();
		TMultiGraph *mg_1 = new TMultiGraph();

		std::string filepath_0  = "output/PU0/output.root";
		std::string filepath_200  = "output/PU200/output.root";

		TFile *file_0 = new TFile(filepath_0.c_str(), "READ");
		std::string graphname_1 = "sigma_EE_r_1.47_1.98";
		std::string graphname_2 = "sigma_EE_r_1.98_2.49";
		std::string graphname_3 = "sigma_EE_r_2.49_3.00";
		std::string graphname_4 = "Sigma_EE_r";

		TGraph *Graph_0_1 = getGraph(file_0, graphname_1.c_str());
		Graph_0_1->SetLineColor(kBlue);
		TGraph *Graph_0_2 = getGraph(file_0, graphname_2.c_str());
		TGraph *Graph_0_3 = getGraph(file_0, graphname_3.c_str());
		TGraph *Graph_0_4 = getGraph(file_0, graphname_4.c_str());
		
		TFile *file_200 = new TFile(filepath_200.c_str(), "READ");

		TGraph *Graph_200_1 = getGraph(file_200, graphname_1.c_str());
		Graph_200_1->SetLineColor(2);
		Graph_200_1->SetTitle(graphname_1.c_str());
		Graph_200_1->SetMarkerColor(2);
		TGraph *Graph_200_2 = getGraph(file_200, graphname_2.c_str());
		Graph_200_2->SetTitle(graphname_2.c_str());
		Graph_200_2->SetLineColor(3);
		Graph_200_2->SetMarkerColor(3);
		TGraph *Graph_200_3 = getGraph(file_200, graphname_3.c_str());
		Graph_200_3->SetTitle(graphname_3.c_str());
		Graph_200_3->SetLineColor(4);
		Graph_200_3->SetMarkerColor(4);
		TGraph *Graph_200_4 = getGraph(file_200, graphname_4.c_str());
		Graph_200_4->SetTitle(graphname_4.c_str());
		Graph_200_4->SetLineColor(5);
		Graph_200_4->SetMarkerColor(5);


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
		mg_1->Draw("ac*");
		setAxis(mg_1);
		c1->BuildLegend(.9, .21, .9, .21); 
}
void ploteta() {
		plot(1,1,1,1,1);
}
