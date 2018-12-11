#include "TROOT.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "stdio.h"

TFile* Loadzero() {
		TFile *f1 = new TFile("output/PU0/testout.root", "READ");
		return f1;
}

TFile* Load200() {
		TFile *f2 = new TFile("output/PU200/testout.root", "READ");
		return f2;
}
TFile* Load( std::string filepath ) {
		TFile *_f = new TFile(filepath.c_str(), "READ");
		return _f;
}
TMultiGraph* SetupGraph() {
		TCanvas *c1 = new TCanvas("c1", "multigraph", 700,500);
		c1->SetGrid();

		TMultiGraph *mg = new TMultiGraph();
		return mg;
}
TGraph* getGraph(TFile *file, std::string name = "Graph") {
		TGraph *graph = (TGraph*) file->Get(name.c_str());
		return graph;
}
void setAxis(TMultiGraph *mg) {
		mg->GetXaxis()->SetTitle("Radius (m) ");
		mg->GetYaxis()->SetTitle("Sigma_E/E ");
}
TMultiGraph* overlay() {
		TMultiGraph *mg = SetupGraph();
		TFile *file_0 = Loadzero();
		TFile *file_200 = Load200();
		TGraph *Graph_0 = getGraph(file_0);
		TGraph *Graph_200 = getGraph(file_200);
		mg->Add(Graph_0);
		mg->Add(Graph_200);
		setAxis(mg);
		mg->Draw("ac*");
		return mg;
}
void plot() {
		TMultiGraph *mg_1 = SetupGraph();

		std::string filepath_0  = "output/PU0/output.root";
		std::string filepath_200  = "output/PU200/output.root";

		TFile *file_0 = new TFile(filepath_0.c_str(), "READ");
		std::string graphname_1 = "sigma_EE_r_1.47_1.98";
		std::string graphname_2 = "sigma_EE_r_1.98_2.49";
		std::string graphname_3 = "sigma_EE_r_2.49_3.00";
		std::string graphname_4 = "Sigma_EE_r";

		TGraph *Graph_0_1 = getGraph(file_0, graphname_1.c_str());
		TGraph *Graph_0_2 = getGraph(file_0, graphname_2.c_str());
		TGraph *Graph_0_3 = getGraph(file_0, graphname_3.c_str());
		TGraph *Graph_0_4 = getGraph(file_0, graphname_4.c_str());
		
		TFile *file_200 = new TFile(filepath_200.c_str(), "READ");

		TGraph *Graph_200_1 = getGraph(file_200, graphname_1.c_str());
		TGraph *Graph_200_2 = getGraph(file_200, graphname_2.c_str());
		TGraph *Graph_200_3 = getGraph(file_200, graphname_3.c_str());
		TGraph *Graph_200_4 = getGraph(file_200, graphname_4.c_str());

		mg_1->Add(Graph_0_1);
		mg_1->Add(Graph_0_2);
		mg_1->Add(Graph_0_3);
		mg_1->Add(Graph_0_4);
		mg_1->Add(Graph_200_1);
		mg_1->Add(Graph_200_2);
		mg_1->Add(Graph_200_3);
		mg_1->Add(Graph_200_4);
		mg_1->Draw("ac*");

}
void plot(int firstbin, int secondbin, int thirdbin, int total) {
		TMultiGraph *mg_1 = SetupGraph();

		std::string filepath_0  = "output/PU0/output.root";
		std::string filepath_200  = "output/PU200/output.root";

		TFile *file_0 = new TFile(filepath_0.c_str(), "READ");
		std::string graphname_1 = "sigma_EE_r_1.47_1.98";
		std::string graphname_2 = "sigma_EE_r_1.98_2.49";
		std::string graphname_3 = "sigma_EE_r_2.49_3.00";
		std::string graphname_4 = "Sigma_EE_r";

		TGraph *Graph_0_1 = getGraph(file_0, graphname_1.c_str());
		TGraph *Graph_0_2 = getGraph(file_0, graphname_2.c_str());
		TGraph *Graph_0_3 = getGraph(file_0, graphname_3.c_str());
		TGraph *Graph_0_4 = getGraph(file_0, graphname_4.c_str());
		
		TFile *file_200 = new TFile(filepath_200.c_str(), "READ");

		TGraph *Graph_200_1 = getGraph(file_200, graphname_1.c_str());
		TGraph *Graph_200_2 = getGraph(file_200, graphname_2.c_str());
		TGraph *Graph_200_3 = getGraph(file_200, graphname_3.c_str());
		TGraph *Graph_200_4 = getGraph(file_200, graphname_4.c_str());

		if (firstbin) {
			mg_1->Add(Graph_0_1);
			mg_1->Add(Graph_200_1);
		}
		if (secondbin) {
			mg_1->Add(Graph_0_2);
			mg_1->Add(Graph_200_2);
		}
		if (thirdbin) {
			mg_1->Add(Graph_0_3);
			mg_1->Add(Graph_200_3);
		}
		if (total) {
			mg_1->Add(Graph_0_4);
			mg_1->Add(Graph_200_4);
		}
		mg_1->Draw("ac*");
