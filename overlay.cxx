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
TMultiGraph* SetupGraph() {
		TCanvas *c1 = new TCanvas("c1", "multigraph", 700,500);
		c1->SetGrid();

		TMultiGraph *mg = new TMultiGraph();
		return mg;
}
TGraph* getGraph(TFile *file, std::string name = "Graph") {
		TGraph *graph = (TGraph*) file->Get("Graph");
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
