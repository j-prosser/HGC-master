#include "TROOT.h"

void Loadzero() {
		TFile *f1 = new TFile("output/PU0/testout.root", "READ");
}

void Load200() {
		TFile *f2 = new TFile("output/PU200/testout.root", "READ");
}
void SetupGraph() {
		TCanvas *c1 = new TCanvas("c1", "multigraph", 700,500);
		c1->SetGrid();

		TMultiGraph *mg = new TMultiGraph();
}
void overlay() {
		gROOT->ProcessLine("TCanvas *c1 = new TCanvas(\"c1\", \"multigraph\", 700,500);");
		gROOT->ProcessLine("c1->SetGrid();");
	  	gROOT->ProcessLine("TMultiGraph *mg = new TMultiGraph();");
		Loadzero();
		gROOT->ProcessLine("mg->Add(Graph);");
		Load200();
		gROOT->ProcessLine("mg->Add(Graph);");
		gROOT->ProcessLine("mg->GetXaxis()->SetTitle(\"Radius (m) \")");
		gROOT->ProcessLine("mg->GetYaxis()->SetTitle(\"Sigma_E/E \")");
		gROOT->ProcessLine("mg->Draw(\"ac*\")");
}
