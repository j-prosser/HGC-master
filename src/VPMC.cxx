#define VPMC_cxx
#include "VPMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector> 
//using namespace std;



VPMC::VPMC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PU0_vito_code_2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PU0_vito_code_2.root");
      }
      f->GetObject("newC3Ds",tree);

   }
   Init(tree);
}

VPMC::~VPMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VPMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VPMC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void VPMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("endcap0", &endcap0_, &b_endcap0_);
   fChain->SetBranchAddress("endcap0.fUniqueID", endcap0_fUniqueID, &b_endcap0_fUniqueID);
   fChain->SetBranchAddress("endcap0.fBits", endcap0_fBits, &b_endcap0_fBits);
   fChain->SetBranchAddress("endcap0._id", endcap0__id, &b_endcap0__id);
   fChain->SetBranchAddress("endcap0._subdet", endcap0__subdet, &b_endcap0__subdet);
   fChain->SetBranchAddress("endcap0._pt", endcap0__pt, &b_endcap0__pt);
   fChain->SetBranchAddress("endcap0._energy", endcap0__energy, &b_endcap0__energy);
   fChain->SetBranchAddress("endcap0._eta", endcap0__eta, &b_endcap0__eta);
   fChain->SetBranchAddress("endcap0._phi", endcap0__phi, &b_endcap0__phi);
   fChain->SetBranchAddress("endcap0._x", endcap0__x, &b_endcap0__x);
   fChain->SetBranchAddress("endcap0._y", endcap0__y, &b_endcap0__y);
   fChain->SetBranchAddress("endcap0._z", endcap0__z, &b_endcap0__z);
   fChain->SetBranchAddress("endcap0._layer", endcap0__layer, &b_endcap0__layer);
   fChain->SetBranchAddress("endcap0._theta", endcap0__theta, &b_endcap0__theta);
   fChain->SetBranchAddress("endcap0._isTrigger", endcap0__isTrigger, &b_endcap0__isTrigger);
   fChain->SetBranchAddress("endcap0._firstLayer", endcap0__firstLayer, &b_endcap0__firstLayer);
   fChain->SetBranchAddress("endcap0._lastLayer", endcap0__lastLayer, &b_endcap0__lastLayer);
   fChain->SetBranchAddress("endcap0._maxLayer", endcap0__maxLayer, &b_endcap0__maxLayer);
   fChain->SetBranchAddress("endcap0._showerLength", endcap0__showerLength, &b_endcap0__showerLength);
   fChain->SetBranchAddress("endcap0._clusters", endcap0__clusters, &b_endcap0__clusters);
   fChain->SetBranchAddress("endcap0._cells", endcap0__cells, &b_endcap0__cells);
   //fChain->SetBranchAddress("endcap0._nearestGen", endcap0__nearestGen, &b_endcap0__nearestGen);
   fChain->SetBranchAddress("endcap0._xNorm", endcap0__xNorm, &b_endcap0__xNorm);
   fChain->SetBranchAddress("endcap0._yNorm", endcap0__yNorm, &b_endcap0__yNorm);
   fChain->SetBranchAddress("endcap1", &endcap1_, &b_endcap1_);
   fChain->SetBranchAddress("endcap1.fUniqueID", endcap1_fUniqueID, &b_endcap1_fUniqueID);
   fChain->SetBranchAddress("endcap1.fBits", endcap1_fBits, &b_endcap1_fBits);
   fChain->SetBranchAddress("endcap1._id", endcap1__id, &b_endcap1__id);
   fChain->SetBranchAddress("endcap1._subdet", endcap1__subdet, &b_endcap1__subdet);
   fChain->SetBranchAddress("endcap1._pt", endcap1__pt, &b_endcap1__pt);
   fChain->SetBranchAddress("endcap1._energy", endcap1__energy, &b_endcap1__energy);
   fChain->SetBranchAddress("endcap1._eta", endcap1__eta, &b_endcap1__eta);
   fChain->SetBranchAddress("endcap1._phi", endcap1__phi, &b_endcap1__phi);
   fChain->SetBranchAddress("endcap1._x", endcap1__x, &b_endcap1__x);
   fChain->SetBranchAddress("endcap1._y", endcap1__y, &b_endcap1__y);
   fChain->SetBranchAddress("endcap1._z", endcap1__z, &b_endcap1__z);
   fChain->SetBranchAddress("endcap1._layer", endcap1__layer, &b_endcap1__layer);
   fChain->SetBranchAddress("endcap1._theta", endcap1__theta, &b_endcap1__theta);
   fChain->SetBranchAddress("endcap1._isTrigger", endcap1__isTrigger, &b_endcap1__isTrigger);
   fChain->SetBranchAddress("endcap1._firstLayer", endcap1__firstLayer, &b_endcap1__firstLayer);
   fChain->SetBranchAddress("endcap1._lastLayer", endcap1__lastLayer, &b_endcap1__lastLayer);
   fChain->SetBranchAddress("endcap1._maxLayer", endcap1__maxLayer, &b_endcap1__maxLayer);
   fChain->SetBranchAddress("endcap1._showerLength", endcap1__showerLength, &b_endcap1__showerLength);
   fChain->SetBranchAddress("endcap1._clusters", endcap1__clusters, &b_endcap1__clusters);
   fChain->SetBranchAddress("endcap1._cells", endcap1__cells, &b_endcap1__cells);
   //fChain->SetBranchAddress("endcap1._nearestGen", endcap1__nearestGen, &b_endcap1__nearestGen);
   fChain->SetBranchAddress("endcap1._xNorm", endcap1__xNorm, &b_endcap1__xNorm);
   fChain->SetBranchAddress("endcap1._yNorm", endcap1__yNorm, &b_endcap1__yNorm);
   Notify();
}

Bool_t VPMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VPMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VPMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void VPMC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L VPMC.C
//      root> VPMC t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

	// UNSET all branches
   //fChain->SetBranchStatus("*",0);
   // Set Branches you want..
   //fChain->SetBranchStatus("",) 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	
	  //std::vector<double> XSX;

	std::cout << "Entry:\t" << jentry << "\n";	
	std::cout << "xNorm:\t" << endcap0__xNorm << "\n";
   }
}



