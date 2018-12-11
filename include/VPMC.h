//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 11 11:01:31 2018 by ROOT version 6.14/04
// from TTree newC3Ds/newC3Ds
// found on file: out.root
//////////////////////////////////////////////////////////

#ifndef VPMC_h
#define VPMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TObject.h"
#include "inc/HGChit.h"
#include "inc/HGCC3D.h"

class VPMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxendcap0 = 1;
   static constexpr Int_t kMaxendcap1 = 2;

   // Declaration of leaf types
   Int_t           endcap0_;
   UInt_t          endcap0_fUniqueID[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0_fBits[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0__id[kMaxendcap0];   //[endcap0_]
   Int_t           endcap0__subdet[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__pt[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__energy[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__eta[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__phi[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__x[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__y[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__z[kMaxendcap0];   //[endcap0_]
   Int_t           endcap0__layer[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__theta[kMaxendcap0];   //[endcap0_]
   Bool_t          endcap0__isTrigger[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0__firstLayer[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0__lastLayer[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0__maxLayer[kMaxendcap0];   //[endcap0_]
   UInt_t          endcap0__showerLength[kMaxendcap0];   //[endcap0_]
   vector<unsigned int> endcap0__clusters[kMaxendcap0];
   vector<unsigned int> endcap0__cells[kMaxendcap0];
   HGCgen          endcap0__nearestGen[kMaxendcap0];
   Float_t         endcap0__xNorm[kMaxendcap0];   //[endcap0_]
   Float_t         endcap0__yNorm[kMaxendcap0];   //[endcap0_]
   Int_t           endcap1_;
   UInt_t          endcap1_fUniqueID[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1_fBits[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1__id[kMaxendcap1];   //[endcap1_]
   Int_t           endcap1__subdet[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__pt[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__energy[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__eta[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__phi[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__x[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__y[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__z[kMaxendcap1];   //[endcap1_]
   Int_t           endcap1__layer[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__theta[kMaxendcap1];   //[endcap1_]
   Bool_t          endcap1__isTrigger[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1__firstLayer[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1__lastLayer[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1__maxLayer[kMaxendcap1];   //[endcap1_]
   UInt_t          endcap1__showerLength[kMaxendcap1];   //[endcap1_]
   vector<unsigned int> endcap1__clusters[kMaxendcap1];
   vector<unsigned int> endcap1__cells[kMaxendcap1];
   HGCgen          endcap1__nearestGen[kMaxendcap1];
   Float_t         endcap1__xNorm[kMaxendcap1];   //[endcap1_]
   Float_t         endcap1__yNorm[kMaxendcap1];   //[endcap1_]

   // List of branches
   TBranch        *b_endcap0_;   //!
   TBranch        *b_endcap0_fUniqueID;   //!
   TBranch        *b_endcap0_fBits;   //!
   TBranch        *b_endcap0__id;   //!
   TBranch        *b_endcap0__subdet;   //!
   TBranch        *b_endcap0__pt;   //!
   TBranch        *b_endcap0__energy;   //!
   TBranch        *b_endcap0__eta;   //!
   TBranch        *b_endcap0__phi;   //!
   TBranch        *b_endcap0__x;   //!
   TBranch        *b_endcap0__y;   //!
   TBranch        *b_endcap0__z;   //!
   TBranch        *b_endcap0__layer;   //!
   TBranch        *b_endcap0__theta;   //!
   TBranch        *b_endcap0__isTrigger;   //!
   TBranch        *b_endcap0__firstLayer;   //!
   TBranch        *b_endcap0__lastLayer;   //!
   TBranch        *b_endcap0__maxLayer;   //!
   TBranch        *b_endcap0__showerLength;   //!
   TBranch        *b_endcap0__clusters;   //!
   TBranch        *b_endcap0__cells;   //!
   TBranch        *b_endcap0__nearestGen;   //!
   TBranch        *b_endcap0__xNorm;   //!
   TBranch        *b_endcap0__yNorm;   //!
   TBranch        *b_endcap1_;   //!
   TBranch        *b_endcap1_fUniqueID;   //!
   TBranch        *b_endcap1_fBits;   //!
   TBranch        *b_endcap1__id;   //!
   TBranch        *b_endcap1__subdet;   //!
   TBranch        *b_endcap1__pt;   //!
   TBranch        *b_endcap1__energy;   //!
   TBranch        *b_endcap1__eta;   //!
   TBranch        *b_endcap1__phi;   //!
   TBranch        *b_endcap1__x;   //!
   TBranch        *b_endcap1__y;   //!
   TBranch        *b_endcap1__z;   //!
   TBranch        *b_endcap1__layer;   //!
   TBranch        *b_endcap1__theta;   //!
   TBranch        *b_endcap1__isTrigger;   //!
   TBranch        *b_endcap1__firstLayer;   //!
   TBranch        *b_endcap1__lastLayer;   //!
   TBranch        *b_endcap1__maxLayer;   //!
   TBranch        *b_endcap1__showerLength;   //!
   TBranch        *b_endcap1__clusters;   //!
   TBranch        *b_endcap1__cells;   //!
   TBranch        *b_endcap1__nearestGen;   //!
   TBranch        *b_endcap1__xNorm;   //!
   TBranch        *b_endcap1__yNorm;   //!

   VPMC(TTree *tree=0);
   virtual ~VPMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef VPMC_cxx
VPMC::VPMC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("out.root");
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
   fChain->SetBranchAddress("endcap0._nearestGen", endcap0__nearestGen, &b_endcap0__nearestGen);
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
   fChain->SetBranchAddress("endcap1._nearestGen", endcap1__nearestGen, &b_endcap1__nearestGen);
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
#endif // #ifdef VPMC_cxx
