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
#include <iostream>
// Figure this out!

//test
//using namespace std;

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
   std::vector<unsigned int> endcap0__clusters[kMaxendcap0];
   std::vector<unsigned int> endcap0__cells[kMaxendcap0];
   //HGCgen          endcap0__nearestGen[kMaxendcap0];
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
   std::vector<unsigned int> endcap1__clusters[kMaxendcap1];
   std::vector<unsigned int> endcap1__cells[kMaxendcap1];
   //HGCgen          endcap1__nearestGen[kMaxendcap1];
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


