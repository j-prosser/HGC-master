#ifndef BuildTreeBase_h
#define BuildTreeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <math.h>
#include <iostream>
#include <vector>

class BuildTreeBase {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain


// Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Int_t           gen_n;
   Int_t           gen_PUNumInt;
   Float_t         gen_TrueNumInt;
   Float_t         vtx_x;
   Float_t         vtx_y;
   Float_t         vtx_z;
   std::vector<float>   *gen_eta;
   std::vector<float>   *gen_phi;
   std::vector<float>   *gen_pt;
   std::vector<float>   *gen_energy;
   std::vector<int>     *gen_charge;
   std::vector<int>     *gen_pdgid;
   std::vector<int>     *gen_status;
   std::vector<std::vector<int> > *gen_daughters;
   std::vector<float>   *genpart_eta;
   std::vector<float>   *genpart_phi;
   std::vector<float>   *genpart_pt;
   std::vector<float>   *genpart_energy;
   std::vector<float>   *genpart_dvx;
   std::vector<float>   *genpart_dvy;
   std::vector<float>   *genpart_dvz;
   std::vector<float>   *genpart_ovx;
   std::vector<float>   *genpart_ovy;
   std::vector<float>   *genpart_ovz;
   std::vector<int>     *genpart_mother;
   std::vector<float>   *genpart_exphi;
   std::vector<float>   *genpart_exeta;
   std::vector<float>   *genpart_exx;
   std::vector<float>   *genpart_exy;
   std::vector<float>   *genpart_fbrem;
   std::vector<int>     *genpart_pid;
   std::vector<int>     *genpart_gen;
   std::vector<int>     *genpart_reachedEE;
   std::vector<bool>    *genpart_fromBeamPipe;
   std::vector<std::vector<float> > *genpart_posx;
   std::vector<std::vector<float> > *genpart_posy;
   std::vector<std::vector<float> > *genpart_posz;
   Int_t           genjet_n;
   std::vector<float>   *genjet_energy;
   std::vector<float>   *genjet_pt;
   std::vector<float>   *genjet_eta;
   std::vector<float>   *genjet_phi;
   std::vector<float>   *gentau_pt;
   std::vector<float>   *gentau_eta;
   std::vector<float>   *gentau_phi;
   std::vector<float>   *gentau_energy;
   std::vector<float>   *gentau_mass;
   std::vector<float>   *gentau_vis_pt;
   std::vector<float>   *gentau_vis_eta;
   std::vector<float>   *gentau_vis_phi;
   std::vector<float>   *gentau_vis_energy;
   std::vector<float>   *gentau_vis_mass;
   std::vector<std::vector<float> > *gentau_products_pt;
   std::vector<std::vector<float> > *gentau_products_eta;
   std::vector<std::vector<float> > *gentau_products_phi;
   std::vector<std::vector<float> > *gentau_products_energy;
   std::vector<std::vector<float> > *gentau_products_mass;
   std::vector<std::vector<int> > *gentau_products_id;
   std::vector<int>     *gentau_decayMode;
   std::vector<int>     *gentau_totNproducts;
   std::vector<int>     *gentau_totNgamma;
   std::vector<int>     *gentau_totNpiZero;
   std::vector<int>     *gentau_totNcharged;
   Int_t           hgcdigi_n;
   std::vector<int>     *hgcdigi_id;
   std::vector<int>     *hgcdigi_subdet;
   std::vector<int>     *hgcdigi_zside;
   std::vector<int>     *hgcdigi_layer;
   std::vector<int>     *hgcdigi_wafer;
   std::vector<int>     *hgcdigi_wafertype;
   std::vector<int>     *hgcdigi_cell;
   std::vector<float>   *hgcdigi_eta;
   std::vector<float>   *hgcdigi_phi;
   std::vector<float>   *hgcdigi_z;
   std::vector<unsigned int> *hgcdigi_data;
   std::vector<int>     *hgcdigi_isadc;
   Int_t           bhdigi_n;
   std::vector<int>     *bhdigi_id;
   std::vector<int>     *bhdigi_subdet;
   std::vector<int>     *bhdigi_zside;
   std::vector<int>     *bhdigi_layer;
   std::vector<int>     *bhdigi_ieta;
   std::vector<int>     *bhdigi_iphi;
   std::vector<float>   *bhdigi_eta;
   std::vector<float>   *bhdigi_phi;
   std::vector<float>   *bhdigi_z;
   std::vector<unsigned int> *bhdigi_data;
   Int_t           tc_n;
   std::vector<unsigned int> *tc_id;
   std::vector<int>     *tc_subdet;
   std::vector<int>     *tc_zside;
   std::vector<int>     *tc_layer;
   std::vector<int>     *tc_wafer;
   std::vector<int>     *tc_wafertype;
   std::vector<int>     *tc_cell;
   std::vector<unsigned int> *tc_data;
   std::vector<float>   *tc_pt;
   std::vector<float>   *tc_mipPt;
   std::vector<float>   *tc_energy;
   std::vector<float>   *tc_eta;
   std::vector<float>   *tc_phi;
   std::vector<float>   *tc_x;
   std::vector<float>   *tc_y;
   std::vector<float>   *tc_z;
   std::vector<unsigned int> *tc_cluster_id;
   std::vector<unsigned int> *tc_multicluster_id;
   std::vector<float>   *tc_multicluster_pt;
   std::vector<int>     *tc_super_n;
   std::vector<double>  *tc_ave_dist;
   Int_t           cl_n;
   std::vector<unsigned int> *cl_id;
   std::vector<float>   *cl_mipPt;
   std::vector<float>   *cl_pt;
   std::vector<float>   *cl_energy;
   std::vector<float>   *cl_eta;
   std::vector<float>   *cl_phi;
   std::vector<int>     *cl_layer;
   std::vector<int>     *cl_subdet;
   std::vector<int>     *cl_cells_n;
   std::vector<std::vector<unsigned int> > *cl_cells_id;
   std::vector<unsigned int> *cl_multicluster_id;
   std::vector<float>   *cl_multicluster_pt;
   Int_t           cl3d_n;
   std::vector<unsigned int> *cl3d_id;
   std::vector<float>   *cl3d_pt;
   std::vector<float>   *cl3d_energy;
   std::vector<float>   *cl3d_eta;
   std::vector<float>   *cl3d_phi;
   std::vector<int>     *cl3d_clusters_n;
   std::vector<std::vector<unsigned int> > *cl3d_clusters_id;
   std::vector<int>     *cl3d_showerlength;
   std::vector<int>     *cl3d_coreshowerlength;
   std::vector<int>     *cl3d_firstlayer;
   std::vector<int>     *cl3d_maxlayer;
   std::vector<float>   *cl3d_seetot;
   std::vector<float>   *cl3d_seemax;
   std::vector<float>   *cl3d_spptot;
   std::vector<float>   *cl3d_sppmax;
   std::vector<float>   *cl3d_szz;
   std::vector<float>   *cl3d_srrtot;
   std::vector<float>   *cl3d_srrmax;
   std::vector<float>   *cl3d_srrmean;
   std::vector<float>   *cl3d_emaxe;
   std::vector<float>   *cl3d_bdteg;
   std::vector<int>     *cl3d_quality;
   Int_t           tower_n;
   std::vector<float>   *tower_pt;
   std::vector<float>   *tower_energy;
   std::vector<float>   *tower_eta;
   std::vector<float>   *tower_phi;
   std::vector<float>   *tower_etEm;
   std::vector<float>   *tower_etHad;
   std::vector<int>     *tower_iEta;
   std::vector<int>     *tower_iPhi;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_gen_n;   //!
   TBranch        *b_gen_PUNumInt;   //!
   TBranch        *b_gen_TrueNumInt;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_energy;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_pdgid;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_daughters;   //!
   TBranch        *b_genpart_eta;   //!
   TBranch        *b_genpart_phi;   //!
   TBranch        *b_genpart_pt;   //!
   TBranch        *b_genpart_energy;   //!
   TBranch        *b_genpart_dvx;   //!
   TBranch        *b_genpart_dvy;   //!
   TBranch        *b_genpart_dvz;   //!
   TBranch        *b_genpart_ovx;   //!
   TBranch        *b_genpart_ovy;   //!
   TBranch        *b_genpart_ovz;   //!
   TBranch        *b_genpart_mother;   //!
   TBranch        *b_genpart_exphi;   //!
   TBranch        *b_genpart_exeta;   //!
   TBranch        *b_genpart_exx;   //!
   TBranch        *b_genpart_exy;   //!
   TBranch        *b_genpart_fbrem;   //!
   TBranch        *b_genpart_pid;   //!
   TBranch        *b_genpart_gen;   //!
   TBranch        *b_genpart_reachedEE;   //!
   TBranch        *b_genpart_fromBeamPipe;   //!
   TBranch        *b_genpart_posx;   //!
   TBranch        *b_genpart_posy;   //!
   TBranch        *b_genpart_posz;   //!
   TBranch        *b_genjet_n;   //!
   TBranch        *b_genjet_energy;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_gentau_pt;   //!
   TBranch        *b_gentau_eta;   //!
   TBranch        *b_gentau_phi;   //!
   TBranch        *b_gentau_energy;   //!
   TBranch        *b_gentau_mass;   //!
   TBranch        *b_gentau_vis_pt;   //!
   TBranch        *b_gentau_vis_eta;   //!
   TBranch        *b_gentau_vis_phi;   //!
   TBranch        *b_gentau_vis_energy;   //!
   TBranch        *b_gentau_vis_mass;   //!
   TBranch        *b_gentau_products_pt;   //!
   TBranch        *b_gentau_products_eta;   //!
   TBranch        *b_gentau_products_phi;   //!
   TBranch        *b_gentau_products_energy;   //!
   TBranch        *b_gentau_products_mass;   //!
   TBranch        *b_gentau_products_id;   //!
   TBranch        *b_gentau_decayMode;   //!
   TBranch        *b_gentau_totNproducts;   //!
   TBranch        *b_gentau_totNgamma;   //!
   TBranch        *b_gentau_totNpiZero;   //!
   TBranch        *b_gentau_totNcharged;   //!
   TBranch        *b_hgcdigi_n;   //!
   TBranch        *b_hgcdigi_id;   //!
   TBranch        *b_hgcdigi_subdet;   //!
   TBranch        *b_hgcdigi_zside;   //!
   TBranch        *b_hgcdigi_layer;   //!
   TBranch        *b_hgcdigi_wafer;   //!
   TBranch        *b_hgcdigi_wafertype;   //!
   TBranch        *b_hgcdigi_cell;   //!
   TBranch        *b_hgcdigi_eta;   //!
   TBranch        *b_hgcdigi_phi;   //!
   TBranch        *b_hgcdigi_z;   //!
   TBranch        *b_hgcdigi_data;   //!
   TBranch        *b_hgcdigi_isadc;   //!
   TBranch        *b_bhdigi_n;   //!
   TBranch        *b_bhdigi_id;   //!
   TBranch        *b_bhdigi_subdet;   //!
   TBranch        *b_bhdigi_zside;   //!
   TBranch        *b_bhdigi_layer;   //!
   TBranch        *b_bhdigi_ieta;   //!
   TBranch        *b_bhdigi_iphi;   //!
   TBranch        *b_bhdigi_eta;   //!
   TBranch        *b_bhdigi_phi;   //!
   TBranch        *b_bhdigi_z;   //!
   TBranch        *b_bhdigi_data;   //!
   TBranch        *b_tc_n;   //!
   TBranch        *b_tc_id;   //!
   TBranch        *b_tc_subdet;   //!
   TBranch        *b_tc_zside;   //!
   TBranch        *b_tc_layer;   //!
   TBranch        *b_tc_wafer;   //!
   TBranch        *b_tc_wafertype;   //!
   TBranch        *b_tc_cell;   //!
   TBranch        *b_tc_data;   //!
   TBranch        *b_tc_pt;   //!
   TBranch        *b_tc_mipPt;   //!
   TBranch        *b_tc_energy;   //!
   TBranch        *b_tc_eta;   //!
   TBranch        *b_tc_phi;   //!
   TBranch        *b_tc_x;   //!
   TBranch        *b_tc_y;   //!
   TBranch        *b_tc_z;   //!
   TBranch        *b_tc_cluster_id;   //!
   TBranch        *b_tc_multicluster_id;   //!
   TBranch        *b_tc_multicluster_pt;   //!
   TBranch        *b_tc_super_n;   //!
   TBranch        *b_tc_ave_dist;   //!
   TBranch        *b_cl_n;   //!
   TBranch        *b_cl_id;   //!
   TBranch        *b_cl_mipPt;   //!
   TBranch        *b_cl_pt;   //!
   TBranch        *b_cl_energy;   //!
   TBranch        *b_cl_eta;   //!
   TBranch        *b_cl_phi;   //!
   TBranch        *b_cl_layer;   //!
   TBranch        *b_cl_subdet;   //!
   TBranch        *b_cl_cells_n;   //!
   TBranch        *b_cl_cells_id;   //!
   TBranch        *b_cl_multicluster_id;   //!
   TBranch        *b_cl_multicluster_pt;   //!
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_id;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_clusters_n;   //!
   TBranch        *b_cl3d_clusters_id;   //!
   TBranch        *b_cl3d_showerlength;   //!
   TBranch        *b_cl3d_coreshowerlength;   //!
   TBranch        *b_cl3d_firstlayer;   //!
   TBranch        *b_cl3d_maxlayer;   //!
   TBranch        *b_cl3d_seetot;   //!
   TBranch        *b_cl3d_seemax;   //!
   TBranch        *b_cl3d_spptot;   //!
   TBranch        *b_cl3d_sppmax;   //!
   TBranch        *b_cl3d_szz;   //!
   TBranch        *b_cl3d_srrtot;   //!
   TBranch        *b_cl3d_srrmax;   //!
   TBranch        *b_cl3d_srrmean;   //!
   TBranch        *b_cl3d_emaxe;   //!
   TBranch        *b_cl3d_bdteg;   //!
   TBranch        *b_cl3d_quality;   //!
   TBranch        *b_tower_n;   //!
   TBranch        *b_tower_pt;   //!
   TBranch        *b_tower_energy;   //!
   TBranch        *b_tower_eta;   //!
   TBranch        *b_tower_phi;   //!
   TBranch        *b_tower_etEm;   //!
   TBranch        *b_tower_etHad;   //!
   TBranch        *b_tower_iEta;   //!
   TBranch        *b_tower_iPhi;   //!

  BuildTreeBase(TTree *tree=0);
   virtual ~BuildTreeBase();
   virtual Int_t    Cut();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BuildTreeBase_cxx
BuildTreeBase::BuildTreeBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {

      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/vols/cms/snwebb/HGC_ntuples/SingleGammaPt25_PU0_threshold-2/ntuples/ntuple_1.root");
      //if (!f || !f->IsOpen()) {
         //f = new TFile("/vols/cms/snwebb/HGC_ntuples/SingleGammaPt25_PU0_threshold-2/ntuples/ntuple_1.root");
      //}
      //TDirectory * dir = (TDirectory*)f->Get("/vols/cms/snwebb/HGC_ntuples/SingleGammaPt25_PU0_threshold-2/ntuples/ntuple_1.root:/hgcalTriggerNtuplizer");
      //dir->GetObject("HGCalTriggerNtuple",tree);
      //TFile *f (TFile*)gROOT->GetListOfFiles()->FindObject("ntuples/ntuple_1.root");
      
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/johannes/host/simple-data/ntuples/ntuple_1.root"); ///home/johannes/host/simple-data/ntuples/ntuple_1.root");
      
	  //std::cout << "IS THIS RUN?  " << std::endl;
	  
	  //if (!f || !f->IsOpen()) {
        //f = new TFile("/home/johannes/host/simple-data/ntuples/ntuple_1.root");
      //}
	  ///home/johannes/host/simple-data/ntuples/ntuple_1.root:/hgcalTriggerNtuplizer
      //TDirectory * dir = (TDirectory*)f->Get("/home/johannes/host/simple-data/ntuples/ntuple_1.root:/hgcalTriggerNtuplizer");
      //dir->GetObject("HGCalTriggerNtuple",tree);

   //}
   //Init(tree);
}

BuildTreeBase::~BuildTreeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BuildTreeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BuildTreeBase::LoadTree(Long64_t entry)
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

void BuildTreeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   gen_eta = 0;
   gen_phi = 0;
   gen_pt = 0;
   gen_energy = 0;
   gen_charge = 0;
   gen_pdgid = 0;
   gen_status = 0;
   gen_daughters = 0;
   genpart_eta = 0;
   genpart_phi = 0;
   genpart_pt = 0;
   genpart_energy = 0;
   genpart_dvx = 0;
   genpart_dvy = 0;
   genpart_dvz = 0;
   genpart_ovx = 0;
   genpart_ovy = 0;
   genpart_ovz = 0;
   genpart_mother = 0;
   genpart_exphi = 0;
   genpart_exeta = 0;
   genpart_exx = 0;
   genpart_exy = 0;
   genpart_fbrem = 0;
   genpart_pid = 0;
   genpart_gen = 0;
   genpart_reachedEE = 0;
   genpart_fromBeamPipe = 0;
   genpart_posx = 0;
   genpart_posy = 0;
   genpart_posz = 0;
   genjet_energy = 0;
   genjet_pt = 0;
   genjet_eta = 0;
   genjet_phi = 0;
   gentau_pt = 0;
   gentau_eta = 0;
   gentau_phi = 0;
   gentau_energy = 0;
   gentau_mass = 0;
   gentau_vis_pt = 0;
   gentau_vis_eta = 0;
   gentau_vis_phi = 0;
   gentau_vis_energy = 0;
   gentau_vis_mass = 0;
   gentau_products_pt = 0;
   gentau_products_eta = 0;
   gentau_products_phi = 0;
   gentau_products_energy = 0;
   gentau_products_mass = 0;
   gentau_products_id = 0;
   gentau_decayMode = 0;
   gentau_totNproducts = 0;
   gentau_totNgamma = 0;
   gentau_totNpiZero = 0;
   gentau_totNcharged = 0;
   hgcdigi_id = 0;
   hgcdigi_subdet = 0;
   hgcdigi_zside = 0;
   hgcdigi_layer = 0;
   hgcdigi_wafer = 0;
   hgcdigi_wafertype = 0;
   hgcdigi_cell = 0;
   hgcdigi_eta = 0;
   hgcdigi_phi = 0;
   hgcdigi_z = 0;
   hgcdigi_data = 0;
   hgcdigi_isadc = 0;
   bhdigi_id = 0;
   bhdigi_subdet = 0;
   bhdigi_zside = 0;
   bhdigi_layer = 0;
   bhdigi_ieta = 0;
   bhdigi_iphi = 0;
   bhdigi_eta = 0;
   bhdigi_phi = 0;
   bhdigi_z = 0;
   bhdigi_data = 0;
   tc_id = 0;
   tc_subdet = 0;
   tc_zside = 0;
   tc_layer = 0;
   tc_wafer = 0;
   tc_wafertype = 0;
   tc_cell = 0;
   tc_data = 0;
   tc_pt = 0;
   tc_mipPt = 0;
   tc_energy = 0;
   tc_eta = 0;
   tc_phi = 0;
   tc_x = 0;
   tc_y = 0;
   tc_z = 0;
   tc_cluster_id = 0;
   tc_multicluster_id = 0;
   tc_multicluster_pt = 0;
   tc_super_n = 0;
   tc_ave_dist = 0;
   cl_id = 0;
   cl_mipPt = 0;
   cl_pt = 0;
   cl_energy = 0;
   cl_eta = 0;
   cl_phi = 0;
   cl_layer = 0;
   cl_subdet = 0;
   cl_cells_n = 0;
   cl_cells_id = 0;
   cl_multicluster_id = 0;
   cl_multicluster_pt = 0;
   cl3d_id = 0;
   cl3d_pt = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_clusters_n = 0;
   cl3d_clusters_id = 0;
   cl3d_showerlength = 0;
   cl3d_coreshowerlength = 0;
   cl3d_firstlayer = 0;
   cl3d_maxlayer = 0;
   cl3d_seetot = 0;
   cl3d_seemax = 0;
   cl3d_spptot = 0;
   cl3d_sppmax = 0;
   cl3d_szz = 0;
   cl3d_srrtot = 0;
   cl3d_srrmax = 0;
   cl3d_srrmean = 0;
   cl3d_emaxe = 0;
   cl3d_bdteg = 0;
   cl3d_quality = 0;
   tower_pt = 0;
   tower_energy = 0;
   tower_eta = 0;
   tower_phi = 0;
   tower_etEm = 0;
   tower_etHad = 0;
   tower_iEta = 0;
   tower_iPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

 if (fChain->GetBranch( "run"				)!=0 )     fChain->SetBranchAddress("run"			, &run, &b_run);
 if (fChain->GetBranch( "event"				)!=0 )     fChain->SetBranchAddress("event"			, &event, &b_event);
 if (fChain->GetBranch( "lumi"				)!=0 )     fChain->SetBranchAddress("lumi"			, &lumi, &b_lumi);
 if (fChain->GetBranch( "gen_n"				)!=0 )     fChain->SetBranchAddress("gen_n"			, &gen_n, &b_gen_n);
 if (fChain->GetBranch( "gen_PUNumInt"			)!=0 )     fChain->SetBranchAddress("gen_PUNumInt"		, &gen_PUNumInt, &b_gen_PUNumInt);
 if (fChain->GetBranch( "gen_TrueNumInt"		)!=0 )     fChain->SetBranchAddress("gen_TrueNumInt"		, &gen_TrueNumInt, &b_gen_TrueNumInt);
 if (fChain->GetBranch( "vtx_x"				)!=0 )     fChain->SetBranchAddress("vtx_x"			, &vtx_x, &b_vtx_x);
 if (fChain->GetBranch( "vtx_y"				)!=0 )     fChain->SetBranchAddress("vtx_y"			, &vtx_y, &b_vtx_y);
 if (fChain->GetBranch( "vtx_z"				)!=0 )     fChain->SetBranchAddress("vtx_z"			, &vtx_z, &b_vtx_z);
 if (fChain->GetBranch( "gen_eta"			)!=0 )     fChain->SetBranchAddress("gen_eta"			, &gen_eta, &b_gen_eta);
 if (fChain->GetBranch( "gen_phi"			)!=0 )     fChain->SetBranchAddress("gen_phi"			, &gen_phi, &b_gen_phi);
 if (fChain->GetBranch( "gen_pt"			)!=0 )     fChain->SetBranchAddress("gen_pt"			, &gen_pt, &b_gen_pt);
 if (fChain->GetBranch( "gen_energy"			)!=0 )     fChain->SetBranchAddress("gen_energy" 		, &gen_energy, &b_gen_energy);
 if (fChain->GetBranch( "gen_charge"			)!=0 )     fChain->SetBranchAddress("gen_charge"         	, &gen_charge, &b_gen_charge);
 if (fChain->GetBranch( "gen_pdgid"			)!=0 )     fChain->SetBranchAddress("gen_pdgid"			, &gen_pdgid, &b_gen_pdgid);
 if (fChain->GetBranch( "gen_status"			)!=0 )     fChain->SetBranchAddress("gen_status" 		, &gen_status, &b_gen_status);
 if (fChain->GetBranch( "gen_daughters"			)!=0 )     fChain->SetBranchAddress("gen_daughters"		, &gen_daughters, &b_gen_daughters);
 if (fChain->GetBranch( "genpart_eta"			)!=0 )     fChain->SetBranchAddress("genpart_eta"		, &genpart_eta, &b_genpart_eta);
 if (fChain->GetBranch( "genpart_phi"			)!=0 )     fChain->SetBranchAddress("genpart_phi"		, &genpart_phi, &b_genpart_phi);
 if (fChain->GetBranch( "genpart_pt"			)!=0 )     fChain->SetBranchAddress("genpart_pt" 		, &genpart_pt, &b_genpart_pt);
 if (fChain->GetBranch( "genpart_energy"		)!=0 )     fChain->SetBranchAddress("genpart_energy"		, &genpart_energy, &b_genpart_energy);
 if (fChain->GetBranch( "genpart_dvx"			)!=0 )     fChain->SetBranchAddress("genpart_dvx"		, &genpart_dvx, &b_genpart_dvx);
 if (fChain->GetBranch( "genpart_dvy"			)!=0 )     fChain->SetBranchAddress("genpart_dvy"		, &genpart_dvy, &b_genpart_dvy);
 if (fChain->GetBranch( "genpart_dvz"			)!=0 )     fChain->SetBranchAddress("genpart_dvz"		, &genpart_dvz, &b_genpart_dvz);
 if (fChain->GetBranch( "genpart_ovx"			)!=0 )     fChain->SetBranchAddress("genpart_ovx"		, &genpart_ovx, &b_genpart_ovx);
 if (fChain->GetBranch( "genpart_ovy"			)!=0 )     fChain->SetBranchAddress("genpart_ovy"		, &genpart_ovy, &b_genpart_ovy);
 if (fChain->GetBranch( "genpart_ovz"			)!=0 )     fChain->SetBranchAddress("genpart_ovz"		, &genpart_ovz, &b_genpart_ovz);
 if (fChain->GetBranch( "genpart_mother"		)!=0 )     fChain->SetBranchAddress("genpart_mother"		, &genpart_mother, &b_genpart_mother);
 if (fChain->GetBranch( "genpart_exphi"			)!=0 )     fChain->SetBranchAddress("genpart_exphi"		, &genpart_exphi, &b_genpart_exphi);
 if (fChain->GetBranch( "genpart_exeta"			)!=0 )     fChain->SetBranchAddress("genpart_exeta"		, &genpart_exeta, &b_genpart_exeta);
 if (fChain->GetBranch( "genpart_exx"			)!=0 )     fChain->SetBranchAddress("genpart_exx"		, &genpart_exx, &b_genpart_exx);
 if (fChain->GetBranch( "genpart_exy"			)!=0 )     fChain->SetBranchAddress("genpart_exy"		, &genpart_exy, &b_genpart_exy);
 if (fChain->GetBranch( "genpart_fbrem"			)!=0 )     fChain->SetBranchAddress("genpart_fbrem"		, &genpart_fbrem, &b_genpart_fbrem);
 if (fChain->GetBranch( "genpart_pid"			)!=0 )     fChain->SetBranchAddress("genpart_pid"		, &genpart_pid, &b_genpart_pid);
 if (fChain->GetBranch( "genpart_gen"			)!=0 )     fChain->SetBranchAddress("genpart_gen"		, &genpart_gen, &b_genpart_gen);
 if (fChain->GetBranch( "genpart_reachedEE"		)!=0 )     fChain->SetBranchAddress("genpart_reachedEE"		, &genpart_reachedEE, &b_genpart_reachedEE);
 if (fChain->GetBranch( "genpart_fromBeamPipe"		)!=0 )     fChain->SetBranchAddress("genpart_fromBeamPipe"	, &genpart_fromBeamPipe, &b_genpart_fromBeamPipe);
 if (fChain->GetBranch( "genpart_posx"			)!=0 )     fChain->SetBranchAddress("genpart_posx"		, &genpart_posx, &b_genpart_posx);
 if (fChain->GetBranch( "genpart_posy"			)!=0 )     fChain->SetBranchAddress("genpart_posy"		, &genpart_posy, &b_genpart_posy);
 if (fChain->GetBranch( "genpart_posz"			)!=0 )     fChain->SetBranchAddress("genpart_posz"		, &genpart_posz, &b_genpart_posz);
 if (fChain->GetBranch( "genjet_n"			)!=0 )     fChain->SetBranchAddress("genjet_n"			, &genjet_n, &b_genjet_n);
 if (fChain->GetBranch( "genjet_energy"			)!=0 )     fChain->SetBranchAddress("genjet_energy"		, &genjet_energy, &b_genjet_energy);
 if (fChain->GetBranch( "genjet_pt"			)!=0 )     fChain->SetBranchAddress("genjet_pt"			, &genjet_pt, &b_genjet_pt);
 if (fChain->GetBranch( "genjet_eta"			)!=0 )     fChain->SetBranchAddress("genjet_eta" 		, &genjet_eta, &b_genjet_eta);
 if (fChain->GetBranch( "genjet_phi"			)!=0 )     fChain->SetBranchAddress("genjet_phi" 		, &genjet_phi, &b_genjet_phi);
 if (fChain->GetBranch( "gentau_pt"			)!=0 )     fChain->SetBranchAddress("gentau_pt"			, &gentau_pt, &b_gentau_pt);
 if (fChain->GetBranch( "gentau_eta"			)!=0 )     fChain->SetBranchAddress("gentau_eta" 		, &gentau_eta, &b_gentau_eta);
 if (fChain->GetBranch( "gentau_phi"			)!=0 )     fChain->SetBranchAddress("gentau_phi" 		, &gentau_phi, &b_gentau_phi);
 if (fChain->GetBranch( "gentau_energy"			)!=0 )     fChain->SetBranchAddress("gentau_energy"		, &gentau_energy, &b_gentau_energy);
 if (fChain->GetBranch( "gentau_mass"			)!=0 )     fChain->SetBranchAddress("gentau_mass"		, &gentau_mass, &b_gentau_mass);
 if (fChain->GetBranch( "gentau_vis_pt"			)!=0 )     fChain->SetBranchAddress("gentau_vis_pt"		, &gentau_vis_pt, &b_gentau_vis_pt);
 if (fChain->GetBranch( "gentau_vis_eta"		)!=0 )     fChain->SetBranchAddress("gentau_vis_eta"		, &gentau_vis_eta, &b_gentau_vis_eta);
 if (fChain->GetBranch( "gentau_vis_phi"		)!=0 )     fChain->SetBranchAddress("gentau_vis_phi"		, &gentau_vis_phi, &b_gentau_vis_phi);
 if (fChain->GetBranch( "gentau_vis_energy"		)!=0 )     fChain->SetBranchAddress("gentau_vis_energy"		, &gentau_vis_energy, &b_gentau_vis_energy);
 if (fChain->GetBranch( "gentau_vis_mass"		)!=0 )     fChain->SetBranchAddress("gentau_vis_mass"		, &gentau_vis_mass, &b_gentau_vis_mass);
 if (fChain->GetBranch( "gentau_products_pt"		)!=0 )     fChain->SetBranchAddress("gentau_products_pt" 	, &gentau_products_pt, &b_gentau_products_pt);
 if (fChain->GetBranch( "gentau_products_eta"		)!=0 )     fChain->SetBranchAddress("gentau_products_eta"	, &gentau_products_eta, &b_gentau_products_eta);
 if (fChain->GetBranch( "gentau_products_phi"		)!=0 )     fChain->SetBranchAddress("gentau_products_phi"	, &gentau_products_phi, &b_gentau_products_phi);
 if (fChain->GetBranch( "gentau_products_energy"	)!=0 )     fChain->SetBranchAddress("gentau_products_energy"	, &gentau_products_energy, &b_gentau_products_energy);
 if (fChain->GetBranch( "gentau_products_mass"		)!=0 )     fChain->SetBranchAddress("gentau_products_mass"	, &gentau_products_mass, &b_gentau_products_mass);
 if (fChain->GetBranch( "gentau_products_id"		)!=0 )     fChain->SetBranchAddress("gentau_products_id" 	, &gentau_products_id, &b_gentau_products_id);
 if (fChain->GetBranch( "gentau_decayMode"		)!=0 )     fChain->SetBranchAddress("gentau_decayMode"		, &gentau_decayMode, &b_gentau_decayMode);
 if (fChain->GetBranch( "gentau_totNproducts"		)!=0 )     fChain->SetBranchAddress("gentau_totNproducts"	, &gentau_totNproducts, &b_gentau_totNproducts);
 if (fChain->GetBranch( "gentau_totNgamma"		)!=0 )     fChain->SetBranchAddress("gentau_totNgamma"		, &gentau_totNgamma, &b_gentau_totNgamma);
 if (fChain->GetBranch( "gentau_totNpiZero"		)!=0 )     fChain->SetBranchAddress("gentau_totNpiZero"		, &gentau_totNpiZero, &b_gentau_totNpiZero);
 if (fChain->GetBranch( "gentau_totNcharged"		)!=0 )     fChain->SetBranchAddress("gentau_totNcharged" 	, &gentau_totNcharged, &b_gentau_totNcharged);
 if (fChain->GetBranch( "hgcdigi_n"			)!=0 )     fChain->SetBranchAddress("hgcdigi_n"			, &hgcdigi_n, &b_hgcdigi_n);
 if (fChain->GetBranch( "hgcdigi_id"			)!=0 )     fChain->SetBranchAddress("hgcdigi_id" 		, &hgcdigi_id, &b_hgcdigi_id);
 if (fChain->GetBranch( "hgcdigi_subdet"		)!=0 )     fChain->SetBranchAddress("hgcdigi_subdet"		, &hgcdigi_subdet, &b_hgcdigi_subdet);
 if (fChain->GetBranch( "hgcdigi_zside"			)!=0 )     fChain->SetBranchAddress("hgcdigi_zside"		, &hgcdigi_zside, &b_hgcdigi_zside);
 if (fChain->GetBranch( "hgcdigi_layer"			)!=0 )     fChain->SetBranchAddress("hgcdigi_layer"		, &hgcdigi_layer, &b_hgcdigi_layer);
 if (fChain->GetBranch( "hgcdigi_wafer"			)!=0 )     fChain->SetBranchAddress("hgcdigi_wafer"		, &hgcdigi_wafer, &b_hgcdigi_wafer);
 if (fChain->GetBranch( "hgcdigi_wafertype"		)!=0 )     fChain->SetBranchAddress("hgcdigi_wafertype"		, &hgcdigi_wafertype, &b_hgcdigi_wafertype);
 if (fChain->GetBranch( "hgcdigi_cell"			)!=0 )     fChain->SetBranchAddress("hgcdigi_cell"		, &hgcdigi_cell, &b_hgcdigi_cell);
 if (fChain->GetBranch( "hgcdigi_eta"			)!=0 )     fChain->SetBranchAddress("hgcdigi_eta"		, &hgcdigi_eta, &b_hgcdigi_eta);
 if (fChain->GetBranch( "hgcdigi_phi"			)!=0 )     fChain->SetBranchAddress("hgcdigi_phi"		, &hgcdigi_phi, &b_hgcdigi_phi);
 if (fChain->GetBranch( "hgcdigi_z"			)!=0 )     fChain->SetBranchAddress("hgcdigi_z"			, &hgcdigi_z, &b_hgcdigi_z);
 if (fChain->GetBranch( "hgcdigi_data"			)!=0 )     fChain->SetBranchAddress("hgcdigi_data"		, &hgcdigi_data, &b_hgcdigi_data);
 if (fChain->GetBranch( "hgcdigi_isadc"			)!=0 )     fChain->SetBranchAddress("hgcdigi_isadc"		, &hgcdigi_isadc, &b_hgcdigi_isadc);
 if (fChain->GetBranch( "bhdigi_n"			)!=0 )     fChain->SetBranchAddress("bhdigi_n"			, &bhdigi_n, &b_bhdigi_n);
 if (fChain->GetBranch( "bhdigi_id"			)!=0 )     fChain->SetBranchAddress("bhdigi_id"			, &bhdigi_id, &b_bhdigi_id);
 if (fChain->GetBranch( "bhdigi_subdet"			)!=0 )     fChain->SetBranchAddress("bhdigi_subdet"		, &bhdigi_subdet, &b_bhdigi_subdet);
 if (fChain->GetBranch( "bhdigi_zside"			)!=0 )     fChain->SetBranchAddress("bhdigi_zside"		, &bhdigi_zside, &b_bhdigi_zside);
 if (fChain->GetBranch( "bhdigi_layer"			)!=0 )     fChain->SetBranchAddress("bhdigi_layer"		, &bhdigi_layer, &b_bhdigi_layer);
 if (fChain->GetBranch( "bhdigi_ieta"			)!=0 )     fChain->SetBranchAddress("bhdigi_ieta"		, &bhdigi_ieta, &b_bhdigi_ieta);
 if (fChain->GetBranch( "bhdigi_iphi"			)!=0 )     fChain->SetBranchAddress("bhdigi_iphi"		, &bhdigi_iphi, &b_bhdigi_iphi);
 if (fChain->GetBranch( "bhdigi_eta"			)!=0 )     fChain->SetBranchAddress("bhdigi_eta" 		, &bhdigi_eta, &b_bhdigi_eta);
 if (fChain->GetBranch( "bhdigi_phi"			)!=0 )     fChain->SetBranchAddress("bhdigi_phi" 		, &bhdigi_phi, &b_bhdigi_phi);
 if (fChain->GetBranch( "bhdigi_z"			)!=0 )     fChain->SetBranchAddress("bhdigi_z"			, &bhdigi_z, &b_bhdigi_z);
 if (fChain->GetBranch( "bhdigi_data"			)!=0 )     fChain->SetBranchAddress("bhdigi_data"		, &bhdigi_data, &b_bhdigi_data);
 if (fChain->GetBranch( "tc_n"				)!=0 )     fChain->SetBranchAddress("tc_n"			, &tc_n, &b_tc_n);
 if (fChain->GetBranch( "tc_id"				)!=0 )     fChain->SetBranchAddress("tc_id"			, &tc_id, &b_tc_id);
 if (fChain->GetBranch( "tc_subdet"			)!=0 )     fChain->SetBranchAddress("tc_subdet"			, &tc_subdet, &b_tc_subdet);
 if (fChain->GetBranch( "tc_zside"			)!=0 )     fChain->SetBranchAddress("tc_zside"			, &tc_zside, &b_tc_zside);
 if (fChain->GetBranch( "tc_layer"			)!=0 )     fChain->SetBranchAddress("tc_layer"			, &tc_layer, &b_tc_layer);
 if (fChain->GetBranch( "tc_wafer"			)!=0 )     fChain->SetBranchAddress("tc_wafer"			, &tc_wafer, &b_tc_wafer);
 if (fChain->GetBranch( "tc_wafertype"			)!=0 )     fChain->SetBranchAddress("tc_wafertype"		, &tc_wafertype, &b_tc_wafertype);
 if (fChain->GetBranch( "tc_cell"			)!=0 )     fChain->SetBranchAddress("tc_cell"			, &tc_cell, &b_tc_cell);
 if (fChain->GetBranch( "tc_data"			)!=0 )     fChain->SetBranchAddress("tc_data"			, &tc_data, &b_tc_data);
 if (fChain->GetBranch( "tc_pt"				)!=0 )     fChain->SetBranchAddress("tc_pt"			, &tc_pt, &b_tc_pt);
 if (fChain->GetBranch( "tc_mipPt"			)!=0 )     fChain->SetBranchAddress("tc_mipPt"			, &tc_mipPt, &b_tc_mipPt);
 if (fChain->GetBranch( "tc_energy"			)!=0 )     fChain->SetBranchAddress("tc_energy"			, &tc_energy, &b_tc_energy);
 if (fChain->GetBranch( "tc_eta"			)!=0 )     fChain->SetBranchAddress("tc_eta"			, &tc_eta, &b_tc_eta);
 if (fChain->GetBranch( "tc_phi"			)!=0 )     fChain->SetBranchAddress("tc_phi"			, &tc_phi, &b_tc_phi);
 if (fChain->GetBranch( "tc_x"				)!=0 )     fChain->SetBranchAddress("tc_x"			, &tc_x, &b_tc_x);
 if (fChain->GetBranch( "tc_y"				)!=0 )     fChain->SetBranchAddress("tc_y"			, &tc_y, &b_tc_y);
 if (fChain->GetBranch( "tc_z"				)!=0 )     fChain->SetBranchAddress("tc_z"			, &tc_z, &b_tc_z);
 if (fChain->GetBranch( "tc_cluster_id"			)!=0 )     fChain->SetBranchAddress("tc_cluster_id"		, &tc_cluster_id, &b_tc_cluster_id);
 if (fChain->GetBranch( "tc_multicluster_id"		)!=0 )     fChain->SetBranchAddress("tc_multicluster_id" 	, &tc_multicluster_id, &b_tc_multicluster_id);
 if (fChain->GetBranch( "tc_multicluster_pt"		)!=0 )     fChain->SetBranchAddress("tc_multicluster_pt" 	, &tc_multicluster_pt, &b_tc_multicluster_pt);
 if (fChain->GetBranch( "tc_super_n"			)!=0 )     fChain->SetBranchAddress("tc_super_n" 		, &tc_super_n, &b_tc_super_n);
 if (fChain->GetBranch( "tc_ave_dist"			)!=0 )     fChain->SetBranchAddress("tc_ave_dist"		, &tc_ave_dist, &b_tc_ave_dist);
 if (fChain->GetBranch( "cl_n"				)!=0 )     fChain->SetBranchAddress("cl_n"			, &cl_n, &b_cl_n);
 if (fChain->GetBranch( "cl_id"				)!=0 )     fChain->SetBranchAddress("cl_id"			, &cl_id, &b_cl_id);
 if (fChain->GetBranch( "cl_mipPt"			)!=0 )     fChain->SetBranchAddress("cl_mipPt"			, &cl_mipPt, &b_cl_mipPt);
 if (fChain->GetBranch( "cl_pt"				)!=0 )     fChain->SetBranchAddress("cl_pt"			, &cl_pt, &b_cl_pt);
 if (fChain->GetBranch( "cl_energy"			)!=0 )     fChain->SetBranchAddress("cl_energy"			, &cl_energy, &b_cl_energy);
 if (fChain->GetBranch( "cl_eta"			)!=0 )     fChain->SetBranchAddress("cl_eta"			, &cl_eta, &b_cl_eta);
 if (fChain->GetBranch( "cl_phi"			)!=0 )     fChain->SetBranchAddress("cl_phi"			, &cl_phi, &b_cl_phi);
 if (fChain->GetBranch( "cl_layer"			)!=0 )     fChain->SetBranchAddress("cl_layer"			, &cl_layer, &b_cl_layer);
 if (fChain->GetBranch( "cl_subdet"			)!=0 )     fChain->SetBranchAddress("cl_subdet"			, &cl_subdet, &b_cl_subdet);
 if (fChain->GetBranch( "cl_cells_n"			)!=0 )     fChain->SetBranchAddress("cl_cells_n" 		, &cl_cells_n, &b_cl_cells_n);
 if (fChain->GetBranch( "cl_cells_id"			)!=0 )     fChain->SetBranchAddress("cl_cells_id"		, &cl_cells_id, &b_cl_cells_id);
 if (fChain->GetBranch( "cl_multicluster_id"		)!=0 )     fChain->SetBranchAddress("cl_multicluster_id" 	, &cl_multicluster_id, &b_cl_multicluster_id);
 if (fChain->GetBranch( "cl_multicluster_pt"		)!=0 )     fChain->SetBranchAddress("cl_multicluster_pt" 	, &cl_multicluster_pt, &b_cl_multicluster_pt);
 if (fChain->GetBranch( "cl3d_n"			)!=0 )     fChain->SetBranchAddress("cl3d_n"			, &cl3d_n, &b_cl3d_n);
 if (fChain->GetBranch( "cl3d_id"			)!=0 )     fChain->SetBranchAddress("cl3d_id"			, &cl3d_id, &b_cl3d_id);
 if (fChain->GetBranch( "cl3d_pt"			)!=0 )     fChain->SetBranchAddress("cl3d_pt"			, &cl3d_pt, &b_cl3d_pt);
 if (fChain->GetBranch( "cl3d_energy"			)!=0 )     fChain->SetBranchAddress("cl3d_energy"		, &cl3d_energy, &b_cl3d_energy);
 if (fChain->GetBranch( "cl3d_eta"			)!=0 )     fChain->SetBranchAddress("cl3d_eta"			, &cl3d_eta, &b_cl3d_eta);
 if (fChain->GetBranch( "cl3d_phi"			)!=0 )     fChain->SetBranchAddress("cl3d_phi"			, &cl3d_phi, &b_cl3d_phi);
 if (fChain->GetBranch( "cl3d_clusters_n"		)!=0 )     fChain->SetBranchAddress("cl3d_clusters_n"		, &cl3d_clusters_n, &b_cl3d_clusters_n);
 if (fChain->GetBranch( "cl3d_clusters_id"		)!=0 )     fChain->SetBranchAddress("cl3d_clusters_id"		, &cl3d_clusters_id, &b_cl3d_clusters_id);
 if (fChain->GetBranch( "cl3d_showerlength"		)!=0 )     fChain->SetBranchAddress("cl3d_showerlength"		, &cl3d_showerlength, &b_cl3d_showerlength);
 if (fChain->GetBranch( "cl3d_coreshowerlength"		)!=0 )     fChain->SetBranchAddress("cl3d_coreshowerlength"	, &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
 if (fChain->GetBranch( "cl3d_firstlayer"		)!=0 )     fChain->SetBranchAddress("cl3d_firstlayer"		, &cl3d_firstlayer, &b_cl3d_firstlayer);
 if (fChain->GetBranch( "cl3d_maxlayer"			)!=0 )     fChain->SetBranchAddress("cl3d_maxlayer"		, &cl3d_maxlayer, &b_cl3d_maxlayer);
 if (fChain->GetBranch( "cl3d_seetot"			)!=0 )     fChain->SetBranchAddress("cl3d_seetot"		, &cl3d_seetot, &b_cl3d_seetot);
 if (fChain->GetBranch( "cl3d_seemax"			)!=0 )     fChain->SetBranchAddress("cl3d_seemax"		, &cl3d_seemax, &b_cl3d_seemax);
 if (fChain->GetBranch( "cl3d_spptot"			)!=0 )     fChain->SetBranchAddress("cl3d_spptot"		, &cl3d_spptot, &b_cl3d_spptot);
 if (fChain->GetBranch( "cl3d_sppmax"			)!=0 )     fChain->SetBranchAddress("cl3d_sppmax"		, &cl3d_sppmax, &b_cl3d_sppmax);
 if (fChain->GetBranch( "cl3d_szz"			)!=0 )     fChain->SetBranchAddress("cl3d_szz"			, &cl3d_szz, &b_cl3d_szz);
 if (fChain->GetBranch( "cl3d_srrtot"			)!=0 )     fChain->SetBranchAddress("cl3d_srrtot"		, &cl3d_srrtot, &b_cl3d_srrtot);
 if (fChain->GetBranch( "cl3d_srrmax"			)!=0 )     fChain->SetBranchAddress("cl3d_srrmax"		, &cl3d_srrmax, &b_cl3d_srrmax);
 if (fChain->GetBranch( "cl3d_srrmean"			)!=0 )     fChain->SetBranchAddress("cl3d_srrmean"		, &cl3d_srrmean, &b_cl3d_srrmean);
 if (fChain->GetBranch( "cl3d_emaxe"			)!=0 )     fChain->SetBranchAddress("cl3d_emaxe" 		, &cl3d_emaxe, &b_cl3d_emaxe);
 if (fChain->GetBranch( "cl3d_bdteg"			)!=0 )     fChain->SetBranchAddress("cl3d_bdteg" 		, &cl3d_bdteg, &b_cl3d_bdteg);
 if (fChain->GetBranch( "cl3d_quality"			)!=0 )     fChain->SetBranchAddress("cl3d_quality"		, &cl3d_quality, &b_cl3d_quality);
 if (fChain->GetBranch( "tower_n"			)!=0 )     fChain->SetBranchAddress("tower_n"			, &tower_n, &b_tower_n);
 if (fChain->GetBranch( "tower_pt"			)!=0 )     fChain->SetBranchAddress("tower_pt"			, &tower_pt, &b_tower_pt);
 if (fChain->GetBranch( "tower_energy"			)!=0 )     fChain->SetBranchAddress("tower_energy"		, &tower_energy, &b_tower_energy);
 if (fChain->GetBranch( "tower_eta"			)!=0 )     fChain->SetBranchAddress("tower_eta"			, &tower_eta, &b_tower_eta);
 if (fChain->GetBranch( "tower_phi"			)!=0 )     fChain->SetBranchAddress("tower_phi"			, &tower_phi, &b_tower_phi);
 if (fChain->GetBranch( "tower_etEm"			)!=0 )     fChain->SetBranchAddress("tower_etEm" 		, &tower_etEm, &b_tower_etEm);
 if (fChain->GetBranch( "tower_etHad"			)!=0 )     fChain->SetBranchAddress("tower_etHad"		, &tower_etHad, &b_tower_etHad);
 if (fChain->GetBranch( "tower_iEta"			)!=0 )     fChain->SetBranchAddress("tower_iEta" 		, &tower_iEta, &b_tower_iEta);
 if (fChain->GetBranch( "tower_iPhi"			)!=0 )     fChain->SetBranchAddress("tower_iPhi" 		, &tower_iPhi, &b_tower_iPhi);
   Notify();
}

Bool_t BuildTreeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BuildTreeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BuildTreeBase::Cut()
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BuildTreeBase_cxx
