//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  4 17:33:05 2015 by ROOT version 5.34/09
// from TTree kfEvent/event information and kfVertex
// found on file: /project/projectdirs/starprod/hft/kfVertex/Run14/AuAu/200GeV/physics2/P15ic/094/15094070/st_physics_15094070_raw_0000007.kfVertex.root
//////////////////////////////////////////////////////////

#ifndef kfEvent_h
#define kfEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class kfEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           mRunId   ;
   Int_t           mEventId ;
   Int_t           mRefMult ;
   Int_t           mGRefMult;
   Float_t         mVx      ;
   Float_t         mVy      ;
   Float_t         mVz      ;
   Float_t         mKfVx    ;
   Float_t         mKfVy    ;
   Float_t         mKfVz    ;

   // List of branches
   TBranch        *b_mRunId;   //!
   TBranch        *b_mEventId;   //!
   TBranch        *b_mRefMult;   //!
   TBranch        *b_mGRefMult;   //!
   TBranch        *b_mVx;   //!
   TBranch        *b_mVy;   //!
   TBranch        *b_mVz;   //!
   TBranch        *b_mKfVx;   //!
   TBranch        *b_mKfVy;   //!
   TBranch        *b_mKfVz;   //!

   kfEvent(TTree *tree=0);
   virtual ~kfEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef kfEvent_cxx
kfEvent::kfEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/project/projectdirs/starprod/hft/kfVertex/Run14/AuAu/200GeV/physics2/P15ic/094/15094070/st_physics_15094070_raw_0000007.kfVertex.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/project/projectdirs/starprod/hft/kfVertex/Run14/AuAu/200GeV/physics2/P15ic/094/15094070/st_physics_15094070_raw_0000007.kfVertex.root");
      }
      f->GetObject("kfEvent",tree);

   }
   Init(tree);
}

kfEvent::~kfEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kfEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kfEvent::LoadTree(Long64_t entry)
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

void kfEvent::Init(TTree *tree)
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

   fChain->SetBranchAddress("mRunId   ", &mRunId   , &b_mRunId);
   fChain->SetBranchAddress("mEventId ", &mEventId , &b_mEventId);
   fChain->SetBranchAddress("mRefMult ", &mRefMult , &b_mRefMult);
   fChain->SetBranchAddress("mGRefMult", &mGRefMult, &b_mGRefMult);
   fChain->SetBranchAddress("mVx      ", &mVx      , &b_mVx);
   fChain->SetBranchAddress("mVy      ", &mVy      , &b_mVy);
   fChain->SetBranchAddress("mVz      ", &mVz      , &b_mVz);
   fChain->SetBranchAddress("mKfVx    ", &mKfVx    , &b_mKfVx);
   fChain->SetBranchAddress("mKfVy    ", &mKfVy    , &b_mKfVy);
   fChain->SetBranchAddress("mKfVz    ", &mKfVz    , &b_mKfVz);
   Notify();
}

Bool_t kfEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kfEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kfEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef kfEvent_cxx
