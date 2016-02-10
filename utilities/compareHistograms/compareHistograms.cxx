/* **************************************************
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "TApplication.h"
#include "TH1D.h"
#include "TH3F.h"
#include "THn.h"
#include "TF1.h"
#include "TFile.h"
#include "TDirectory.h"

#include "StyleUtilities.h"

using namespace std;

void draw(TH1* h1, TH1* h2)
{
  TCanvas* cv = canvas();
  cv->cd();
  setStyle(h1,20,1);
  setStyle(h1,20,2);
  h1->Draw();
  h2->Draw("sames");
}

vector<double> getChi2OverNdf(TObject* obj1, TObject* obj2)
{
  if(!obj1 || !obj2)
  {
    cout << "An object doesn't exit!" << endl;
    return {std::numeric_limits<double>::quiet_NaN()};
  }

  if(TH1* h1 = dynamic_cast<TH1*>(obj1))
  {
    TH1* h2 = static_cast<TH1*>(obj2);
    double chi2OverNdf = h1->Chi2Test(h2,"Chi2/NDF");

    if(chi2OverNdf) draw(h1,h2);
    return {chi2OverNdf};
  }
  else if(THn* hn1 = dynamic_cast<THn*>(obj1))
  {
    THn* hn2 = static_cast<THn*>(obj2);

    vector<double> chi2OverNdf;
    // for(int id=0; id < hn1->GetNdimensions(); ++id)
    for(int id=0; id < 1; ++id)
    {
      TH1* h1 = (TH1*)hn1->Projection(id);
      TH1* h2 = (TH1*)hn2->Projection(id);
      chi2OverNdf.push_back(h1->Chi2Test(h2,"Chi2/NDF"));

      if(chi2OverNdf.back()) draw(h1,h2);
      // if(chi2OverNdf.back()) cout << h1->GetName() << endl;

      // delete h1;
      // delete h2;
    }

    return chi2OverNdf;
  }

  return {};
}

TH1F* getChi2OverNdf(string const f1Name, string const f2Name)
{
  TFile* f1 = new TFile(f1Name.c_str());
  TFile* f2 = new TFile(f2Name.c_str());

  TList* listOfKeys=f1->GetListOfKeys();

  if(!listOfKeys->GetEntries())
  {
    cout<< f1->GetName() << " is empty. Exiting ..."<<endl;
    delete listOfKeys;
    return nullptr;
  }

  //Loop over list of keys
  vector<double> chi2OverNdf;
  TH1F* hChi2OverNdf = new TH1F("hChi2OverNdf","hChi2Ndf",100,0,10);

  for(int iKey=0; iKey<listOfKeys->GetEntries(); iKey++)
  {
    for(auto x : getChi2OverNdf(f1->Get(listOfKeys->At(iKey)->GetName()),
                                f2->Get(listOfKeys->At(iKey)->GetName())))
    {
      chi2OverNdf.push_back(x);
      hChi2OverNdf->Fill(x);
    }
  }

  TGraph* grChi2OverNdf = new TGraph(chi2OverNdf.size());
  for(size_t ie = 0; ie < chi2OverNdf.size(); ++ie)
  {
    grChi2OverNdf->SetPoint(ie,ie,chi2OverNdf[ie]);
  }

  return hChi2OverNdf;
}

int main(int argc, char* argv[])
{
  if(argc<3)
  {
    std::cout<< "Usage: "<< argv[0] << " file0.root file1.root" << std::endl;
    return 0;
  }

  TH1::AddDirectory(false);
  string fileName0{argv[1]};
  string fileName1{argv[2]};


  TApplication theApp("tapp", &argc,argv);
  TH1F* hChi2OverNdf = getChi2OverNdf(fileName0,fileName1);

  TCanvas* cv = canvas();
  cv->cd();
  setStyle(hChi2OverNdf,20,1);
  hChi2OverNdf->Draw("P");

  theApp.Run();
}
