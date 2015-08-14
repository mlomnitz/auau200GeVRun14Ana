#include <iostream>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TFile.h"

#include "TPDFManager.h"
#include "StyleUtilities.h"

using namespace std;

int main()
{
  // create dummy histograms
  TH2F *h2 = new TH2F("h2","h2",100,-5,5,100,-5,5);
  TF2 *xyg = new TF2("xyg","xygaus",-5,5,-5,5);
  xyg->SetParameters(1,0,2,0,0.5);  //amplitude, meanx,sigmax,meany,sigmay
  h2->FillRandom("xyg",1e4);

  TH1F* h1 = new TH1F("h1","h1",100,-1,1);
  h1->FillRandom("gaus",1e4);

  setStyle(h1,24,2);
  setStyle(h2);

  setGraphicsStyle();

  TPDFManager* pdf = new TPDFManager("test");

  // plot events statistics
  pdf->newPage(2,1,"First page"); // divide the page into 2x1 canvases
  pdf->draw(h1);
  pdf->draw(h2,"colz",false,false,false,true);

  pdf->newPage(2,1); // divide the page into 2x1 canvases
  pdf->draw(h2->ProjectionX(),h2->ProjectionY()); // plot two histograms in the same canvas

  gStyle->SetOptStat(0);
  pdf->newLegend("testing legend"); // create a new legend
  pdf->draw(h2->ProjectionX(),h1,"",true); // plot two histograms in the same canvas and add legend

  pdf->close();

  return 0;
}
