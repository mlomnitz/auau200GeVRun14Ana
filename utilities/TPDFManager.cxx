#include "TCanvas.h"
#include "TPad.h"
#include "TPDF.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1.h"

#include "TPDFManager.h"

TPDFManager::TPDFManager(TString outFileName): mCanvas(NULL),
  mPad(NULL), mPDF(NULL),
  mPageHeader(NULL), mLegend(NULL), mNCanvasX(1), mNCanvasY(1), mSubCanvasIndex(0), mMaxCanvasPerPage(1)
{
    mCanvas = new TCanvas("mCanvas","",800,450);
    mCanvas->cd();
    mPad = new TPad("mPad","",0,0,1,0.95);
    mPad->Draw();
    mPDF = new TPDF(Form("%s.pdf",outFileName.Data()));
}

void TPDFManager::newPage(int const nX,int const nY,TString const headerText)
{
  mNCanvasX = nX;
  mNCanvasY = nY;

  if(mSubCanvasIndex != 0)
  {
    mPad->Draw();
    mCanvas->cd();
    mCanvas->Update();
    mPDF->NewPage();
    delete mPad;
    mPad = new TPad("mPad","",0,0,1,0.95);
    mPad->Draw();
  }

  if(headerText!="HEADER")
  {
    delete mPageHeader;
    mPageHeader = new TPaveText(0,0.95,1,1);
    mPageHeader->SetFillColor(0);
    mPageHeader->SetTextColor(1);
    mPageHeader->SetTextFont(22);
    mPageHeader->AddText(headerText);
    mPageHeader->AddLine(0.0,0.95,1.0,0.95);
  }

  mPageHeader->Draw();
  mPad->cd();
  mPad->Divide(mNCanvasX,mNCanvasY);
  mPad->SetGrid(0,0);
  mMaxCanvasPerPage = mNCanvasX*mNCanvasY;
  mSubCanvasIndex = 0;
}

void TPDFManager::draw(std::vector<TH1*>& hists,Option_t* drawOpt,bool legend,bool logx,bool logy,bool logz)
{
  if(mSubCanvasIndex + 1 > mMaxCanvasPerPage)
    newPage(mNCanvasX,mNCanvasY);

  mSubCanvasIndex +=1;
  mPad->cd(mSubCanvasIndex);
  mPad->cd(mSubCanvasIndex)->Draw();

  hists.at(0)->Draw(drawOpt);
  if(legend) mLegend->AddEntry(hists.at(0),hists.at(0)->GetTitle(),"LP");

  for(size_t iHist=1;iHist<hists.size();++iHist)
  {
    hists.at(iHist)->Draw(Form("sames:%s",drawOpt));
    if(legend) mLegend->AddEntry(hists.at(iHist),hists.at(iHist)->GetTitle(),"LP");
  }

  if(legend) mLegend->Draw();
  if(logx)   mPad->cd(mSubCanvasIndex)->SetLogx();
  if(logy)   mPad->cd(mSubCanvasIndex)->SetLogy();
  if(logz)   mPad->cd(mSubCanvasIndex)->SetLogz();

  mCanvas->cd();
  mCanvas->Update();
  mPad->Draw();
}

void TPDFManager::draw(TH1* hist,Option_t* drawOpt,bool legend,bool logx,bool logy,bool logz)
{
    std::vector<TH1*> temp;
    temp.push_back(hist);

    draw(temp,drawOpt,legend,logx,logy,logz);
}

void TPDFManager::newLegend(TString header,float lx,float ly,float ux,float uy)
{
  mLegend = new TLegend(lx,ly,ux,uy);
  mLegend->SetFillColor(0);
  mLegend->SetBorderSize(0);
  mLegend->SetTextFont(132);
  mLegend->SetTextSize(0.035);
  if(header.Length()) mLegend->SetHeader(header.Data());
}

void TPDFManager::close()
{
  mPad->Draw();
  mCanvas->cd();
  mCanvas->Update();
  mPDF->Close();
}
