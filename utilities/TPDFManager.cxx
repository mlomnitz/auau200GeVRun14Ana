#include <iostream>
#include <algorithm>
#include <memory>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TPad.h"
#include "TPDF.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1.h"
#include "TLine.h"
#include "TGraph.h"

#include "TPDFManager.h"

TPDFManager::TPDFManager(TString outFileName)
{
   mCanvas = new TCanvas("mCanvas", "", 800, 450);
   mCanvas->cd();
   mPad = new TPad("mPad", "", 0, 0, 1, 0.95);
   mPad->Draw();
   mPDF = new TPDF(Form("%s.pdf", outFileName.Data()));
}

void TPDFManager::drawFirstPage(TString input)
{
   mPad->cd();
   mPad->Draw();
   mCanvas->cd();
   mCanvas->Update();
   mSubCanvasIndex = -1;

   /*mPageHeader = new TPaveText(0,0.95,1,1);
   mPageHeader->SetFillColor(0);
   mPageHeader->SetTextColor(1);
   mPageHeader->SetTextFont(42);
   mPageHeader->AddText(headerText);
   mPageHeader->AddLine(0.0,0.95,1.0,0.95);
   mPageHeader->Draw();
   */

   TPaveText* body = new TPaveText(0.05, 0.05, 0.95, 0.9);
   body->AddText("");
   body->SetFillColor(0);
   body->SetTextColor(1);
   body->SetTextFont(42);
   body->SetTextSize(0.055);
   body->SetBorderSize(0);

   std::ifstream in(input.Data());
   std::string line;
   getline(in, line);

   while (line != "----")
   {
      body->AddText(line.c_str());
      getline(in, line);
   }

   // mPad->cd();
   body->Draw();
}

void TPDFManager::newPage(int const nX, int const nY, TString const headerText)
{
   mNCanvasX = nX;
   mNCanvasY = nY;

   if (mSubCanvasIndex != 0)
   {
      mPad->Draw();
      mCanvas->cd();
      mCanvas->Update();
      mPDF->NewPage();
      delete mPad;
      mPad = new TPad("mPad", "", 0, 0, 1, 0.95);
      mPad->Draw();
   }

   if (headerText != "HEADER")
   {
      delete mPageHeader;
      mPageHeader = new TPaveText(0, 0.95, 1, 1);
      mPageHeader->SetFillColor(0);
      mPageHeader->SetTextColor(1);
      mPageHeader->SetTextFont(42);
      mPageHeader->AddText(headerText);
      mPageHeader->AddLine(0.0, 0.95, 1.0, 0.95);
   }

   mPageHeader->Draw();
   mPad->cd();
   mPad->Divide(mNCanvasX, mNCanvasY);
   mPad->SetGrid(0, 0);
   mMaxCanvasPerPage = mNCanvasX * mNCanvasY;
   mSubCanvasIndex = 0;
}

void TPDFManager::draw(std::vector<TH1*>& objects, Option_t* drawOpt, TLegend* legend, bool logx, bool logy, bool logz)
{
   if (mSubCanvasIndex + 1 > mMaxCanvasPerPage)
      newPage(mNCanvasX, mNCanvasY);

   mSubCanvasIndex += 1;
   mPad->cd(mSubCanvasIndex);
   mPad->cd(mSubCanvasIndex)->Draw();

   objects.at(0)->Draw(drawOpt);
   if (legend)
   {
      TString header = legend->GetHeader();
      legend->Clear();
      legend->SetHeader(header.Data());
      legend->AddEntry(objects.at(0), objects.at(0)->GetTitle(), "LP");
   }

   for (size_t iHist = 1; iHist < objects.size(); ++iHist)
   {
      objects.at(iHist)->Draw(Form("sames%s", drawOpt));
      if (legend) legend->AddEntry(objects.at(iHist), objects.at(iHist)->GetTitle(), "LP");
   }

   if (legend) legend->Draw();
   if (logx)   mPad->cd(mSubCanvasIndex)->SetLogx();
   if (logy)   mPad->cd(mSubCanvasIndex)->SetLogy();
   if (logz)   mPad->cd(mSubCanvasIndex)->SetLogz();

   mCanvas->cd();
   mCanvas->Update();
   mPad->Draw();
}

/*
void TPDFManager::draw(std::vector<TH1*>& hists,Option_t* drawOpt,TLegend* legend,bool logx,bool logy,bool logz)
{
  std::vector<TObject*> objects(hists.size());

  for(auto it = hists.begin(); it!=hists.end();++it)
  {
    objects.push_back(*it);
  }

  draw(objects,drawOpt,legend,logx,logy,logz);
}
*/

void TPDFManager::drawOnSamePad(std::vector<TObject*>& obj, Option_t* drawOpt)
{
   mPad->cd(mSubCanvasIndex);

   obj.at(0)->Draw(drawOpt);

   for (size_t iHist = 1; iHist < obj.size(); ++iHist)
   {
      obj.at(iHist)->Draw(Form("sames%s", drawOpt));
   }

   mCanvas->cd();
   mCanvas->Update();
   mPad->Draw();
}

void TPDFManager::drawOnSamePad(TObject* obj, Option_t* drawOpt)
{
   mPad->cd(mSubCanvasIndex);

   obj->Draw(drawOpt);

   mCanvas->cd();
   mCanvas->Update();
   mPad->Draw();
}

void TPDFManager::drawLine(float x1, float y1, float x2, float y2)
{
   mLines.push_back(new TLine(x1, y1, x2, y2));
   mLines.back()->SetLineStyle(2);
   mLines.back()->SetLineWidth(0.8);
   mPad->cd(mSubCanvasIndex);
   mLines.back()->Draw("same");

   mCanvas->cd();
   mCanvas->Update();
   mPad->Draw();
}

void TPDFManager::drawShadedArea(float x1, float y1, float x2, float y2)
{
   mShadedAreas.push_back(new TGraph(4));
   mShadedAreas.back()->SetPoint(0, x1, y1);
   mShadedAreas.back()->SetPoint(1, x2, y1);
   mShadedAreas.back()->SetPoint(2, x2, y2);
   mShadedAreas.back()->SetPoint(3, x1, y2);
   mShadedAreas.back()->SetLineWidth(0.8);
   mShadedAreas.back()->SetLineColor(1);
   mShadedAreas.back()->SetFillColor(1);
   mShadedAreas.back()->SetFillStyle(3853);

   mPad->cd(mSubCanvasIndex);
   mShadedAreas.back()->Draw("same:F");

   mCanvas->cd();
   mCanvas->Update();
   mPad->Draw();
}

void TPDFManager::draw(TH1* hist, Option_t* drawOpt, TLegend* legend, bool logx, bool logy, bool logz)
{
   std::vector<TH1*> temp;
   temp.push_back(hist);

   draw(temp, drawOpt, legend, logx, logy, logz);
}

void TPDFManager::draw(TH1* hist0, TH1* hist1, Option_t* drawOpt, TLegend* legend, bool logx, bool logy, bool logz)
{
   std::vector<TH1*> temp;
   temp.push_back(hist0);
   temp.push_back(hist1);

   draw(temp, drawOpt, legend, logx, logy, logz);
}

TLegend* TPDFManager::newLegend(TString header, float lx, float ly, float ux, float uy)
{
   TLegend* leg = new TLegend(lx, ly, ux, uy);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.055);
   if (header.Length()) leg->SetHeader(header.Data());

   mLegends.push_back(leg);

   return leg;
}

void TPDFManager::close()
{
   mPad->Draw();
   mCanvas->cd();
   mCanvas->Update();
   mPDF->Close();
   std::for_each(mLegends.begin(), mLegends.end(), std::default_delete<TLegend>());
   std::for_each(mLines.begin(), mLines.end(), std::default_delete<TLine>());
   std::for_each(mShadedAreas.begin(), mShadedAreas.end(), std::default_delete<TGraph>());
}
