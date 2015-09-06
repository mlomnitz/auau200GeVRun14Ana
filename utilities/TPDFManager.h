 /* **************************************************
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#ifndef TPDFManager_H
#define TPDFManager_H

#include <vector>
#include "TString.h"

class TCanvas;
class TPad;
class TPDF;
class TPaveText;
class TLegend;
class TH1;
class TObject;
class TLine;
class TGraph;

class TPDFManager
{
  public:
    TPDFManager(TString outFileName);
    void newPage(int nX=1,int nY=1,TString headerText="HEADER");
    void drawFirstPage(TString inputFile);
    void draw(TH1* hist,Option_t* drawOpt="",TLegend* legend=NULL,bool logx=false,bool logy=false,bool logz=false);
    void draw(TH1* hist0,TH1* hist1,Option_t* drawOpt="",TLegend* legend=NULL,bool logx=false,bool logy=false,bool logz=false);
    void draw(std::vector<TH1*>& hists,Option_t* drawOpt="",TLegend* legend=NULL,bool logx=false,bool logy=false,bool logz=false);
    // void draw(std::vector<TObject*>& objects,Option_t* drawOpt="",TLegend* legend=NULL,bool logx=false,bool logy=false,bool logz=false);
    void drawOnSamePad(std::vector<TObject*>& obj,Option_t* drawOpt="");
    void drawOnSamePad(TObject* obj,Option_t* drawOpt="same");
    void drawLine(float x1,float y1,float x2,float y2);
    void drawShadedArea(float x1,float y1,float x2,float y2);
    TLegend* newLegend(TString header="",float lx=0.15,float ly=0.6,float ux=0.4,float uy=0.85);
    void close();

  private:
    TCanvas*   mCanvas = nullptr;
    TPad*      mPad = nullptr;
    TPDF*      mPDF = nullptr;
    TPaveText* mPageHeader = nullptr;
    // TLegend*   mLegend = nullptr;
    std::vector<TLegend*> mLegends;
    std::vector<TLine*> mLines;
    std::vector<TGraph*> mShadedAreas;
    int        mNCanvasX = 1;
    int        mNCanvasY = 1;
    int        mSubCanvasIndex = 0;
    int        mMaxCanvasPerPage = 0;
};
#endif
