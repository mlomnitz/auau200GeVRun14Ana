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

class TPDFManager
{
  public:
    TPDFManager(TString outFileName);
    void newPage(int nX=1,int nY=1,TString headerText="HEADER");
    void draw(TH1* hist,Option_t* drawOpt="",bool legend=false,bool logx=false,bool logy=false,bool logz=false);
    void draw(TH1* hist0,TH1* hist1,Option_t* drawOpt="",bool legend=false,bool logx=false,bool logy=false,bool logz=false);
    void draw(std::vector<TH1*>& hists,Option_t* drawOpt="",bool legend=false,bool logx=false,bool logy=false,bool logz=false);
    void newLegend(TString header="",float lx=0.15,float ly=0.6,float ux=0.4,float uy=0.85);
    void close();

  private:
    TCanvas*   mCanvas = NULL;
    TPad*      mPad = NULL;
    TPDF*      mPDF = NULL;
    TPaveText* mPageHeader = NULL;
    TLegend*   mLegend = NULL;
    int        mNCanvasX = 1;
    int        mNCanvasY = 1;
    int        mSubCanvasIndex = 0;
    int        mMaxCanvasPerPage = 0;
};
#endif
