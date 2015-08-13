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
    void draw(std::vector<TH1*>& hists,Option_t* drawOpt="",bool legend=false,bool logx=false,bool logy=false,bool logz=false);
    void draw(TH1* hist,Option_t* drawOpt="",bool legend=false,bool logx=false,bool logy=false,bool logz=false);
    void newLegend(TString header="",float lx=0.8,float ly=0.8,float ux=0.9,float uy=0.88);
    void close();

  private:
    TCanvas*   mCanvas;
    TPad*      mPad;
    TPDF*      mPDF;
    TPaveText* mPageHeader;
    TLegend*   mLegend;
    int        mNCanvasX;
    int        mNCanvasY;
    int        mSubCanvasIndex;
    int        mMaxCanvasPerPage;
};
#endif
