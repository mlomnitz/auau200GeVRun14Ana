/*
 * File:   saveRootFileContent.C
 * Author: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * This macro is to save all histograms stored in a root file.
 * Created on October 8, 2011, 3:05 PM
 */
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "TUnixSystem.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

using namespace std;

void saveDirContent(TDirectory* dir,string path);
void saveHistogram(TH1* hist,string path);
void setGraphicsStyle();


string plotsExtension="png";


void saveRootFileContent(string fileName)
{
    setGraphicsStyle();

    /// Open file
    TFile* inputFile = new TFile(fileName.c_str());

    if (!inputFile)
    {
        cout << fileName << " does not exist. Existing ..." << endl;
        return;
    }
    
    /// Make a directory named after the file name
    size_t found = fileName.find(".root");
    string histsDirName = fileName.substr(0,found) + "_hists";

    if(gSystem->cd(histsDirName.c_str()))
    {
        cout<<"Directory "<<histsDirName<<" exists. Please remove it and try again. Exiting ..."<<endl;
        return;
    }

    cout<<"Creating directory "<<histsDirName<<endl;
    gSystem->mkdir(histsDirName.c_str());
    
    saveDirContent(dynamic_cast<TDirectory*>(inputFile),histsDirName);

    inputFile->Close();

    /// Make tarball and delete directory
    string command;

    cout<<"Making tarball "<<histsDirName<<".tar"<<endl;
    command = "tar -cf "+ histsDirName +".tar " + histsDirName;
    gSystem->Exec(command.c_str());
    cout<<"Removing "<<histsDirName<<" directory"<<endl;
    command = "rm -r "+ histsDirName;
    gSystem->Exec(command.c_str());

    cout<<"Done. All the generated plots are in "<<histsDirName<<".tar"<<endl;
    cout<<"GoodBye."<<endl;
}
/////////////////////////////////////////////////////////////
void saveDirContent(TDirectory* dir,string path)
{
    if(!dir) return;

    // Get list of keys in the directory
    TList* listOfKeys=dir->GetListOfKeys();

    if(!listOfKeys->GetEntries())
    {
        cout<< dir->GetName() << " is empty. Exiting ..."<<endl;
	delete listOfKeys;
        return;
    }

    //Loop over list of keys

    int nKeys = listOfKeys->GetEntries();
    TKey* key=0;
    TObject* obj=0;

    for(int iKey=0; iKey<nKeys; iKey++)
    {
        key = dynamic_cast<TKey*>(listOfKeys->At(iKey));
        obj = dynamic_cast<TObject*>(key->ReadObj());

        //Check if the obj is a directory, a hist, or something else
        if(obj->IsA()->InheritsFrom(TDirectory::Class()))
        {
            stringstream s;
            s << key->GetName();
            string dirPath = path + "/" + s.str();
            cout<<"Creating directory "<<dirPath<<endl;
            gSystem->mkdir(dirPath.c_str());
            saveDirContent(dynamic_cast<TDirectory*>(key->ReadObj()),dirPath);
        }
        else if(obj->IsA()->InheritsFrom(TH1::Class()))
        {
            saveHistogram(dynamic_cast<TH1*>(key->ReadObj()),path);
        }
        else cout<<"Object "<<key->GetName()<<" is neither a directory, nor a histogram"<<endl;
    }

    if(!dir->IsA()->InheritsFrom(TFile::Class())) delete listOfKeys;
}
/////////////////////////////////////////////////////////////
void saveHistogram(TH1* hist,string path)
{
    if(!hist) return;
    
    TCanvas* canvas=new TCanvas("canvas","canvas");

    //If hist is not two dimensional do the following. If it is two dimemnsional you shouldn't have SetFillColor(0)
    if(!hist->IsA()->InheritsFrom(TH2::Class()))
    {
        gStyle->SetFillColor(0);
        canvas->SetFillColor(0);
        canvas->SetGrid();
    }

    stringstream s;
    s << hist->GetName();
    string histName = path + "/" + s.str() + "." + plotsExtension;

    if(hist->IsA()->InheritsFrom(TH2::Class())) 
    {
	    canvas->SetLogz();
	    hist->Draw("colz");
    }
    else hist->Draw();

    canvas->SaveAs(histName.c_str());

    delete canvas;
    delete hist;
}
/////////////////////////////////////////////////////////////
void setGraphicsStyle()
{
    // **************************** Set graphic style ***************************************
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.05,"Y");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
    gStyle->SetNumberContours(NCont);
    // **************************************************************************************
}

