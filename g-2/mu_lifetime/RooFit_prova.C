#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFit.h"
#include <vector>


using namespace std;
using namespace RooFit;

void RooFit_prova(){

    double time;
    double cal_time;

    TFile *file = new TFile("myFile.root", "RECREATE");

    TTree *tree = new TTree("myTree", "");

    tree->Branch("cal_time_var", &cal_time);

    string line;

    double pendenza = 6.66658278e+01;
    double intercetta= 9.73509944e-03;

    //cout << 0.1*pendenza+intercetta << endl;

    ifstream infile("output.txt");
    while(getline(infile,line))
    {
        stringstream ss(line);
        ss >> time;

        if(time>7.){

            cal_time=(time-intercetta)/pendenza;
            tree->Fill(); 
        }

    }

    file->Write();

    infile.close();

    file->Close();

   /* struct RooFit_prova
    {
        double time;
    };

    RooFit_prova tempo;
    */

    TFile *File = TFile::Open("myFile.root", "READ");

    TTree *data_tree = (TTree*) File->Get("myTree");

    TBranch *branch = data_tree->GetBranch("cal_time_var");

    //branch->SetAddress(&tempo.time);

    RooRealVar cal_time_var("cal_time_var", "cal_time_var", 0.5, 15);

    //RooDataSet data("data", "data", data_tree, RooArgSet(cal_time_var));

    RooDataSet data("data", "data", RooArgSet(cal_time_var), Import(*data_tree));

    //fit parameters
    RooRealVar tau("tau", "mean lifetime",2.2, 2.0, 3.); 

    RooRealVar bl("bl", "baseline", 1, 0, 5);

    RooRealVar Nf("Nf", "Nf", 4200, 3500, 5500);

    //fit pdfs
    RooPolynomial pol0("pol0","pol0", cal_time_var ,bl);

    RooGenericPdf exp("exponential", "exp(-cal_time_var*(1/tau))", RooArgSet(cal_time_var,tau));

    //fraction of the first pdf (second one's fraction is 1-fraction)
    RooRealVar  fraction("fraction", "fraction", 0.9,0,1);

    RooAddPdf model("model", "model", RooArgList(exp, pol0), RooArgList(fraction));

    model.fitTo(data);

    RooPlot *cal_time_frame = cal_time_var.frame();

    data.plotOn(cal_time_frame);
    model.plotOn(cal_time_frame);

    cal_time_frame->Draw();

    /*int j=0;
    int N = data_tree->GetEntries();
    for(int i =0; i<N; i++)
    {
        branch->GetEntry(i);
        if(tempo.time<1.08 && tempo.time>0.5){j=j+1;}
    }
    cout << endl << endl << j << "  " << N <<  endl;

    cal_time_var.setRange("range", 0.5, 1.08);
    RooAbsReal *intModel= model.createIntegral(cal_time_var, NormSet(cal_time_var), Range("range"));
    
    double integrale = intModel->getVal();

    RooAbsReal *totalIntegral= model.createIntegral(cal_time_var);
    double intTOT= totalIntegral->getVal(); //to check that it gives 1 

    cout << endl <<  integrale << "  " << integrale*N  <<  endl;
    */
}