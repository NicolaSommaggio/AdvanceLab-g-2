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

    ifstream infile("/home/nicola/Scrivania/LAB/AdvanceLab-g-2/g-2/Data/Fondo/output.txt");
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

    struct RooFit_prova
    {
        double time;
    };

    RooFit_prova tempo;
    

    TFile *File = TFile::Open("myFile.root", "READ");

    TTree *data_tree = (TTree*) File->Get("myTree");

    TBranch *branch = data_tree->GetBranch("cal_time_var");

    branch->SetAddress(&tempo.time);

    RooRealVar cal_time_var("cal_time_var", "cal_time_var", 0.5, 15);

    //RooDataSet data("data", "data", data_tree, RooArgSet(cal_time_var));

    RooDataSet data("data", "data", RooArgSet(cal_time_var), Import(*data_tree));

    //fit parameters
    RooRealVar tau("tau", "mean lifetime",2.2, 1.5, 3.5); 

    RooRealVar bl("bl", "baseline", 138., 100., 200.);

    RooRealVar lambda("lambda", "lambda", 1.7e-6, 1e-6, 3e-6); 

    //fit pdfs
    RooPolynomial pol0("pol0","pol0",bl);

    RooGenericPdf dark("background", "lambda*exp(-lambda*cal_time_var)", RooArgSet(cal_time_var,lambda));

    RooGenericPdf exp("exponential", "(1/tau)*exp(-cal_time_var*(1/tau))", RooArgSet(cal_time_var,tau));

    //fraction of the first pdf (second one's fraction is 1-fraction)
    RooRealVar  fraction("fraction", "fraction", 0.7,0,1);

    RooAddPdf model("model", "model", RooArgList(exp, dark), RooArgList(fraction));

    model.fitTo(data);

    RooPlot *cal_time_frame = cal_time_var.frame();

    data.plotOn(cal_time_frame, Name("data"));
    model.plotOn(cal_time_frame, Name("Model"));
    exp.plotOn(cal_time_frame, LineColor(kRed));
    dark.plotOn(cal_time_frame, LineColor(kGreen), Name("exp background"));
    //pol0.plotOn(cal_time_frame, LineColor(kMagenta), Name("Linear background"));

    TLegend *legend = new TLegend();
    legend->AddEntry(cal_time_frame->getObject(0), "Data", "p");
    legend->AddEntry(cal_time_frame->getObject(1), "Total fit", "l");
    legend->AddEntry(cal_time_frame->getObject(3), "Background", "l");
    

   //frame_pull->GetYaxis()->SetTitle("Pulls");

   cout << "chi^2 = " << cal_time_frame->chiSquare() << endl;

      // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
   // -------------------------------------------------------
 
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hresid = cal_time_frame->residHist("data", "model");
 
   // Construct a histogram with the pulls of the data w.r.t the curve
   RooHist *hpull = cal_time_frame->pullHist("data", "Model");
 
   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame2 = cal_time_var.frame(Title("Residual Distribution"), Bins(10));
   frame2->addPlotable(hresid, "P");
 
   // Create a new frame to draw the pull distribution and add the distribution to the frame
   RooPlot *frame3 = cal_time_var.frame(Title("Pull Distribution"), Bins(10));
   frame3->addPlotable(hpull, "P");

    cal_time_frame->Draw();
    legend->Draw();
    new TCanvas();
    frame2->Draw();
    new TCanvas();
    frame3->Draw();

}