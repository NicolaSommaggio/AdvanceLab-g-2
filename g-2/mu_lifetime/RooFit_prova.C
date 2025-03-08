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

void RooFit_prova(string inputfile){

    double time;
    double cal_time;

    TFile *file = new TFile("myFile.root", "RECREATE");

    TTree *tree = new TTree("myTree", "");

    tree->Branch("cal_time_var", &cal_time);

    string line;

    double pendenza = 1;
    double intercetta= 0;

    //cout << 0.1*pendenza+intercetta << endl;

    ifstream infile(inputfile.c_str());
    while(getline(infile,line))
    {
        stringstream ss(line);
        ss >> time;

        if(time>0.1){

            cal_time=(time-intercetta)/pendenza;
            tree->Fill(); 
        }

    }

    file->Write();

    infile.close();

    file->Close();

    /*
    struct RooFit_prova
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
    RooRealVar tau("tau", "mean lifetime",2.19, 1.4, 3.7);
    
   //RooRealVar tau_minus("tau_minus", "mean lifetime minus",0.163, 0.1, 0.2);

    //RooRealVar bl("bl", "baseline", 0.5, 1e-3, 50.);

    //RooRealVar theta("theta", "theta", 1, 1e-3, 10);

    RooRealVar lambda("lambda", "lambda", 1.7e-6, 1e-6, 2.5e-6);

    //RooRealVar slope("slope", "slope", -1., -1e-3, -1e3);

    //RooRealVar alpha("alpha", "alpha", 1. , 0.1, 10. );

    //fit pdfs
    RooPolynomial pol0("pol0", "pol0", cal_time_var);
    
    
    //RooPolynomial pol1("pol1", "pol1", cal_time_var, slope);

    RooGenericPdf bkg("bkg", "lambda*exp(-lambda*cal_time_var)", RooArgSet(cal_time_var,lambda));

    RooGenericPdf exp("exp", "(1/tau)*exp(-cal_time_var*(1/tau))", RooArgSet(cal_time_var,tau));

    //RooGenericPdf gamma("gamma", "pow(lambda,alpha)/tgamma(alpha)*pow(cal_time_var, alpha-1)*exp(-lambda*cal_time_var)", RooArgSet(cal_time_var,lambda,alpha));

    //RooGenericPdf exp_tot("exp_tot", "0.445*(1/tau_minus)*exp(-cal_time_var*(1/tau_minus)) + 0.555*(1/tau)*exp(-cal_time_var*(1/tau))", RooArgSet(cal_time_var,tau_minus,tau));
    //RooGenericPdf exp_minus("exponential_minus", "(1/tau_minus)*exp(-cal_time_var*(1/tau_minus))", RooArgSet(cal_time_var,tau_minus));

    //fraction of the first pdf (second one's fraction is 1-fraction)

    

    RooRealVar  fraction("fraction", "fraction", 0.5,0,1);
    

    //RooGenericPdf real("real", "real", "fraction*1/2.47*exp(-1/2.47*cal_time_var)", RooArgSet(cal_time_var, fraction));

    RooAddPdf model("model", "model", RooArgList(exp, pol0), RooArgList(fraction));

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooMsgService::instance().setSilentMode(true);

    const int n=3;

    double tau_values[n];
    double tau_val_errs[n];
    double starts[n];
    double xerrors[n];

    double start=0;
    double tau_result=0;
    double tau_result_error=0;

    for(int i =0; i<n; i++)
    {

    start=0.5+i*0.5;

    RooFitResult* r= model.fitTo(data, Range(start,15), Save(), PrintLevel(-1), PrintEvalErrors(-1));
    TCanvas *c5 = new TCanvas();
    r->Print();

    tau_result = tau.getVal();
    tau_result_error = tau.getError();
    //cout << tau_result << "  " << tau_result_error << endl;

    starts[i]=start;
    tau_values[i]=tau_result;
    tau_val_errs[i]=tau_result_error;
    xerrors[i]=0;    

    RooPlot *cal_time_frame = cal_time_var.frame();

    data.plotOn(cal_time_frame, Name("data"));
    model.plotOn(cal_time_frame, Name("Model"));
    model.plotOn(cal_time_frame, Components("exp"), LineColor(kGreen));
    model.plotOn(cal_time_frame, Components("pol0"), LineStyle(kDashed), LineColor(kMagenta));

    exp .paramOn(cal_time_frame);

    TLegend *legend = new TLegend();
    legend->AddEntry(cal_time_frame->getObject(0), "Data", "p");
    legend->AddEntry(cal_time_frame->getObject(1), "Total fit", "l");
    legend->AddEntry(cal_time_frame->getObject(2), "#mu exponential", "l");
    legend->AddEntry(cal_time_frame->getObject(3), "Background", "l");
    

   //frame_pull->GetYaxis()->SetTitle("Pulls");

   cout << "chi^2 = " << cal_time_frame->chiSquare() << endl;

      // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
   // -------------------------------------------------------
 
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hresid = cal_time_frame->residHist("data", "Model");
 
   // Construct a histogram with the pulls of the data w.r.t the curve
   RooHist *hpull = cal_time_frame->pullHist("data", "Model");
 
   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame2 = cal_time_var.frame(Title("Residual Distribution"), Bins(10));
   frame2->addPlotable(hresid, "P");
 
   // Create a new frame to draw the pull distribution and add the distribution to the frame
   RooPlot *frame3 = cal_time_var.frame(Title("Pull Distribution"), Bins(10));
   frame3->addPlotable(hpull, "P");

    TCanvas *c1 = new TCanvas();
    c1->Divide(1,2);
    c1->cd(1);
    cal_time_frame->Draw();
    legend->Draw();
    c1->cd(2);
    frame3->Draw();
   
    //new TCanvas();
    //frame2->Draw();
    //new TCanvas();
    //frame3->Draw();
    }

    TCanvas* Finale=new TCanvas();
    Finale->cd();
    TGraphErrors* graph= new TGraphErrors(n, starts, tau_values, xerrors, tau_val_errs);
    graph->Draw("APE");
    graph->SetMarkerStyle(20);
    graph->GetXaxis()->SetTitle("fit range start values");
    graph->GetYaxis()->SetTitle("tau");

}




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

void G2fit(string inputfile){

    double time;
    double cal_time;

    TFile *file = new TFile("myFile.root", "RECREATE");

    TTree *tree = new TTree("myTree", "");

    tree->Branch("cal_time_var", &cal_time);

    string line;

    double pendenza = 1;
    double intercetta= 0;

    //cout << 0.1*pendenza+intercetta << endl;

    ifstream infile(inputfile.c_str());
    while(getline(infile,line))
    {
        stringstream ss(line);
        ss >> time;

        if(time>0.1){

            cal_time=(time-intercetta)/pendenza;
            tree->Fill(); 
        }

    }

    file->Write();

    infile.close();

    file->Close();

    /*
    struct RooFit_prova
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
    RooRealVar tau("tau", "mean lifetime",2.19, 1.4, 3.7);
    RooRealVar polar("polar", "polarizazioje",0, 0.46, 1);
    RooRealVar pulsazione("pulsazione", "pulsazione",3.5, 4.6, 5.5);
    RooRealVar phi("phi", "fase",0, 0.1, 6.8);
    
   //RooRealVar tau_minus("tau_minus", "mean lifetime minus",0.163, 0.1, 0.2);

    //RooRealVar bl("bl", "baseline", 0.5, 1e-3, 50.);

    //RooRealVar theta("theta", "theta", 1, 1e-3, 10);

    RooRealVar lambda("lambda", "lambda", 1.7e-6, 1e-6, 2.5e-6);

    //RooRealVar slope("slope", "slope", -1., -1e-3, -1e3);

    //RooRealVar alpha("alpha", "alpha", 1. , 0.1, 10. );

    //fit pdfs
    RooPolynomial pol0("pol0", "pol0", cal_time_var);
    
    
    //RooPolynomial pol1("pol1", "pol1", cal_time_var, slope);

    RooGenericPdf bkg("bkg", "lambda*exp(-lambda*cal_time_var)", RooArgSet(cal_time_var,lambda));

    RooGenericPdf exp("exp", "(1/tau) * exp(-cal_time_var / tau) * (1 + polar * sin(pulsazione * cal_time_var + phi))", 
    RooArgSet(cal_time_var, tau, polar, pulsazione, phi));
    //RooGenericPdf gamma("gamma", "pow(lambda,alpha)/tgamma(alpha)*pow(cal_time_var, alpha-1)*exp(-lambda*cal_time_var)", RooArgSet(cal_time_var,lambda,alpha));

    //RooGenericPdf exp_tot("exp_tot", "0.445*(1/tau_minus)*exp(-cal_time_var*(1/tau_minus)) + 0.555*(1/tau)*exp(-cal_time_var*(1/tau))", RooArgSet(cal_time_var,tau_minus,tau));
    //RooGenericPdf exp_minus("exponential_minus", "(1/tau_minus)*exp(-cal_time_var*(1/tau_minus))", RooArgSet(cal_time_var,tau_minus));

    //fraction of the first pdf (second one's fraction is 1-fraction)

    

    RooRealVar  fraction("fraction", "fraction", 0.5,0,1);
    

    //RooGenericPdf real("real", "real", "fraction*1/2.47*exp(-1/2.47*cal_time_var)", RooArgSet(cal_time_var, fraction));

    RooAddPdf model("model", "model", RooArgList(exp, pol0), RooArgList(fraction));

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooMsgService::instance().setSilentMode(true);

    const int n=10;

    double tau_values[n];
    double tau_val_errs[n];
    double starts[n];
    double xerrors[n];

    double start=0;
    double tau_result=0;
    double tau_result_error=0;

    for(int i =0; i<n; i++)
    {

    start=0.5+i*0.5;

    RooFitResult* r= model.fitTo(data, Range(start,15), Save(), PrintLevel(-1), PrintEvalErrors(-1));
    TCanvas *c5 = new TCanvas();
    r->Print();

    tau_result = tau.getVal();
    tau_result_error = tau.getError();
    //cout << tau_result << "  " << tau_result_error << endl;

    starts[i]=start;
    tau_values[i]=tau_result;
    tau_val_errs[i]=tau_result_error;
    xerrors[i]=0;    

    RooPlot *cal_time_frame = cal_time_var.frame();

    data.plotOn(cal_time_frame, Name("data"));
    model.plotOn(cal_time_frame, Name("Model"));
    model.plotOn(cal_time_frame, Components("exp"), LineColor(kGreen));
    model.plotOn(cal_time_frame, Components("pol0"), LineStyle(kDashed), LineColor(kMagenta));

    exp .paramOn(cal_time_frame);

    TLegend *legend = new TLegend();
    legend->AddEntry(cal_time_frame->getObject(0), "Data", "p");
    legend->AddEntry(cal_time_frame->getObject(1), "Total fit", "l");
    legend->AddEntry(cal_time_frame->getObject(2), "#mu decay", "l");
    legend->AddEntry(cal_time_frame->getObject(3), "Background", "l");
    

   //frame_pull->GetYaxis()->SetTitle("Pulls");

   cout << "chi^2 = " << cal_time_frame->chiSquare() << endl;

      // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
   // -------------------------------------------------------
 
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hresid = cal_time_frame->residHist("data", "Model");
 
   // Construct a histogram with the pulls of the data w.r.t the curve
   RooHist *hpull = cal_time_frame->pullHist("data", "Model");
 
   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame2 = cal_time_var.frame(Title("Residual Distribution"), Bins(10));
   frame2->addPlotable(hresid, "P");
 
   // Create a new frame to draw the pull distribution and add the distribution to the frame
   RooPlot *frame3 = cal_time_var.frame(Title("Pull Distribution"), Bins(10));
   frame3->addPlotable(hpull, "P");

    TCanvas *c1 = new TCanvas();
    c1->Divide(1,2);
    c1->cd(1);
    cal_time_frame->Draw();
    legend->Draw();
    c1->cd(2);
    frame3->Draw();
   
    //new TCanvas();
    //frame2->Draw();
    //new TCanvas();
    //frame3->Draw();
    }

    TCanvas* Finale=new TCanvas();
    Finale->cd();
    TGraphErrors* graph= new TGraphErrors(n, starts, tau_values, xerrors, tau_val_errs);
    graph->Draw("APE");
    graph->SetMarkerStyle(20);
    graph->GetXaxis()->SetTitle("fit range start values");
    graph->GetYaxis()->SetTitle("tau");

}