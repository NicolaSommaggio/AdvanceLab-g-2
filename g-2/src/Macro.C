#include <TF1.h>
#include <TSpectrum.h>
#include <fstream>
#include <TMath.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <string>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TStyle.h>


#include <stdint.h>
//typedef  char int8_t;
// typedef  short int16_t;
// typedef  unsigned short uint16_t;
// typedef  int int32_t;
// typedef  unsigned int uint32_t;
// typedef  unsigned long long uint64_t;
// typedef  long long int64_t;

#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "InfoAcq.cc"
#include "Event.cc"

/*
Function descriptions:
1) adc_to_mv converts the data from arbitrary units to mV
2) HISTO plot Histogram for Gain
3)Read Tree
4) class FitFunction NOT USED JUST IGNORE IT
5) fitfun model to fit poisson
6) fitfunexp another model just used to test NOT USED
7) FindMaxinumSiPM: find the gain of a given dataset
8) FindMaximumSiPMEXP: the same using fitfunexp NOT USED
9) findmaxlight: a lighter and faster version of FindMaximumSiPM
10)FitGain: linear fit of the GAIN
11)DarkCounts: Prototype to test the darkcounts
12)DarkCountslight: final and faster version of DarkCounts
13)grafico: Plot the gainfit and the darkcountrate for a given SiPM

*/


float adc_to_mv(int16_t raw, int16_t rangeIndex, int16_t maxADCValue)
{
	uint16_t inputRanges [12] = {
			10,
			20,
			50,
			100,
			200,
			500,
			1000,
			2000,
			5000,
			10000,
			20000,
			50000};

	return (raw * inputRanges[rangeIndex])*1. / maxADCValue;
}

TH1F* HISTO(const char *fileName, bool negative, int binNumber)
{

	// dichiaro le struct
	InfoAcq::chSettings chSet1;
	InfoAcq::chSettings chSet2;
	InfoAcq::samplingSettings sampSet;
	
	// dichiaro le variabili dell'evento
/*	uint64_t ID;
	uint32_t samplesStored;
	int64_t triggerInstant;
	int16_t timeUnit;
	int16_t* sample;

	uint64_t waveformInBlock;
	uint64_t elapsedTime;
	uint64_t waveformStored;
*/
	unsigned long long ID;
	int samplesStored;
	long long triggerInstant;
	short timeUnit;
	short* sample;


	// apro il file in sola lettura
	TFile *input_file = new TFile(fileName,"READ");

	// leggo i trees
	TTree *treeCh = (TTree*)input_file->Get("Channels");
	TTree *treeSamp = (TTree*)input_file->Get("SampSets");
	TTree *treeEvt = (TTree*)input_file->Get("Event");
	// TFile->Get() restituisce un oggetto generico che va
	// convertito esplicitamente anteponendo (TTree*)

	// prelevo i branch con le info e li associo alle struct
	treeCh->SetBranchAddress("Ch1",&chSet1.enabled);
	treeCh->SetBranchAddress("Ch2",&chSet2.enabled);
	treeSamp->SetBranchAddress("Settings",&sampSet.max_adc_value);

	// leggo le entries
	// dichiaro l'oggetto InfoAcq e lo riempio
	InfoAcq* info = new InfoAcq();
	cout << "Riempio l'oggetto INFO\n";
	treeCh->GetEntry(0);
	treeSamp->GetEntry(0);
	info->FillSettings(&chSet1,&chSet2,&sampSet);

	// imposto i branches per gli eventi
	sample = new short[sampSet.samplesStoredPerEvent];
	treeEvt->SetBranchAddress("ID",&ID);
	treeEvt->SetBranchAddress("nSamp",&samplesStored);
	treeEvt->SetBranchAddress("Instant",&triggerInstant);
	treeEvt->SetBranchAddress("TimeUnit",&timeUnit);
	treeEvt->SetBranchAddress("WaveformsA",&sample[0]);


	Long64_t nEvt = treeEvt->GetEntries();
	float maximum = 0.0;
	float minimum = 0.0;
	
	// spettro in energia
	float xmin= negative? adc_to_mv(sampSet.max_adc_value,chSet1.range,-1*sampSet.max_adc_value) : 0 ; 
	float xmax = negative? 0 : adc_to_mv(sampSet.max_adc_value,chSet1.range,sampSet.max_adc_value) ;
	TH1F* spectrumMaximum = new TH1F( "hMax", "Maxima Distribution Spectrum", binNumber, xmin,xmax );
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	TH1F* spectrumMaximum = new TH1F( "hMax", "Maximum Spectrum", 256, xmin,xmax );
	cerr<<" Range   : "<<xmax-xmin<<" mV "<<endl;
	cerr<<" #Events : "<<nEvt<<endl;
	cerr<<" #Samples: "<<sampSet.samplesStoredPerEvent<<endl;
	cerr<<" Timestamp: "<<sampSet.timeIntervalNanoseconds<<" ns"<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (Long64_t index=0; index<nEvt; index++) 
	{
		treeEvt->GetEntry(index);
		for (int ii=0; ii<sampSet.samplesStoredPerEvent; ii++)
		{
			float value =  adc_to_mv(sample[ii],chSet1.range,sampSet.max_adc_value);
			if (value > maximum) maximum = value ;
			if (value < minimum) minimum = value ;
			
		}
		spectrumMaximum->Fill(negative?minimum-maximum:maximum-minimum);
		maximum = 0.0;
		minimum = 0.0;
	}
	return spectrumMaximum;
}


double FindEntries(const char *filename, int Nbins, double xmin, double xmax){

HISTO(filename, false , Nbins );
TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");	

int firstbin = hMax->FindBin(xmin);
int secondbin = hMax->FindBin(xmax);

double integral=hMax->Integral(firstbin,secondbin);

cout << TMath::Sqrt(integral) << endl;

cout << "efficiency    " << (hMax->GetEntries()-integral)/(hMax->GetEntries())*100 << "%" << endl;

hMax->Draw();
return integral;

}

void linear(const char *filename, int N, double xmin, double xmax){

    ifstream file(filename);
    if (file.fail()){
        cout << "the file does not open correctly!" << endl;
    }

    double x[N];
    double y[N];
    double yerror[N];
    double xerror[N];

    for(int i=0; i<N; i++){
      xerror[i]=0;
    }

    int i=0;
    while(true){
      file >> x[i] >> y[i] >> yerror[i];
      i++;
      if ( file.eof() ) break;
    }

    //file.close();

    TCanvas *c1= new TCanvas();

    TGraphErrors *graph = new TGraphErrors (N,x,y,xerror,yerror); //Npoints, x ,y , ex, ey

    	// Define a function which fits the points
    	TF1 *fitfun = new TF1("fitfun","[0]+[1]*x",xmin,xmax);
      fitfun->SetParNames("Intercept", "Slope");
      fitfun->SetParameter(0,0);
      fitfun->SetParameter(1,1);
    	// Fit

    	graph->Fit(fitfun, "RQN+");

		auto frame = c1->DrawFrame(-10, -10, 310, 800);
		frame->GetXaxis()->SetTitle("displacement [mm]");
		frame->GetYaxis()->SetTitle("no signal counts");
		graph->SetMarkerColor(kAzure-3);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(1.1);
		graph->Draw("p");
		fitfun->Draw("SAME");

    	// Get the parameters
      double m, q, sigma_m, sigma_q;
    	m = fitfun->GetParameter(1);
    	q = fitfun->GetParameter(0);
      sigma_m= fitfun->GetParError(1);
      sigma_q= fitfun->GetParError(0);

    	cout << "Slope = " << m << " +/- " << sigma_m << endl << "Intercept = " << q << " +/- " << sigma_q << endl;
      //auto st1 = new TStyle("st1","my style");
      //st1->SetOptFit(0111);
      //st1->cd();


}


void ReadTree(const char *fileName, bool negative)
{

	// dichiaro le struct
	InfoAcq::chSettings chSet1;
	InfoAcq::chSettings chSet2;
	InfoAcq::samplingSettings sampSet;
	
	// dichiaro le variabili dell'evento
/*	uint64_t ID;
	uint32_t samplesStored;
	int64_t triggerInstant;
	int16_t timeUnit;
	int16_t* sample;

	uint64_t waveformInBlock;
	uint64_t elapsedTime;
	uint64_t waveformStored;
*/
	unsigned long long ID;
	int samplesStored;
	long long triggerInstant;
	short timeUnit;
	short* sample;


	// apro il file in sola lettura
	TFile *input_file = new TFile(fileName,"READ");

	// leggo i trees
	TTree *treeCh = (TTree*)input_file->Get("Channels");
	TTree *treeSamp = (TTree*)input_file->Get("SampSets");
	TTree *treeEvt = (TTree*)input_file->Get("Event");
	// TFile->Get() restituisce un oggetto generico che va
	// convertito esplicitamente anteponendo (TTree*)

	// prelevo i branch con le info e li associo alle struct
	treeCh->SetBranchAddress("Ch1",&chSet1.enabled);
	treeCh->SetBranchAddress("Ch2",&chSet2.enabled);
	treeSamp->SetBranchAddress("Settings",&sampSet.max_adc_value);

	// leggo le entries
	// dichiaro l'oggetto InfoAcq e lo riempio
	InfoAcq* info = new InfoAcq();
	cout << "Riempio l'oggetto INFO\n";
	treeCh->GetEntry(0);
	treeSamp->GetEntry(0);
	info->FillSettings(&chSet1,&chSet2,&sampSet);

	// imposto i branches per gli eventi
	sample = new short[sampSet.samplesStoredPerEvent];
	treeEvt->SetBranchAddress("ID",&ID);
	treeEvt->SetBranchAddress("nSamp",&samplesStored);
	treeEvt->SetBranchAddress("Instant",&triggerInstant);
	treeEvt->SetBranchAddress("TimeUnit",&timeUnit);
	treeEvt->SetBranchAddress("WaveformsA",&sample[0]);


	Long64_t nEvt = treeEvt->GetEntries();
	float maximum = 0.0;
	float minimum = 0.0;
	
	// spettro in energia
	float xmin= negative? adc_to_mv(sampSet.max_adc_value,chSet1.range,-1*sampSet.max_adc_value) : 0 ; 
	float xmax = negative? 0 : adc_to_mv(sampSet.max_adc_value,chSet1.range,sampSet.max_adc_value) ;
	TH1F* spectrumMaximum = new TH1F( "hMax", "Maxima Distribution Spectrum", 128, xmin,xmax );
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	TH1F* spectrumMaximum = new TH1F( "hMax", "Maximum Spectrum", 256, xmin,xmax );
	cerr<<" Range   : "<<xmax-xmin<<" mV "<<endl;
	cerr<<" #Events : "<<nEvt<<endl;
	cerr<<" #Samples: "<<sampSet.samplesStoredPerEvent<<endl;
	cerr<<" Timestamp: "<<sampSet.timeIntervalNanoseconds<<" ns"<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (Long64_t index=0; index<nEvt; index++) 
	{
		treeEvt->GetEntry(index);
		for (int ii=0; ii<sampSet.samplesStoredPerEvent; ii++)
		{
			float value =  adc_to_mv(sample[ii],chSet1.range,sampSet.max_adc_value);
			if (value > maximum) maximum = value ;
			if (value < minimum) minimum = value ;
		}
		spectrumMaximum->Fill(negative?minimum-maximum:maximum-minimum);
		maximum = 0.0;
		minimum = 0.0;
	}
	
}

class FitFunction
{
	private:
		int n;
  //  float _bw;
	public:
  //  set_bw(float bw){_bw=bw;}
	FitFunction(int n): n(n){};
	double operator()(double* x, double *p)
	{
		//const int n = 5; // n gaussiane
		// 2 + 2n parametri

		double offset = p[0];
		double gain = p[1];
		double norm = p[3];
		double mu = p[2];
		//double *area = &p[2]; // n ampiezze
		double *sigma = &p[4]; // n sigma
		double X = x[0];
		double result = 0;
		cout<<"called"<<endl;
		for (int i =0; i < n; i++)
		{
			result += norm*TMath::Poisson(i, mu)*TMath::Gaus(X, offset-i*gain, sigma[i], kTRUE);
		}
		return result;
	}
};


double fitfun(double* x, double *p)
	{
		//const int n = 5; // n gaussiane
		// 2 + 2n parametri

		double offset = p[0];
		double gain = p[1];
		double norm = p[3];
		double mu = p[2];
		//double *area = &p[2]; // n ampiezze
		double *sigma = &p[4]; // n sigma
		double X = x[0];
		double result = 0;
		int n=4;
		for (int i =0; i < n; i++)
		{
			result += norm*TMath::Poisson(i, mu)*TMath::Gaus(X, offset-i*gain, sigma[i], kTRUE);
		}
		return result;
	}
	
double fitfunexp(double* x, double *p)
	{
		//const int n = 5; // n gaussiane
		// 2 + 2n parametri

		double offset = p[0];
		double gain = p[1];
		double norm = p[3];
		double mu = p[2];
		//double *area = &p[2]; // n ampiezze
		double *sigma = &p[4]; // n sigma
		double X = x[0];
		double result = 0;
		int n=4;
		for (int i =0; i < n; i++)
		{
			result += norm*TMath::Exp(mu*-1*i)*TMath::Gaus(X, offset-i*gain, sigma[i], kTRUE);
		}
		return result;
	}

int findMaximumSiPM(string filename1, string filename2, string filename3 /* ,int nsigma*/, int nPeaks=4)
{

string filename = filename1+filename2+filename3;

float offset,  eOffset, gain,  eGain;

//	gSystem->Load("ReadTree_C.so");
//	gROOT->ProcessLineSync(".L ReadTree.C+");
	ReadTree(filename.c_str(), true);
	TSpectrum t;
	TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
	hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");
	if (hMax == 0) return 1;
//
	t.Search(hMax,2,"",0.01);

	//	int nPeaks = t.GetNPeaks() + 1;

	if (nPeaks <= 1) return 1;
	//FitFunction* ff = new FitFunction(nPeaks);
	double xmax = hMax->GetXaxis()->GetXmax();
	double xmin = hMax->GetXaxis()->GetXmin() + (xmax - hMax->GetXaxis()->GetXmin())/90.0;

	//TF1* f = new TF1("fitFun", ff, xmin, xmax, nPeaks + 4, "ff"); 
	TF1* f = new TF1("fitfun", fitfun, xmin, xmax, nPeaks + 4); 
	f->SetNpx(1000);

	f->SetParNames("offset","gain", "mu", "norm");
	float mu=t.GetPositionY()[1]/t.GetPositionY()[0];
	if (mu>1) mu=1;
	f->SetParameters(t.GetPositionX()[0], 
			 t.GetPositionX()[0] - t.GetPositionX()[1],
			 mu, 
			 hMax->Integral("width"));
			 //			 hMax->Integral(2,90,"width"));
	cout << t.GetPositionX()[0] << " "<< t.GetPositionX()[0] - t.GetPositionX()[1] << " "<< 
		t.GetPositionX()[0]/t.GetPositionX()[1] << " "<< hMax->Integral("width") << endl;
	f->SetParLimits(0, -40, 0);
	f->SetParLimits(1, 1, 40);
	f->SetParLimits(2, 0, 1);
	f->SetParLimits(3, hMax->Integral("width")/3.0, hMax->Integral("width")*3.0);
	for (int i = 0; i < nPeaks; i++)
	{
		f->SetParameter(4+i, 1+0.15*i);
		f->SetParLimits(4+i, 0.05,4);
	}


	hMax->Fit(f, "EMR+");
	TCanvas c;
	hMax->Draw();
	string imgname = "graphs/"+filename1+filename2  + ".pdf";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/"+filename1+filename2  + ".C";
	c.SaveAs(imgname.c_str());
/*
	imgname = "graphs/" + filename + ".tex";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/" + filename + ".C";
	c.SaveAs(imgname.c_str());
*/

	offset = f->GetParameter(0);
	eOffset = f->GetParError(0);
	gain = f->GetParameter(1);
	eGain = f->GetParError(1);
	//file << f->GetParameter(0) << " " << f->GetParError(0) << " " << f->GetParameter(1) << " " << f->GetParError(1) << endl;

	// *** New ***
	cerr<< "### From Fit : "<<endl;
	cerr<< "\n *** Gain   : "<<gain<<" +- "<<eGain<<endl;
	float mu_fit = f->GetParameter(2);
	float sigma_mu_fit = f->GetParError(2);
	cerr<<   " *** mu_fit : "<<mu_fit<<" +- "<<sigma_mu_fit<<"\n\n"<<endl;

	float *areas = new float[nPeaks-1];
	for(int i=1;	i<nPeaks;	i++)
	{
	cerr<< " *** Peak "<<i<<endl;
	float sigma = f->GetParameter(3+i);
	cerr<< "   # Sigma : "<< sigma <<endl;
	float centroid = -1.*i*gain;

///////////////////////////////
// Half gain for integral
//////////////////////////////
	float half_gain = gain / 2.0;
        int binMin = hMax->GetXaxis()->FindBin(centroid - half_gain );
        int binMax = hMax->GetXaxis()->FindBin(centroid + half_gain );
	areas[i-1] = hMax->Integral(binMin,binMax);

/*
//////////////////////////////
// Nsigmas for integral
/////////////////////////////
        int binMin = hMax->GetXaxis()->FindBin(centroid - nsigma*sigma);
        int binMax = hMax->GetXaxis()->FindBin(centroid + nsigma*sigma);
	areas[i-1] = hMax->Integral(binMin,binMax);
	cerr<< "   # Area (from  "<<binMin<< " to "<<binMax<<" bins) : "<<areas[i-1]<<"\n"<<endl;
*/
	}


	cerr<<"	*** Half gain *** "<<endl;

// Poissonian mean value
/*
	double *mu_s = new double[nPeaks-2];
	double *sigma_mu_s = new double[nPeaks-2];
	cerr<< "\n\n *** Poissonian expectation value *** "<<endl;
	for(int i=2;	i<nPeaks;	i++)
	{
	float tmp_mu = (areas[i-1]/areas[i-2])*(float)i;
	float tmp_sigma_mu = (1./areas[i-2])* TMath::Sqrt( areas[i-1]* ( 1.+ (areas[i-1]/areas[i-2])  )	);
		mu_s[i-2] = tmp_mu;
		sigma_mu_s[i-2] = tmp_sigma_mu;
		cerr<< " /mu_{"<<(i-1)<<"} : "<< tmp_mu<<" +- "<<tmp_sigma_mu<<endl;
	}
	double sigma_mu_mean ;
	double mu_mean = weighted_mean_with_error(mu_s,sigma_mu_s,nPeaks-2,sigma_mu_mean);
*/

cerr<< "\n\n *** Poissonian expectation value *** "<<endl;
        float   mu_s = (areas[1]/areas[0])* 2.;
        float   sigma_mu_s = (1./areas[0])* TMath::Sqrt( areas[1]* ( 1.+ (areas[1]/areas[0])  ) );
                cerr<< " /mu  : "<< mu_s<<" +- "<<sigma_mu_s<<endl;

//ofstream fout; fout.open("Gain_log.txt",ios::app);
ofstream fout; fout.open("Mu_log.txt",ios::app);
//fout<<filename2<<'\t'<<gain<<'\t'<<eGain<<'\t'<<mu_mean<<'\t'<<sigma_mu_mean<<endl;
fout<<filename2<<'\t'<<gain<<'\t'<<eGain<<'\t'<<mu_s<<'\t'<<sigma_mu_s<<endl;
fout.close();
	return 0;
}

int findMaximumSiPMEXP(string filename1, string filename2, string filename3 /* ,int nsigma*/, int nPeaks=4)
{

string filename = filename1+filename2+filename3;

float offset,  eOffset, gain,  eGain;

//	gSystem->Load("ReadTree_C.so");
//	gROOT->ProcessLineSync(".L ReadTree.C+");
	ReadTree(filename.c_str(), true);
	TSpectrum t;
	TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
	hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");
	if (hMax == 0) return 1;
//
	t.Search(hMax,2,"",0.01);

	//	int nPeaks = t.GetNPeaks() + 1;

	if (nPeaks <= 1) return 1;
	//FitFunction* ff = new FitFunction(nPeaks);
	double xmax = hMax->GetXaxis()->GetXmax();
	double xmin = hMax->GetXaxis()->GetXmin() + (xmax - hMax->GetXaxis()->GetXmin())/90.0;

	//TF1* f = new TF1("fitFun", ff, xmin, xmax, nPeaks + 4, "ff"); 
	TF1* f = new TF1("fitfunexp", fitfun, xmin, xmax, nPeaks + 4); 
	f->SetNpx(1000);

	f->SetParNames("offset","gain", "mu", "norm");
	float mu=t.GetPositionY()[1]/t.GetPositionY()[0];
	if (mu>1) mu=1;
	f->SetParameters(t.GetPositionX()[0], 
			 t.GetPositionX()[0] - t.GetPositionX()[1],
			 mu, 
			 hMax->Integral("width"));
			 //			 hMax->Integral(2,90,"width"));
	cout << t.GetPositionX()[0] << " "<< t.GetPositionX()[0] - t.GetPositionX()[1] << " "<< 
		t.GetPositionX()[0]/t.GetPositionX()[1] << " "<< hMax->Integral("width") << endl;
	f->SetParLimits(0, -40, 0);
	f->SetParLimits(1, 1, 40);
	f->SetParLimits(2, 0, 1);
	f->SetParLimits(3, hMax->Integral("width")/3.0, hMax->Integral("width")*3.0);
	for (int i = 0; i < nPeaks; i++)
	{
		f->SetParameter(4+i, 1+0.15*i);
		f->SetParLimits(4+i, 0.05,4);
	}


	hMax->Fit(f, "EMR+");
	TCanvas c;
	hMax->Draw();
	string imgname = "graphs/"+filename1+filename2  + ".pdf";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/"+filename1+filename2  + ".C";
	c.SaveAs(imgname.c_str());
/*
	imgname = "graphs/" + filename + ".tex";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/" + filename + ".C";
	c.SaveAs(imgname.c_str());
*/

	offset = f->GetParameter(0);
	eOffset = f->GetParError(0);
	gain = f->GetParameter(1);
	eGain = f->GetParError(1);
	//file << f->GetParameter(0) << " " << f->GetParError(0) << " " << f->GetParameter(1) << " " << f->GetParError(1) << endl;

	// *** New ***
	cerr<< "### From Fit : "<<endl;
	cerr<< "\n *** Gain   : "<<gain<<" +- "<<eGain<<endl;
	float mu_fit = f->GetParameter(2);
	float sigma_mu_fit = f->GetParError(2);
	cerr<<   " *** mu_fit : "<<mu_fit<<" +- "<<sigma_mu_fit<<"\n\n"<<endl;

	float *areas = new float[nPeaks-1];
	for(int i=1;	i<nPeaks;	i++)
	{
	cerr<< " *** Peak "<<i<<endl;
	float sigma = f->GetParameter(3+i);
	cerr<< "   # Sigma : "<< sigma <<endl;
	float centroid = -1.*i*gain;

///////////////////////////////
// Half gain for integral
//////////////////////////////
	float half_gain = gain / 2.0;
        int binMin = hMax->GetXaxis()->FindBin(centroid - half_gain );
        int binMax = hMax->GetXaxis()->FindBin(centroid + half_gain );
	areas[i-1] = hMax->Integral(binMin,binMax);

/*
//////////////////////////////
// Nsigmas for integral
/////////////////////////////
        int binMin = hMax->GetXaxis()->FindBin(centroid - nsigma*sigma);
        int binMax = hMax->GetXaxis()->FindBin(centroid + nsigma*sigma);
	areas[i-1] = hMax->Integral(binMin,binMax);
	cerr<< "   # Area (from  "<<binMin<< " to "<<binMax<<" bins) : "<<areas[i-1]<<"\n"<<endl;
*/
	}


	cerr<<"	*** Half gain *** "<<endl;

// Poissonian mean value
/*
	double *mu_s = new double[nPeaks-2];
	double *sigma_mu_s = new double[nPeaks-2];
	cerr<< "\n\n *** Poissonian expectation value *** "<<endl;
	for(int i=2;	i<nPeaks;	i++)
	{
	float tmp_mu = (areas[i-1]/areas[i-2])*(float)i;
	float tmp_sigma_mu = (1./areas[i-2])* TMath::Sqrt( areas[i-1]* ( 1.+ (areas[i-1]/areas[i-2])  )	);
		mu_s[i-2] = tmp_mu;
		sigma_mu_s[i-2] = tmp_sigma_mu;
		cerr<< " /mu_{"<<(i-1)<<"} : "<< tmp_mu<<" +- "<<tmp_sigma_mu<<endl;
	}
	double sigma_mu_mean ;
	double mu_mean = weighted_mean_with_error(mu_s,sigma_mu_s,nPeaks-2,sigma_mu_mean);
*/

cerr<< "\n\n *** Poissonian expectation value *** "<<endl;
        float   mu_s = (areas[1]/areas[0])* 2.;
        float   sigma_mu_s = (1./areas[0])* TMath::Sqrt( areas[1]* ( 1.+ (areas[1]/areas[0])  ) );
                cerr<< " /mu  : "<< mu_s<<" +- "<<sigma_mu_s<<endl;

//ofstream fout; fout.open("Gain_log.txt",ios::app);
ofstream fout; fout.open("Mu_log.txt",ios::app);
//fout<<filename2<<'\t'<<gain<<'\t'<<eGain<<'\t'<<mu_mean<<'\t'<<sigma_mu_mean<<endl;
fout<<filename2<<'\t'<<gain<<'\t'<<eGain<<'\t'<<mu_s<<'\t'<<sigma_mu_s<<endl;
fout.close();
	return 0;
}



std::array<double, 4>  findmaxlight(string filename1, int nPeaks=4)
{

string filename = filename1;

float offset,  eOffset, gain,  eGain;

	ReadTree(filename.c_str(), true);
	TSpectrum t;
	TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
	hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");
	
//
	t.Search(hMax,2,"",0.01);

	//	int nPeaks = t.GetNPeaks() + 1;


	//FitFunction* ff = new FitFunction(nPeaks);
	double xmax = hMax->GetXaxis()->GetXmax();
	double xmin = hMax->GetXaxis()->GetXmin() + (xmax - hMax->GetXaxis()->GetXmin())/90.0;

	//TF1* f = new TF1("fitFun", ff, xmin, xmax, nPeaks + 4, "ff"); 
	TF1* f = new TF1("fitfun", fitfun, xmin, xmax, nPeaks + 4); 
	f->SetNpx(1000);

	f->SetParNames("offset","gain", "mu", "norm");
	float mu=t.GetPositionY()[1]/t.GetPositionY()[0];
	if (mu>1) mu=1;
	f->SetParameters(t.GetPositionX()[0], 
			 t.GetPositionX()[0] - t.GetPositionX()[1],
			 mu, 
			 hMax->Integral("width"));
			 //			 hMax->Integral(2,90,"width"));
	cout << t.GetPositionX()[0] << " "<< t.GetPositionX()[0] - t.GetPositionX()[1] << " "<< 
		t.GetPositionX()[0]/t.GetPositionX()[1] << " "<< hMax->Integral("width") << endl;
	f->SetParLimits(0, -40, 0);
	f->SetParLimits(1, 1, 40);
	f->SetParLimits(2, 0, 1);
	f->SetParLimits(3, hMax->Integral("width")/3.0, hMax->Integral("width")*3.0);
	for (int i = 0; i < nPeaks; i++)
	{
		f->SetParameter(4+i, 1+0.15*i);
		f->SetParLimits(4+i, 0.05,4);
	}


	hMax->Fit(f, "EMRQ+");
	TCanvas c;
	hMax->Draw();

	offset = f->GetParameter(0);
	eOffset = f->GetParError(0);
	gain = f->GetParameter(1);
	eGain = f->GetParError(1);

	float mu_fit = f->GetParameter(2);
	float sigma_mu_fit = f->GetParError(2);

	std::array<double, 4> out;

	out[0]=offset;
	out[1]=eOffset;
	out[2]=gain;
	out[3]=eGain;

	return out;
}

double Landau(const char *nomefile, int Nbins)
{

float location, elocation, sigma, esigma, norm, mpv; //in the root landau function mu is not the mpv

	//ReadTree(filename.c_str(), false);  
	TSpectrum t;
	HISTO(nomefile, false , Nbins );
	TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
	hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");

	t.Search(hMax,2,"",0.01);

	//	int nPeaks = t.GetNPeaks() + 1;


	//FitFunction* ff = new FitFunction(nPeaks);
	double xmax = hMax->GetXaxis()->GetXmax();
	double xmin = hMax->GetXaxis()->GetXmin() + (xmax - hMax->GetXaxis()->GetXmin())/90.0;
	
	//TF1* f = new TF1("fitFun", ff, xmin, xmax, nPeaks + 4, "ff"); 
	TF1* f = new TF1("fitfun", "[0]*TMath::Landau(x,[1],[2])", xmin, xmax); 
	f->SetNpx(1000);

	f->SetParNames("norm","location", "sigma");
	f->SetParameters(60.0,95.0,70.0); 
			 //			 hMax->Integral(2,90,"width"));

	hMax->Fit(f, "EMR+");
	hMax->Draw();

	location = f->GetParameter(1);
	elocation = f->GetParError(1);
	sigma = f->GetParameter(2);
	esigma = f->GetParError(2);
	norm = f-> GetParameter(0);

	mpv = f->GetMaximumX();
	double* fitres=(double*) new double[3];
	fitres[0]=norm;
	fitres[1]=location;
	fitres[2]=sigma;
	/*TArrow* arrow = new TArrow(mpv, 140, mpv, 0);
	arrow -> SetLineWidth(1);
	arrow -> SetFillColor(0);
	arrow -> Draw("SAME");*/

	return mpv;
}

int npeaks=15;
//definition of f for multi peaks (not reliable)
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0]*TMath::Landau(x[0],par[1],par[2]);
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+3]; // "height"
      Double_t mean  = par[3*p+4];
      Double_t sigma = par[3*p+5];

      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

void MultiPeaksLandaunotreliable(const char *filename, int Nbins) {
   npeaks = 15;

   HISTO(filename, false , Nbins );
   TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
   hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");

   Double_t par[3000];
   TF1 *f = new TF1("f",fpeaks,0,350,3+3*npeaks);
   f->SetNpx(1000);
   f->SetParameters(par);

   // Use TSpectrum to find the peak candidates
   TSpectrum *s = new TSpectrum(2*npeaks);
   Int_t nfound = s->Search(hMax,1,"",0.30);
   printf("Found %d candidate peaks to fit\n",nfound);
 
   // Loop on all found peaks
   par[0] = 100;
   par[1] = 100;
   par[2] = 100;
   npeaks = 0;
   Double_t *xpeaks;
   xpeaks = s->GetPositionX(); //The number of found peaks and their positions are written into the members fNpeaks and fPositionX
   for (int p=0;p<nfound;p++) {
      Double_t xp = xpeaks[p];
      Int_t bin = hMax->GetXaxis()->FindBin(xp);
      Double_t yp = hMax->GetBinContent(bin); // to find the "y" value of a peak to initialize the height value for the peak
      par[3*npeaks+3] = yp; // "height"
      par[3*npeaks+4] = xp; // "mean"
      par[3*npeaks+5] = 3; // "sigma"
      npeaks++;
   }
   printf("Now fitting: Be patient\n");
   TF1 *fit = new TF1("fit",fpeaks,0,350,3+3*npeaks);
   fit->SetParameters(par);
   fit->SetNpx(1000);
   hMax->Fit("fit", "EMRQ+");
   hMax->Draw();
   cout << "approx mpv=" << fit->GetParameter(1) << "+/-" << fit->GetParError(1) << endl;

   

}
Double_t fpeaksJDP(Double_t *x, Double_t *par) {
	Double_t result=0;
	Double_t Lnorm=par[0];
	Double_t Llocation=par[1];
	Double_t Lsigma=par[2];
	
   Double_t gain=par[3];
   Double_t offgain=par[4];
   Double_t sigmagrowth=par[5];
   double_t sigmaoffset=par[6];
   double_t npeaks=par[7];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t mean  = gain*p+offgain;
	  Double_t norm  = Lnorm*TMath::Landau(mean,Llocation,Lsigma);
      Double_t sigma = sigmagrowth*sqrt(p)+sigmaoffset;

      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}
void MultiPeaksLandau(const char *filename, int Nbins=300,double xmin=5,double xmax=60) {
   
   HISTO(filename, false , Nbins );
   TH1F* hMax = (TH1F*)gDirectory->FindObject("hMax");
   hMax->GetXaxis()->SetTitle("Pulse Height[mV] ");
   Double_t par[3000];
   TF1 *f = new TF1("f",fpeaksJDP,xmin,xmax,8);
   f->SetNpx(1000);
   f->SetParameters(par);

	TSpectrum* t=new TSpectrum();
	t->Search(hMax,1,"",0.1);
	int npeaks=0;
	for(int i=0;i<t->GetNPeaks();i++)
	{
		if(t->GetPositionX()[i]>xmin &&t->GetPositionX()[i]<xmax)
		npeaks=npeaks+1;
	}
   // Loop on all found peaks
   par[0] = 200;
   f->SetParLimits(0,100,300);
   par[1] = 70;
   f->SetParLimits(1,40,80);
   par[2] = 50;
   f->SetParLimits(2,10,80);
   par[3]= t->GetPositionX()[1]-t->GetPositionX()[0];
   f->SetParLimits(3,5,15);
   par[4]= t->GetPositionX()[0];
   f->SetParLimits(4,5,20);
   par[5]=0.1;
   f->SetParLimits(5,0,1);
   par[6]=1.8;
   f->SetParLimits(7,0.5,10);
   par[7]=npeaks;
   f->SetParLimits(7,5,25);
   printf("Now fitting: Be patient\n");
   
   f->SetParameters(par);
   f->SetNpx(1000);
TCanvas* c1=new TCanvas();
c1->cd();
   hMax->Fit("f", "EMR+");
   hMax->Draw();

}
   


TGraphErrors fitGain(string files,int mute=0)
{
	ifstream fin(files.c_str());
	ofstream debug("log.txt");
	string filename;
	vector<float> offset, eOffset, gain,eGain, V,eV;
	float tmpOffset, tmpEOffset, tmpGain, tmpEGain, tmpV;
	gStyle->SetOptFit();
	
	std::array<double, 4> data;
	int npeaks;
	//PREPARE THE FILES
	while(fin>> filename>>tmpV >> npeaks)
	{
		data=findmaxlight(filename,npeaks);
		debug << data[0] << ' ' << data[1] << endl;
		offset.push_back(data[0]);
		eOffset.push_back(data[1]);
		gain.push_back(data[2]);
		eGain.push_back(data[3]);
		V.push_back(tmpV);
		eV.push_back(0.04);
		
	}
	//LINEAR FIT
	TGraphErrors g = TGraphErrors(offset.size(),&V[0] , &gain[0], &eV[0], &eGain[0]);
	
	TF1 f = TF1("linfit", "pol1");
	g.Fit(&f, "EM");
	
	TCanvas c5= TCanvas();
	c5.cd();
	g.GetXaxis()->SetTitle("V");
	g.GetYaxis()->SetTitle("gain ");
	g.Draw("AP");
	
	
	//INVERTED FIT FOR THE BREAKDOWN
	TGraphErrors c =TGraphErrors(offset.size(),&gain[0] , &V[0], &eGain[0], &eV[0]);
	
	TF1 d =TF1("linfit", "pol1");
	c.Fit(&d, "EM");

	cout<<"Breakdown V:	"<<d.GetParameter(0)<<"	+-	"<<d.GetParError(0)<<endl;
	TLegend *leg=new TLegend(0.2,0.2,0.4,0.4);
	c5.SetGrid();

	//GRAPHICAL: move the legend in the top left positions
		gPad->Update();
		TPaveStats *st = (TPaveStats*)g.FindObject("stats");
		
		st->SetX1NDC(0.15); //new x start position
		st->SetX2NDC(0.5); //new x end position
		st->SetY1NDC(0.68); //new x start position
		st->SetY2NDC(0.88); //new x end position
		
		c5.Draw();
	
	

	if (mute==1) { c5.Close(); }

	return g;
}

double DarkCount(const char *fileName, bool negative=0, int nwave=0)
{

	// dichiaro le struct
	InfoAcq::chSettings chSet1;
	InfoAcq::chSettings chSet2;
	InfoAcq::samplingSettings sampSet;
	
	// dichiaro le variabili dell'evento
/*	uint64_t ID;
	uint32_t samplesStored;
	int64_t triggerInstant;
	int16_t timeUnit;
	int16_t* sample;

	uint64_t waveformInBlock;
	uint64_t elapsedTime;
	uint64_t waveformStored;
*/
	unsigned long long ID;
	int samplesStored;
	long long triggerInstant;
	short timeUnit;
	short* sample;


	// apro il file in sola lettura
	TFile *input_file = new TFile(fileName,"READ");

	// leggo i trees
	TTree *treeCh = (TTree*)input_file->Get("Channels");
	TTree *treeSamp = (TTree*)input_file->Get("SampSets");
	TTree *treeEvt = (TTree*)input_file->Get("Event");
	// TFile->Get() restituisce un oggetto generico che va
	// convertito esplicitamente anteponendo (TTree*)

	// prelevo i branch con le info e li associo alle struct
	treeCh->SetBranchAddress("Ch1",&chSet1.enabled);
	treeCh->SetBranchAddress("Ch2",&chSet2.enabled);
	treeSamp->SetBranchAddress("Settings",&sampSet.max_adc_value);

	// leggo le entries
	// dichiaro l'oggetto InfoAcq e lo riempio
	InfoAcq* info = new InfoAcq();
	cout << "Riempio l'oggetto INFO\n";
	treeCh->GetEntry(0);
	treeSamp->GetEntry(0);
	info->FillSettings(&chSet1,&chSet2,&sampSet);

	// imposto i branches per gli eventi
	sample = new short[sampSet.samplesStoredPerEvent];
	treeEvt->SetBranchAddress("ID",&ID);
	treeEvt->SetBranchAddress("nSamp",&samplesStored);
	treeEvt->SetBranchAddress("Instant",&triggerInstant);
	treeEvt->SetBranchAddress("TimeUnit",&timeUnit);
	treeEvt->SetBranchAddress("WaveformsA",&sample[0]);


	Long64_t nEvt = treeEvt->GetEntries();
	float maximum = 0.0;
	float minimum = 0.0;
	
	// spettro in energia
	float xmin= negative? adc_to_mv(sampSet.max_adc_value,chSet1.range,-1*sampSet.max_adc_value) : 0 ; 
	float xmax = negative? 0 : adc_to_mv(sampSet.max_adc_value,chSet1.range,sampSet.max_adc_value) ;
	TH1D* Waveform = new TH1D( "Waveform", "Waveform", sampSet.samplesStoredPerEvent, 0,sampSet.samplesStoredPerEvent );
	TH1D* Waveformtemp = new TH1D( "Waveformtemp", "Waveformtemp", sampSet.samplesStoredPerEvent, 0,sampSet.samplesStoredPerEvent );

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//	TH1F* spectrumMaximum = new TH1F( "hMax", "Maximum Spectrum", 256, xmin,xmax );
	cerr<<" Range   : "<<xmax-xmin<<" mV "<<endl;
	cerr<<" #Events : "<<nEvt<<endl;
	cerr<<" #Samples: "<<sampSet.samplesStoredPerEvent<<endl;
	cerr<<" Timestamp: "<<sampSet.timeIntervalNanoseconds<<" ns"<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	TSpectrum t;
	int nPeaks;
	double mean=0, dev=0, max=0, supermax=0,count=0;
	double mean2=0, dev2=0;
	int temp=0;
	double peaksvector[nEvt];
	double peaksvector2[nEvt];
	double Wave[sampSet.samplesStoredPerEvent];
	double Wavesuper[sampSet.samplesStoredPerEvent];
	for (Long64_t index=0; index<nEvt; index++) 
	{
		treeEvt->GetEntry(index);
		for (int ii=0; ii<sampSet.samplesStoredPerEvent; ii++)
		{
			float value =  adc_to_mv(sample[ii],chSet1.range,sampSet.max_adc_value);
			if (value<max)
			{
				max=value;
			}
			
			if(nwave==index)
			{
				Waveform->SetBinContent(ii,-1*value);
				if (value<supermax)
				{
				supermax=value;
				}
				
			}
				Wave[ii]=-1*value;
				if(ii>3)
				{
					if(Wave[ii]>5&&Wave[ii]>Wave[ii-1]&&Wave[ii-1]>Wave[ii-2]&&Wave[ii-2]>Wave[ii-3]&&Wave[ii-3]>Wave[ii-4]&&temp<0)
					{
						count=count+1;
						temp=20;
						if(nwave==index)
						{
							Wavesuper[ii]=100;
						}
					}
					else
					{
						if(nwave==index)
						{
							Wavesuper[ii]=0;
						}
					}
					
					
				}
				
				temp=temp-1;
				Waveformtemp->SetBinContent(ii,-1*value);
		}
		t.Search(Waveformtemp,20,"",-8/max);
		nPeaks = t.GetNPeaks();
		peaksvector[index]=nPeaks;
		peaksvector2[index]=count;
		mean2=mean2+count;
		mean=mean+peaksvector[index];
		Waveformtemp->Reset();
		max=0;
		count=0;
		
	}
	mean=mean/nEvt;
	mean2=mean2/nEvt;
	for(int i=0;i<nEvt;i++)
	{
		dev=dev+pow(mean-peaksvector[i],2);
		dev2=dev2+pow(mean2-peaksvector2[i],2);
	}
	dev=sqrt(dev/(nEvt-1))/sqrt(nEvt);
	dev2=sqrt(dev2/(nEvt-1))/sqrt(nEvt);
		t.Search(Waveform,20,"",-8/supermax);


	 nPeaks = t.GetNPeaks();

	Waveform->DrawClone();
	for(int i=0;i<sampSet.samplesStoredPerEvent;i++)
	{
		Waveformtemp->SetBinContent(i,Wavesuper[i]);
	}
	Waveformtemp->SetLineColor(2);
	Waveformtemp->DrawClone("same");
	double time=sampSet.timeIntervalNanoseconds*sampSet.samplesStoredPerEvent*pow(10,-9);
	cout<<mean<<"	"<<dev<<endl;
	cout<<mean2<<"	"<<dev2<<endl;
	cout<<time<<endl;
	cout<<mean/time<<"	"<<dev/time<<endl;
	cout<<mean2/time<<"	"<<dev2/time<<endl;
return nPeaks;
}

std::array <double,2> DarkCountlight(const char *fileName, bool negative=0)
{

	//UNTIL YOU READ THE STOP COMMENT THIS IS JUST A SLITLY MODIFIED VERSION OF READ FILE
	//
	//
	
	// dichiaro le struct
	InfoAcq::chSettings chSet1;
	InfoAcq::chSettings chSet2;
	InfoAcq::samplingSettings sampSet;
	
	// dichiaro le variabili dell'evento
/*	uint64_t ID;
	uint32_t samplesStored;
	int64_t triggerInstant;
	int16_t timeUnit;
	int16_t* sample;

	uint64_t waveformInBlock;
	uint64_t elapsedTime;
	uint64_t waveformStored;
*/
	unsigned long long ID;
	int samplesStored;
	long long triggerInstant;
	short timeUnit;
	short* sample;


	// apro il file in sola lettura
	TFile *input_file = new TFile(fileName,"READ");

	// leggo i trees
	TTree *treeCh = (TTree*)input_file->Get("Channels");
	TTree *treeSamp = (TTree*)input_file->Get("SampSets");
	TTree *treeEvt = (TTree*)input_file->Get("Event");
	// TFile->Get() restituisce un oggetto generico che va
	// convertito esplicitamente anteponendo (TTree*)

	// prelevo i branch con le info e li associo alle struct
	treeCh->SetBranchAddress("Ch1",&chSet1.enabled);
	treeCh->SetBranchAddress("Ch2",&chSet2.enabled);
	treeSamp->SetBranchAddress("Settings",&sampSet.max_adc_value);

	// leggo le entries
	// dichiaro l'oggetto InfoAcq e lo riempio
	InfoAcq* info = new InfoAcq();
	cout << "Riempio l'oggetto INFO\n";
	treeCh->GetEntry(0);
	treeSamp->GetEntry(0);
	info->FillSettings(&chSet1,&chSet2,&sampSet);

	// imposto i branches per gli eventi
	sample = new short[sampSet.samplesStoredPerEvent];
	treeEvt->SetBranchAddress("ID",&ID);
	treeEvt->SetBranchAddress("nSamp",&samplesStored);
	treeEvt->SetBranchAddress("Instant",&triggerInstant);
	treeEvt->SetBranchAddress("TimeUnit",&timeUnit);
	treeEvt->SetBranchAddress("WaveformsA",&sample[0]);


	Long64_t nEvt = treeEvt->GetEntries();
	float maximum = 0.0;
	float minimum = 0.0;
	
	// spettro in energia
	float xmin= negative? adc_to_mv(sampSet.max_adc_value,chSet1.range,-1*sampSet.max_adc_value) : 0 ; 
	float xmax = negative? 0 : adc_to_mv(sampSet.max_adc_value,chSet1.range,sampSet.max_adc_value) ;

	//STOP READ FILE
	//
	//
	int nPeaks;
	double count=0;
	double mean2=0, dev2=0;
	int temp=0;
	double peaksvector2[nEvt];
	double Wave[sampSet.samplesStoredPerEvent];

	

	//READ ALL FILE AND COUNT THE NUMBER OF EVENT (SIMULATED TRIGGER)
	for (Long64_t index=0; index<nEvt; index++) 
	{
		treeEvt->GetEntry(index);
		for (int ii=0; ii<sampSet.samplesStoredPerEvent; ii++)
		{
			float value =  adc_to_mv(sample[ii],chSet1.range,sampSet.max_adc_value);

				Wave[ii]=-1*value;
				if(ii>3)
				{
					if(Wave[ii]>5&&Wave[ii]>Wave[ii-1]&&Wave[ii-1]>Wave[ii-2]&&Wave[ii-2]>Wave[ii-3]&&Wave[ii-3]>Wave[ii-4]&&temp<0)
					{
						count=count+1;
						temp=20;
					}
						
					
				}
				
				temp=temp-1;
		}
		peaksvector2[index]=count;
		mean2=mean2+count;
		count=0;
		
	}
	//COMPUTE MEAN DEV AND RATE
	mean2=mean2/nEvt;
	for(int i=0;i<nEvt;i++)
	{
		dev2=dev2+pow(mean2-peaksvector2[i],2);
	}
	dev2=sqrt(dev2/(nEvt-1))/sqrt(nEvt);

	double time=sampSet.timeIntervalNanoseconds*sampSet.samplesStoredPerEvent*pow(10,-9);
	cout<<"Mean Events: "<<mean2<<"	+-	"<<dev2<<endl;
	cout<<"Time: "<<time<<endl;
	cout<<"Rate: "<<mean2/time<<"	+-	"<<dev2/time<<endl;
	std:array <double,2> rate={mean2/time,dev2/time};
return rate;
}


void grafico(string gainfile,string darkfile)
{
	ifstream fin(darkfile.c_str());
	string filename;
	vector<float> rate, erate, V,eV;
	float tmpV;
	gStyle->SetOptFit();
	
	std::array<double, 2> data;
	while(fin >> filename >> tmpV)
	{
		cout<<"1"<<endl;
		data=DarkCountlight(filename.c_str(),0);
		cout<<"2"<<endl;
		rate.push_back(data[0]);
		erate.push_back(data[1]);
		V.push_back(tmpV);
		eV.push_back(0.04);
		cout<<data[0]<<"	"<<data[1]<<endl;
	}
	cout<<"3"<<endl;
	TGraphErrors* dark = new TGraphErrors(V.size(),&V[0] , &rate[0], &eV[0], &erate[0]);

	TCanvas* canvas=new TCanvas("questo","quello");
	TGraphErrors gainG=fitGain(gainfile,1);
	canvas->cd();
	double maxrate=0;
	for(int i=0;i<rate.size();i++)
	{
		if(rate[i]>maxrate)
		{maxrate=rate[i];}
	}
	auto frame = canvas->DrawFrame(30.3, 5, 33.2, 20);
	frame->GetXaxis()->SetTitle("V");
	frame->GetYaxis()->SetTitle("gain (mV)");
	gainG.DrawClone("p");
	Float_t rightmax = 1.1*maxrate;
    Float_t scale = gPad->GetUymax()/rightmax;
	cout<<rightmax<<"	"<<gPad->GetUymax()<<endl;
	dark->SetLineColor(3);
	dark->Scale(scale);

	dark->Draw("same P");
	
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),
	0,rightmax,510,"+L");
   axis->SetLineColor(3);
   axis->SetLabelColor(3);
   axis->SetTitleColor(3);
   axis->SetTextSize(2);
   axis->SetTitle("Dark Count Rate (Hz)");
   axis->Draw();
   
/*g->GetYaxis()->SetTitle("V");
	g->GetXaxis()->SetTitle("gain ");
	g->Draw("AP");
	cout<<"Breakdown V:	"<<f->GetParameter(0)<<"	+-	"<<f->GetParError(0)<<endl;
	TLegend *leg=new TLegend(0.2,0.2,0.4,0.4);
	c5->SetGrid();

	gPad->Update();
	TPaveStats *st = (TPaveStats*)g->FindObject("stats");
	st->SetX1NDC(0.15); //new x start position
	st->SetX2NDC(0.5); //new x end position
	st->SetY1NDC(0.68); //new x start position
	st->SetY2NDC(0.88); //new x end position
	c5->Draw();

	string imgname = "graphs/gain.png";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/gain.tex";
	c.SaveAs(imgname.c_str());
	imgname = "graphs/gain.C";
	c.SaveAs(imgname.c_str());
	*/
	
}