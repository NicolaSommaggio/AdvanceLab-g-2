struct variabili
{
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

void plt(const char *namefile, int Nbin, double min, double max)

{

	TCanvas* c1 = new TCanvas();
	//c1->SetLogy();

	variabili data;
	TFile *infile = new TFile(namefile); 
	TTree *intree = (TTree*)infile->Get("acq_tree_0");

	TBranch *inbranch = intree->GetBranch("acq_ch0");

	inbranch->SetAddress(&data.timetag);

	int entries= inbranch->GetEntries();
	TH1D *h = new TH1D("histo", "", Nbin, min , max);


	//histogram filling
	for (int j=0; j<entries; j++) {
		inbranch->GetEntry(j);
		h->Fill(data.qlong);
	}
  

  h->Draw();
	//h->GetYaxis()->SetTitle("");
	//h->GetXaxis()->SetTitle("");
  //h->SetTitle("Calibrated spectrum Det1 Ch0");

}
