#include "tdrstyle.C"


int uniqueglbl=1;
std::pair<TH1D*,TH1D*> extractFit( TCanvas* hCanv, double fitLow, double fitUp, int etabin, double _syst0){
	TList* t = hCanv->GetListOfPrimitives();
	TH2D* h = (TH2D*) hCanv->GetPrimitive(t->At(1)->GetName());
	int nbinsX = h->GetNbinsX();	

	//declare variables to be used inside
	double diff,sig;
	double prob1,prob2;
	int hbinX;
	double sum2=0.;
	int i = etabin;
	//loop and get projection
	
	TH1D* proj = h->ProjectionY(("hprojY"+std::to_string(uniqueglbl)).c_str(),i,i);
		//loop through and adjust errors according to systematics
	for(int j=1; j<=proj->GetNbinsX(); j++){
		proj->SetBinError(j,  std::sqrt( proj->GetBinError(j)*proj->GetBinError(j) + _syst0*_syst0 ) );
	}

	proj->Fit("pol2","","",fitLow,fitUp);
	TF1* fit2 = proj->GetFunction("pol2");
		
	double X[1] = {0.};		
	prob2 = fit2->GetProb();
	TF1* fit=fit2;
		
	std::cout<<"pol2 P="<<prob2<<std::endl;		
	
	if( prob2 < 0.01){
		std::cout<<"Warning >> pol2 with P<1%"<<std::endl;
	}
	TH1D *hint = new TH1D(("hint"+std::to_string(uniqueglbl)).c_str(),"Fitted CI", 20, 3, 20);
   	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
	hint->SetFillColor(kRed);
	hint->SetMarkerStyle(0);
	uniqueglbl++;

	return std::make_pair(proj,hint);



}
TCanvas* getCanvas(TFile* f,std::string fname, std::string histPath ){
	f = TFile::Open(fname.c_str());
	f->cd(histPath.c_str());
	std::string name = (gDirectory->GetListOfKeys()->At(0)->GetName());
	return (TCanvas*) f->Get((histPath+name).c_str());
}
void plotextrap(int option){
setTDRStyle();
	TFile* f_data;
	TFile* f_mc;
	TCanvas* c_data;
	TCanvas* c_mc;
	double _syst0;
	std::string name1;
	std::string name2;
	std::string path1;
	std::string path2;
	std::string legtitle;
	path1 = "tpTree/Medium_pt_eta/fit_eff_plots/";
	path2 = "tpTree/Medium_pt_eta/fit_eff_plots/";
	
	//1-3 iso, 4-6 sip
	if(option == 1){
		//read in set of files 2016
		name1 = "./extrapfiles/TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2016 Isolation | Medium Id";
	}
	if(option == 2){
		//read in set of files 2017
		name1 = "./extrapfiles/TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2017 Isolation | Medium Id";
	}
	if(option == 3){
		//read in set of files 2018
		name1 = "./extrapfiles/TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2018 Isolation | Medium Id";
	}



	if(option == 4){
		//read in set of files 2016
		name1 = "./extrapfiles/TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2016 SIP | Med. Id #cap Iso.";
	}
	if(option == 5){
		//read in set of files 2017
		name1 = "./extrapfiles/TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2017 SIP | Med. Id #cap Iso.";
	}
	if(option == 6){
		//read in set of files 2018
		name1 = "./extrapfiles/TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root";
		name2 = "./extrapfiles/TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root";
		_syst0 = 0.0057;
		legtitle= "2018 SIP | Med. Id #cap Iso.";
	}

	
	c_data = getCanvas(f_data, name1, path1);
	c_mc = getCanvas(f_mc, name2, path2);
	std::pair<TH1D*,TH1D*> h_data_bar, h_data_end, h_mc_bar, h_mc_end;
	
	h_data_bar = extractFit( c_data, 3,50, 1, _syst0);
	h_data_end = extractFit( c_data, 3,50, 2, _syst0);
	
	h_mc_bar = extractFit( c_mc, 3,50,1,_syst0);
	h_mc_end = extractFit( c_mc, 3,50,2,_syst0);
	
	gStyle->SetHatchesSpacing(.4);
	gStyle->SetHatchesLineWidth(1);

	TCanvas* c1 = new TCanvas();
	

	h_data_bar.first->GetXaxis()->SetRangeUser(1.,50.);
	h_data_bar.first->SetMinimum(.75);
	h_data_bar.first->SetMaximum(1.);
	h_data_bar.first->GetXaxis()->SetTitleFont(42);
	h_data_bar.first->GetXaxis()->SetLabelFont(42);
	h_data_bar.first->GetXaxis()->SetTitle("p_{T} [GeV]");
	h_data_bar.first->GetYaxis()->SetTitle("#varepsilon");
	



	h_data_bar.first->Draw();

	c1->Update();
	TPaveStats* stats =(TPaveStats*)c1->GetPrimitive("stats");
	stats->SetName("h1stats");
	stats->SetY1NDC(.152);
	stats->SetY2NDC(.252);
	stats->SetX1NDC(0.96);
	stats->SetX2NDC(.71);
	stats->SetLineWidth(2);
	stats->SetLineColor(kBlue);
	//stats->Draw();

	
	//h_data_bar.first->Draw("SAMES");
	h_data_bar.second->Draw("SAMES E3");
	h_data_bar.second->SetFillStyle(3352);
	h_data_bar.first->SetMarkerStyle(20);	
	h_data_bar.first->GetFunction("pol2")->SetLineColor(kBlue);
	h_data_bar.first->GetFunction("pol2")->SetLineWidth(2);
	h_data_bar.second->SetFillColor(kBlue);
	h_data_bar.first->SetMarkerColor(kBlue);
	h_data_bar.first->SetMarkerSize(1.3);


///////////////////////////////////////////////////////////////////
	
	h_mc_bar.first->Draw("SAMES");

	c1->Update();
	TPaveStats* stats2 =(TPaveStats*)c1->GetPrimitive("stats");
	stats2->SetName("h2stats");
	stats2->SetY1NDC(.255);
	stats2->SetY2NDC(.355);
	stats2->SetX1NDC(0.96);
	stats2->SetX2NDC(.71);
	stats2->SetLineWidth(2);
	stats2->SetLineColor(kRed);


	h_mc_bar.second->Draw("SAMES E3");
	h_mc_bar.second->SetFillStyle(3325);
	h_mc_bar.first->SetMarkerStyle(22);
	h_mc_bar.first->GetFunction("pol2")->SetLineColor(kRed);
	h_mc_bar.first->GetFunction("pol2")->SetLineWidth(2);
	h_mc_bar.second->SetFillColor(kRed);
	h_mc_bar.first->SetMarkerColor(kRed);
	h_mc_bar.first->SetMarkerSize(1.3);

	
/////////////////////////////////////////////////////////////////////
	
	h_data_end.first->Draw("SAMES");

	c1->Update();
	TPaveStats* stats3 =(TPaveStats*)c1->GetPrimitive("stats");
	stats3->SetName("h3stats");
	stats3->SetY1NDC(.152);
	stats3->SetY2NDC(.252);
	stats3->SetX1NDC(.455);
	stats3->SetX2NDC(.705);
	stats3->SetLineWidth(2);
	stats3->SetLineColor(kMagenta+2);
	//stats3->Draw();


	h_data_end.second->Draw("E3 SAMES");
	h_data_end.second->SetFillStyle(3352);	
        h_data_end.first->SetMarkerStyle(21);
	h_data_end.first->GetFunction("pol2")->SetLineColor(kMagenta+2);
	h_data_end.first->GetFunction("pol2")->SetLineWidth(2);
	h_data_end.second->SetFillColor(kMagenta+2);
	h_data_end.first->SetMarkerColor(kMagenta+2);
	h_data_end.first->SetMarkerSize(1.3);

	
///////////////////////////////////////////////////////////////////////
	
	h_mc_end.first->Draw("SAMES");

	c1->Update();
	TPaveStats* stats4 =(TPaveStats*)c1->GetPrimitive("stats");
	stats4->SetName("h4stats");
	stats4->SetY1NDC(.255);
	stats4->SetY2NDC(.355);
	stats4->SetX1NDC(.455);
	stats4->SetX2NDC(.705);
	stats4->SetLineWidth(2);
	stats4->SetLineColor(kYellow+1);

	h_mc_end.second->Draw("E3 SAMES");
	h_mc_end.second->SetFillStyle(3325);	
	h_mc_end.first->SetMarkerStyle(23);
	h_mc_end.first->GetFunction("pol2")->SetLineColor(kYellow+1);
	h_mc_end.first->GetFunction("pol2")->SetLineWidth(2);
	h_mc_end.second->SetFillColor(kYellow+1);
	h_mc_end.first->SetMarkerColor(kYellow+1);
	h_mc_end.first->SetMarkerSize(1.3);

	//TList* t;
	//t = cnew->GetListOfPrimitives();
	//t->Print();

	auto legend = new TLegend(0.55,0.362,0.85,0.455);
   	legend->SetHeader(legtitle.c_str(),"C"); // option "C" allows to center the header
   	legend->AddEntry(h_data_bar.first,"Data |#eta|<1.2","fp");
   	legend->AddEntry(h_mc_bar.first,"MC |#eta|<1.2","fp");
   	legend->AddEntry(h_data_end.first,"Data |#eta|#geq 1.2","fp");
	legend->AddEntry(h_mc_end.first,"MC |#eta|#geq 1.2","fp");
   	legend->Draw();

	
	c1->SaveAs(("./extrapPlot/effFit_"+to_string(option)+".pdf").c_str());

}
