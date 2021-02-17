#include "tdrstyle.C"
#include <algorithm>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TRatioPlot.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "TROOT.h"

void convert(TH1F* h, TGraphAsymmErrors* g){
        auto nPoints = g->GetN();
        for(int i=0; i < nPoints; ++i){
                double x,y;
                g->GetPoint(i,x,y);
                h->Fill(x,y);
                int binx = h->FindBin(x);
                h->SetBinError(binx, g->GetErrorY(i) );
		if(g->GetErrorY(i) > 0.05){
			h->SetBinError(binx, g->GetErrorYlow(i) );
		}
		//std::cout<<g->GetErrorY(i)<<" "<<std::endl;;
        }
}

void formatratio(TRatioPlot* h){
h->GetLowerRefYaxis()->SetTitle("Data/MC");
h->GetLowerRefYaxis()->SetTitleSize(0.03);
h->GetLowerRefYaxis()->SetTitleOffset(1.48);
h->GetLowerRefYaxis()->SetLabelSize(0.03);

h->GetLowerRefXaxis()->SetLabelSize(0.04);
h->GetLowerRefXaxis()->SetLabelOffset(0.004);
h->GetLowerRefXaxis()->SetTitleSize(0.04);

h->GetUpperRefYaxis()->SetTitleSize(0.04);
h->GetUpperRefYaxis()->SetLabelSize(0.04);
h->GetUpperRefYaxis()->SetLabelOffset(0.001);


}

void makeratio(TCanvas* data, TCanvas* mc, std::string axis_label, int bins, double* edges, std::string legend_header, std::string uniqueID, std::string savename){


		

	TGraphAsymmErrors* graph_data = (TGraphAsymmErrors*) data->GetPrimitive("hxy_fit_eff");
	TGraphAsymmErrors* graph_mc = (TGraphAsymmErrors*) mc->GetPrimitive("hxy_fit_eff");

	TH1F* h_data = new TH1F(("h_data"+uniqueID).c_str(),("Data;"+axis_label).c_str(), bins, edges);
	TH1F* h_mc = new TH1F(("h_mc"+uniqueID).c_str(),("MC;"+axis_label).c_str(), bins, edges);

	convert( h_data, graph_data);
	convert( h_mc, graph_mc);

	//h_mc->SetFillColor(46);
	h_mc->SetMarkerColor(kRed);
	h_mc->SetLineColor(kRed);

	h_data->SetMarkerStyle(8);



	TCanvas* crat = new TCanvas(("canvas"+uniqueID).c_str(), ("canvas"+uniqueID).c_str());
	TRatioPlot* rat = new TRatioPlot(h_data,h_mc);
	rat->SetH2DrawOpt("E");
	rat->SetH1DrawOpt("E");
	rat->Draw();
	crat->SetLeftMargin(33);
	formatratio(rat);
	rat->GetUpperRefYaxis()->SetRangeUser(.8,1.);
	rat->GetLowerRefYaxis()->SetRangeUser(.95,1.05);
	rat->GetUpperRefYaxis()->SetLabelSize(0.03);
	crat->Update();
	TPad *pad = rat->GetUpperPad();
	auto legend = pad->BuildLegend();
	legend->SetHeader(legend_header.c_str(),"C");
	legend->SetX1NDC(0.01);
	legend->SetX2NDC(3);
	crat->Update();
	crat->SaveAs(("./plots/canvas"+uniqueID+savename+".pdf").c_str());
}

void efficiency_plotter(int proc, int year){
setTDRStyle();////////

	////////////// ratio plotting for ID
	int uniqueId =0; 
	double* edges;
	int nbins;
	std::string proc_label;
	std::string id_data;
	std::string id_mc;

	std::string isomed_data;
	std::string isomed_mc;

	std::string sipisomed_data;
	std::string sipisomed_mc;

	std::string trigger="";
	std::string canvasPath_barrel;
	std::string canvasPath_endcap;
	std::string savename;
	if(proc==1){//Z
	    proc_label = "Z";
	    nbins = 5;
	    Double_t ztemp[5+1] = {10, 20, 30, 40, 60, 100};
	    edges = new double[nbins+1];
	    std::copy(ztemp, ztemp+nbins+1, &edges[0]);
	    if(year == 2016){
		        proc_label = "2016 "+proc_label;
			id_mc = "TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__Medium_A_pt1.root";
			id_data = "TnP_MuonID_data2016_IsoTkMu22_veryLoose__Medium_A_pt1.root";
			
			isomed_mc = "TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root";
			isomed_data = "TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root";

			sipisomed_mc = "TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root";
			sipisomed_data = "TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root";

			canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";//_tag_IsoTkMu22_pass";
			canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";//_tag_IsoTkMu22_pass";
			trigger="_tag_IsoTkMu22_pass";
			savename = "Z2016";
	    }
	    if(year == 2017){
			proc_label = "2017 "+proc_label;
			id_mc = "TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root";
			id_data = "TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__Medium_A_pt1.root";
		
			isomed_mc = "TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root";
			isomed_data = "TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root";
		
			sipisomed_mc = "TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root";
			sipisomed_data = "TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root";
		
			canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";
			canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";
			trigger="_tag_IsoMu24_eta2p1_pass";
			savename = "Z2017";
	    }
	    if(year == 2018){
		    	proc_label = "2018 "+proc_label;
			id_mc = "TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root";
			id_data = "TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__Medium_A_pt1.root";

			isomed_mc = "TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root";
			isomed_data = "TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root";

			sipisomed_mc = "TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root";
			sipisomed_data = "TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root";

			canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";
			canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";
			trigger ="_tag_IsoMu24_eta2p1_pass";
			savename = "Z2018";
	    }
	}    
	if(proc==2){//J
	    proc_label = "J/#psi";
            nbins = 7;
            Double_t ztemp[7+1] = {3.0, 4.0,  5.0, 6.0, 7.0, 9.0, 14.0, 20.0};
            edges = new double[nbins+1];
            std::copy(ztemp, ztemp+nbins+1, &edges[0]);
	    if(year == 2016){
	        proc_label = "2016 "+proc_label;
		id_mc = "TnPJ_MuonID_mc2016_weight_Mu7p5Tk2_veryLoose__Medium_E.root";
		id_data = "TnPJ_MuonID_data2016_Mu7p5Tk2_veryLoose__Medium_E.root";
		canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";
                canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";
                trigger = "_tag_Mu7p5_Track2_Jpsi_MU_pass";
		savename = "J2016";
	    }
	    if(year == 2017){
		proc_label = "2017 "+proc_label;
                id_mc = "TnPJ_MuonID_mc2017_weight_Mu7p5Tk2_veryLoose__Medium_E.root";
                id_data = "TnPJ_MuonID_data2017_Mu7p5Tk2_veryLoose__Medium_E.root";
		canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";
                canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";
                trigger = "_tag_Mu7p5_Track2_Jpsi_MU_pass";
		savename = "J2017";
            }
	    if(year == 2018){
		proc_label = "2018 "+proc_label;
                id_mc = "TnPJ_MuonID_mc2018_weight_Mu7p5Tk2_veryLoose__Medium_E.root";
                id_data = "TnPJ_MuonID_data2018_Mu7p5Tk2_veryLoose__Medium_E.root";
		canvasPath_barrel = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin0_&";
                canvasPath_endcap = "tpTree/Medium_pt_eta/fit_eff_plots/pt_PLOT_abseta_bin1_&";
                trigger = "_tag_Mu7p5_Track2_Jpsi_MU_pass";
		savename = "J2018";
 	    }
	}
	
	//for(int i=0; i<6;i++){
	//	std::cout<<edges[i]<<" ";
	//}	

		
////////////////////id plots	
	TFile* fid_data = TFile::Open(id_data.c_str());
	TFile* fid_mc = TFile::Open(id_mc.c_str());

	//do barrel
	TCanvas* canv_data = (TCanvas*) fid_data->Get((canvasPath_barrel+trigger).c_str());
	TCanvas* canv_mc = (TCanvas*) fid_mc->Get((canvasPath_barrel+trigger).c_str());
	

	makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{ID}", nbins, edges, proc_label+" Muons |#eta|<1.2", std::to_string(uniqueId++), savename+"b");
	//do endcap

	canv_data = (TCanvas*) fid_data->Get((canvasPath_endcap+trigger).c_str());
	canv_mc = (TCanvas*) fid_mc->Get((canvasPath_endcap+trigger).c_str());
	makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{ID}", nbins, edges, proc_label+" Muons |#eta|#geq1.2",std::to_string(uniqueId++), savename+"e");

	if(proc == 2) return;//we only make medium id for jpsi

///////////////////isomed plots
	TFile* fisomed_data = TFile::Open(isomed_data.c_str());
	TFile* fisomed_mc = TFile::Open(isomed_mc.c_str());
	canv_data  = (TCanvas*) fisomed_data->Get((canvasPath_barrel+"_Medium_pass_&"+trigger).c_str());
	canv_mc = (TCanvas*) fisomed_mc->Get((canvasPath_barrel+"_Medium_pass_&"+trigger).c_str());
	
	makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{Isolated|ID}", nbins, edges, proc_label+" Muons |#eta|<1.2",std::to_string(uniqueId++), savename+"b");

	canv_data  = (TCanvas*) fisomed_data->Get((canvasPath_endcap+"_Medium_pass_&"+trigger).c_str());
        canv_mc = (TCanvas*) fisomed_mc->Get((canvasPath_endcap+"_Medium_pass_&"+trigger).c_str());

	makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{Isolated|ID}", nbins, edges, proc_label+" Muons |#eta|#geq1.2",std::to_string(uniqueId++),savename+"e");

//////////////////sipisomed plots
	TFile* fsipisomed_data = TFile::Open(sipisomed_data.c_str());
	TFile* fsipisomed_mc = TFile::Open(sipisomed_mc.c_str());
	canv_data = (TCanvas*) fsipisomed_data->Get((canvasPath_barrel+"_Medium_pass_&"+trigger).c_str());
	canv_mc = (TCanvas*) fsipisomed_mc->Get((canvasPath_barrel+"_Medium_pass_&"+trigger).c_str());

	makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{Prompt|ID #cap Isolated}", nbins, edges, proc_label+" Muons |#eta|<1.2",std::to_string(uniqueId++),savename+"b");

	canv_data = (TCanvas*) fsipisomed_data->Get((canvasPath_endcap+"_Medium_pass_&"+trigger).c_str());
        canv_mc = (TCanvas*) fsipisomed_mc->Get((canvasPath_endcap+"_Medium_pass_&"+trigger).c_str());

        makeratio(canv_data, canv_mc, "p_{T} [GeV];Efficiency #varepsilon_{Prompt|ID #cap Isolated}", nbins, edges, proc_label+" Muons |#eta|#geq1.2",std::to_string(uniqueId++),savename+"e");

}
