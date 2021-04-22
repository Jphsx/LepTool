
#include "E_bank.h"
#include "E_bank_fit.h"

class E_bank_manager{

	public:
	//hold all the ebanks
	//muons
	E_bank* id_Jmu_MC;
	E_bank* id_Jmu_Data;
	
	E_bank* id_Zmu_MC;
	E_bank_fit* iso_med_Zmu_MC;
	E_bank_fit* sip_isomed_Zmu_MC;

	E_bank* id_Zmu_Data;
        E_bank_fit* iso_med_Zmu_Data;
        E_bank_fit* sip_isomed_Zmu_Data;

	E_bank_fit* vl_Zmu_Data;
	E_bank_fit* vl_Zmu_MC;
	
	
	//electrons //need to get from andres
	E_bank* id_Zel_MC;
	E_bank* iso_med_Zel_MC;
	E_bank* sip_isomed_Zel_MC;

	E_bank* id_Zel_Data;
        E_bank* iso_med_Zel_Data;
        E_bank* sip_isomed_Zel_Data;

	E_bank* vl_Zel_Data;
	E_bank* vl_Zel_MC;


	E_bank_manager();
	
	double getMCIdEfficiency(double pt, double eta, int pdg, int year);
	double getDataIdEfficiency(double pt, double eta, int pdg, int year);

	double getMCIsoEfficiency(double pt, double eta, int pdg, int year);
	double getDataIsoEfficiency(double pt, double eta, int pdg, int year);
	
	double getMCSipEfficiency(double pt, double eta, int pdg, int year);
	double getDataSipEfficiency(double pt, double eta, int pdg, int year);

	double getMCIdError(double pt, double eta, int pdg, int year);
	double getDataIdError(double pt, double eta, int pdg, int year);

	double getMCIsoError(double pt, double eta, int pdg, int year);
	double getDataIsoError(double pt, double eta, int pdg, int year);
	
	double getMCSipError(double pt, double eta, int pdg, int year);
	double getDataSipError(double pt, double eta, int pdg, int year);

	std::pair<double,double> getMCIdPair(double pt, double eta, int pdg, int year);
	std::pair<double,double> getDataIdPair(double pt, double eta, int pdg, int year);

	std::pair<double,double> getMCIsoPair(double pt, double eta, int pdg, int year);
	std::pair<double,double> getDataIsoPair(double pt, double eta, int pdg, int year);
	
	std::pair<double,double> getMCSipPair(double pt, double eta, int pdg, int year);
	std::pair<double,double> getDataSipPair(double pt, double eta, int pdg, int year);

	double getSF(double effdata, double effmc);
	double getSFerr(double effdata, double effmc, double errdata, double errmc);
	std::pair<double,double> getSFpair(double effdata, double effmc, double errdata, double errmc );
	std::pair<double,double> getSFpair(std::pair<double,double> data, std::pair<double,double> mc );

	std::pair<double,double> getGold(std::pair<double,double> id, std::pair<double,double> iso, std::pair<double,double> sip);
	std::pair<double,double> getSilver(std::pair<double,double> id, std::pair<double,double> iso, std::pair<double,double> sip);
	std::pair<double,double> getBronze(std::pair<double,double> id, std::pair<double,double> iso);
	
	double getGSerr(double v1, double v2, double v3, double e1, double e2, double e3 );
	
	double getBerr(double v1, double v2, double e1, double e2 );

	int verbose =1;
};


E_bank_manager::E_bank_manager(){
	//hardcoded systematics
	//	

	//big collection of paths
	//muon path
	//path for my local computer
	//std::string path = "../";
	//path on unl
	std::string path="/home/t3-ku/janguian/TnP_Muon_Output/";
	std::string pathEl="/home/t3-ku/janguian/TnP_Electron_Output/";

	std::cout<<"Beginning Electron init"<<std::endl;
	std::cout<<"Zel id MC"<<std::endl;
	id_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Tight.root",
				pathEl+"TnPZ_susyID_MC2017_Tight.root",
				pathEl+"TnPZ_susyID_MC2018_tight.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002});
	id_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001});

	std::cout<<"Zel sip MC"<<std::endl;
	iso_med_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Iso.root",
				pathEl+"TnPZ_susyID_MC2017_Iso.root",
				pathEl+"TnPZ_susyID_MC2018_Iso.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005});
        iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001});


	std::cout<<"Zel iso MC"<<std::endl;
	sip_isomed_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Sip3D.root",
				pathEl+"TnPZ_susyID_MC2017_Sip3D.root",
				pathEl+"TnPZ_susyID_MC2018_Sip3D.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002});
        sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001});


	std::cout<<"Zel id Data"<<std::endl;
	id_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Tight.root",
				pathEl+"TnPZ_susyID_data2016_Tight.root",//PLACEHOLDER NO 2017 AVAILABLE
				pathEl+"TnPZ_susyID_data2018_tight.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002});
        id_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001});


	std::cout<<"Zel iso Data"<<std::endl;
	iso_med_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Iso.root",
				pathEl+"TnPZ_susyID_data2016_Iso.root",//PLACEHOLDER NO 2017 AVAILABLE
				pathEl+"TnPZ_susyID_data2018_Iso.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005});
        iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001});

	std::cout<<"Zel sip Data"<<std::endl;
	sip_isomed_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Sip3D.root",
				pathEl+"TnPZ_susyID_data2016_Sip3D.root",//PLACEHOLDER NO 2017 AVAILABLE
				pathEl+"TnPZ_susyID_data2018_Sip3D.root",
				 "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002});
        sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001});

	std::cout<<"Zel VL MC"<<std::endl;
	vl_Zel_MC = new E_bank(2016,
			pathEl+"TnPZ_susyID_MC2016_veryLoose.root",
			pathEl+"TnPZ_susyID_MC2017_veryLoose.root",
			pathEl+"TnPZ_susyID_MC2018_veryLoose.root",
			"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	vl_Zel_MC->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01});
        vl_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003});

	std::cout<<"Zel VL Data"<<std::endl;
	vl_Zel_Data = new E_bank(2016,
			pathEl+"TnPZ_susyID_data2016_veryLoose.root",
			pathEl+"TnPZ_susyID_data2016_veryLoose.root",//PLACEHOLDER NO 2017 AVAILABLE
			pathEl+"TnPZ_susyID_data2018_veryLoose.root",
			"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	vl_Zel_Data->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01});
        vl_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003});

	std::cout<<"Beginning Muon init"<<std::endl;
	/////MUON ID/////////
	std::cout<<"J id data .."<<std::endl;
	id_Jmu_Data = new E_bank(2016,
				path+"TnPJ_MuonID_data2016_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2017_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2018_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Jmu_Data->applySystematics(0.0151);	
	id_Jmu_Data->applySystematics(std::vector<double>{0.001,0.001});

	std::cout<<"Z id data .."<<std::endl;
	id_Zmu_Data = new E_bank(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Zmu_Data->applySystematics(0.0054);
	id_Zmu_Data->applySystematics(std::vector<double>{0.001,0.0003});

	std::cout<<"J id mc .. "<<std::endl;
	id_Jmu_MC = new E_bank(2016,
				path+"TnPJ_MuonID_mc2016_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2017_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2018_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Jmu_MC->applySystematics(0.0151);	
	id_Jmu_MC->applySystematics(std::vector<double>{0.001,0.001});

	std::cout<<"Z id mc .. "<<std::endl;
	id_Zmu_MC = new E_bank(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Zmu_MC->applySystematics(0.0054);
	id_Zmu_MC->applySystematics(std::vector<double>{0.001,0.0003});

	/////MUON ISO///////////

	std::cout<<"Z iso data .. "<<std::endl;
	iso_med_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//iso_med_Zmu_Data->applySystematics(0.0057);	
	iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.002});
	iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.008,0.004});
	iso_med_Zmu_Data->doLowPtFit(3.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//iso_med_Zmu_Data->applySystematicsLow(std::vector<double>{0.008,0.004});
		
	std::cout<<"Z iso MC .. "<<std::endl;		
	iso_med_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//iso_med_Zmu_MC->applySystematics(0.0057);
	iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.002});
	iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.008,0.004});
	iso_med_Zmu_MC->doLowPtFit(3.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//iso_med_Zmu_MC->applySystematicsLow(std::vector<double>{0.008,0.004});	


	//////MUON SIP//////////
	std::cout<<"Z sip data .. "<<std::endl;
	sip_isomed_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//sip_isomed_Zmu_Data->applySystematics(0.0057);
	sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002});
	sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.005,0.003});
	sip_isomed_Zmu_Data->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//sip_isomed_Zmu_Data->applySystematicsLow(std::vector<double>{0.005,0.003});
				
	std::cout<<"Z sip mc .. "<<std::endl;
	sip_isomed_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//sip_isomed_Zmu_MC->applySystematics(0.0057);
	sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002});
	sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.005,0.003});
	sip_isomed_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//sip_isomed_Zmu_MC->applySystematicsLow(std::vector<double>{0.005,0.003});

	/////////MUON VERY LOOSE//////////
	std::cout<<" Z VL data "<<std::endl;
	vl_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//vl_Zmu_Data->applySystematics(0.0054);
	vl_Zmu_Data->applySystematics(std::vector<double>{0.001, 0.0003});	
	vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.001,0.001});
	vl_Zmu_Data->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//vl_Zmu_Data->applySystematicsLow(std::vector<double>{0.001,0.001});

	std::cout<<" Z VL mc "<<std::endl;
	vl_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//vl_Zmu_MC->applySystematics(0.0054);
	vl_Zmu_MC->applySystematics(std::vector<double>{0.001,0.0003});
	vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.001,0.001});
	vl_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//vl_Zmu_MC->applySystematics(std::vector<double>{0.001,0.001});
		
}
///////////eff def
double E_bank_manager::getMCIdEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_MC->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_MC->getEfficiency(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_MC->getEfficiency(pt,eta,year);
		}
		
	}
	return -1;
}
double E_bank_manager::getDataIdEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_Data->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_Data->getEfficiency(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_Data->getEfficiency(pt,eta,year);
		}
		
	}
	return -1;
}
double E_bank_manager::getMCIsoEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_MC->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataIsoEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_Data->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getMCSipEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_MC->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataSipEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_Data->getEfficiency(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getEfficiency(pt,eta,year);
	}
	return -1;
}
/////////////error def
double E_bank_manager::getMCIdError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_MC->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_MC->getError(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_MC->getError(pt,eta,year);
		}
		
	}
	return -1;
}	
double E_bank_manager::getDataIdError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_Data->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_Data->getError(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_Data->getError(pt,eta,year);
		}
		
	}
	return -1;
}
double E_bank_manager::getMCIsoError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_MC->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getError(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataIsoError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_Data->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getError(pt,eta,year);
	}
	return -1;
}	
double E_bank_manager::getMCSipError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_MC->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getError(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataSipError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_Data->getError(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getError(pt,eta,year);
	}
	return -1;
}
////////////////pair def
std::pair<double,double> E_bank_manager::getMCIdPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_MC->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_MC->getPair(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_MC->getPair(pt,eta,year);
		}
		
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getDataIdPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return id_Zel_Data->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		//threshold defaults to 20GeV for jpsi and z
		if(pt < 20.){
			return id_Jmu_Data->getPair(pt,eta,year);
		}
		if(pt > 20.){
			return id_Zmu_Data->getPair(pt,eta,year);
		}
		
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getMCIsoPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_MC->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getDataIsoPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return iso_med_Zel_Data->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getMCSipPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_MC->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getDataSipPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		return sip_isomed_Zel_Data->getPair(pt,eta,year);
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}

double E_bank_manager::getSF(double effdata, double effmc){ return effdata/effmc; }
double E_bank_manager::getSFerr(double effdata, double effmc, double errdata, double errmc){ 
		
	return (effdata/effmc)*std::sqrt( (errdata*errdata)/(effdata*effdata) + (errmc*errmc)/(effmc*effmc) ); 
}
double E_bank_manager::getGSerr(double v1, double v2, double v3, double e1, double e2, double e3 ){
	//gold and silver have the same form
	return (v1*v2*v3)*std::sqrt( (e1*e1)/(v1*v1) + (e2*e2)/(v2*v2) + (e3*e3)/(v3*v3) );	
}
double E_bank_manager::getBerr(double v1, double v2, double e1, double e2 ){
	//bronze has same form as G and S but less parameters
	return (v1*v2)*std::sqrt( (e1*e1)/(v1*v1) + (e2*e2)/(v2*v2)  );
}

std::pair<double,double> E_bank_manager::getSFpair(double effdata, double effmc, double errdata, double errmc ){
	return std::make_pair( getSF(effdata,effmc) , getSFerr( effdata,effmc, errdata,errmc) );
}
std::pair<double,double> E_bank_manager::getSFpair(std::pair<double,double> data, std::pair<double,double> mc ){
	return std::make_pair( getSF(data.first, mc.first), getSFerr(data.first,mc.first, data.second, mc.second));
}
std::pair<double,double> E_bank_manager::getGold(std::pair<double,double> id, std::pair<double,double> iso, std::pair<double,double> sip){
	return std::make_pair( id.first * iso.first * sip.first,  getGSerr(id.first,iso.first,sip.first, id.second,iso.second,sip.second ) );
}
std::pair<double,double> E_bank_manager::getSilver(std::pair<double,double> id, std::pair<double,double> iso, std::pair<double,double> sip){
	return std::make_pair( id.first * iso.first * (1.-sip.first),  getGSerr(id.first,iso.first,(1.-sip.first), id.second,iso.second,sip.second ) );
}
std::pair<double,double> E_bank_manager::getBronze(std::pair<double,double> id, std::pair<double,double> iso){
	return std::make_pair( (1.-(id.first*iso.first)),  getBerr( id.first, iso.first, id.second, iso.second) );
}



int main(){

	E_bank_manager* e1 = new E_bank_manager();
	
	//hax, can do this to do electrons!
	
	E_bank_fit* id_Zmu_Data = (E_bank_fit*) e1->id_Zmu_Data;
	E_bank_fit* id_Zmu_MC = (E_bank_fit*) e1->id_Zmu_MC;
	
	std::cout<<"j16"<<std::endl;
	e1->id_Jmu_Data->printMap(e1->id_Jmu_Data->_map16);
	std::cout<<"zid16"<<std::endl;
	id_Zmu_Data->printMap(id_Zmu_Data->_map16);
	std::cout<<"zid16"<<std::endl;
	e1->id_Zmu_Data->printMap(e1->id_Zmu_Data->_map16);

	std::cout<<"ziso16"<<std::endl;
	e1->iso_med_Zmu_Data->printMap(e1->iso_med_Zmu_Data->_fitmap16);
	e1->iso_med_Zmu_Data->printMap(e1->iso_med_Zmu_Data->_map16);

	e1->id_Zel_MC->printMap(e1->id_Zel_MC->_map16);
	e1->id_Zel_MC->printMap(e1->id_Zel_MC->_map17);
	e1->id_Zel_MC->printMap(e1->id_Zel_MC->_map18);	


}


