
#include "E_bank.h"
#include "E_bank_fit.h"
#include "TText.h"
#include "tdrstyle.C"
#include "TLegend.h"
#include "TStyle.h"
#include "THStack.h"

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
	
	
	E_bank* id_Zel_MC;
	E_bank* iso_med_Zel_MC;
	E_bank* sip_isomed_Zel_MC;

	E_bank* id_Zel_Data;
        E_bank* iso_med_Zel_Data;
        E_bank* sip_isomed_Zel_Data;

	E_bank* vl_Zel_Data;
	E_bank* vl_Zel_MC;

	//fullsim fast sim 
	E_bank* id_SimF_el;
	E_bank* iso_med_SimF_el;
	E_bank* sip_isomed_SimF_el;
	E_bank* vl_SimF_el;

	E_bank* id_SimF_mu;
        E_bank* iso_med_SimF_mu;
        E_bank* sip_isomed_SimF_mu;
        E_bank* vl_SimF_mu;


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

	double getIdSimSF(double pt, double eta, int pdg, int year);
	double getIdSimSFerr(double pt, double eta, int pdg, int year);
	std::pair<double, double> getIdSimSFpair( double pt, double eta, int pdg, int year);

	double getIsoSimSF(double pt, double eta, int pdg, int year);
        double getIsoSimSFerr(double pt, double eta, int pdg, int year);
        std::pair<double, double> getIsoSimSFpair( double pt, double eta, int pdg, int year);

	double getSipSimSF(double pt, double eta, int pdg, int year);
        double getSipSimSFerr(double pt, double eta, int pdg, int year);
        std::pair<double, double> getSipSimSFpair( double pt, double eta, int pdg, int year);

	double getVlSimSF(double pt, double eta, int pdg, int year);
        double getVlSimSFerr(double pt, double eta, int pdg, int year);
        std::pair<double, double> getVlSimSFpair( double pt, double eta, int pdg, int year);




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
	std::string path="./TnP_Muon_Output/";
	std::string pathEl="./TnP_Electron_Output/";
	std::string pathSimF="./FullSim-FastSim_SF/";


	std::cout<<"Beginning Electron init"<<std::endl;
	std::cout<<"Zel id MC"<<std::endl;
	id_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Tight.root",
				pathEl+"TnPZ_susyID_MC2017_Tight.root",
				pathEl+"TnPZ_susyID_MC2018_tight.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002}, 2018);
	id_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001}, 2018);
	id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.001, 0.005}, 2017);
	id_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
	id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.003, 0.003}, 2016);
	id_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.002}, 2016);

	std::cout<<"Zel iso MC"<<std::endl;
	iso_med_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Iso.root",
				pathEl+"TnPZ_susyID_MC2017_Iso.root",
				pathEl+"TnPZ_susyID_MC2018_Iso.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002}, 2018);
    iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001}, 2018);
    iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.003, 0.003}, 2017);
    iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
    iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.004, 0.003}, 2016);
    iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2016);


	std::cout<<"Zel sip MC"<<std::endl;
	sip_isomed_Zel_MC = new E_bank(2016,
				pathEl+"TnPZ_susyID_MC2016_Sip3D.root",
				pathEl+"TnPZ_susyID_MC2017_Sip3D.root",
				pathEl+"TnPZ_susyID_MC2018_Sip3D.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005}, 2018);
    sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2018);
	sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.006, 0.004, 0.007}, 2017);
    sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.0006}, 2017);
    sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.001, 0.007, 0.011}, 2016);
    sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.005, 0.012, 0.003}, 2016);

	std::cout<<"Zel id Data"<<std::endl;
	id_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Tight.root",
				pathEl+"TnPZ_susyID_data2017_Tight.root",
				pathEl+"TnPZ_susyID_data2018_tight.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002}, 2018);
	id_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001}, 2018);
	id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.001, 0.005}, 2017);
	id_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
	id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.003, 0.003}, 2016);
	id_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.002}, 2016);;


	std::cout<<"Zel iso Data"<<std::endl;
	iso_med_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Iso.root",
				pathEl+"TnPZ_susyID_data2017_Iso.root",
				pathEl+"TnPZ_susyID_data2018_Iso.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002}, 2018);
    iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001}, 2018);
    iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.003, 0.003}, 2017);
    iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
    iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.004, 0.003}, 2016);
    iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2016);

	std::cout<<"Zel sip Data"<<std::endl;
	sip_isomed_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Sip3D.root",
				pathEl+"TnPZ_susyID_data2017_Sip3D.root",
				pathEl+"TnPZ_susyID_data2018_Sip3D.root",
				 "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005}, 2018);
    sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2018);
	sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.006, 0.004, 0.007}, 2017);
    sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.0006}, 2017);
    sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.001, 0.007, 0.011}, 2016);
    sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.005, 0.012, 0.003}, 2016);

	std::cout<<"Zel VL MC"<<std::endl;
	vl_Zel_MC = new E_bank(2016,
			pathEl+"TnPZ_susyID_MC2016_veryLoose.root",
			pathEl+"TnPZ_susyID_MC2017_veryLoose.root",
			pathEl+"TnPZ_susyID_MC2018_veryLoose.root",
			"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	vl_Zel_MC->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01}, 2018);
    vl_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003}, 2018);
    vl_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.007, 0.03}, 2017);
    vl_Zel_MC->applySystematic_ptRange(20., 999., {0.003, 0.0001, 0.0007}, 2017);
    vl_Zel_MC->applySystematic_ptRange(0., 20., {0.02, 0.04, 0.02}, 2016);
    vl_Zel_MC->applySystematic_ptRange(20., 999., {0.003, 0.005, 0.003}, 2016);

	std::cout<<"Zel VL Data"<<std::endl;
	vl_Zel_Data = new E_bank(2016,
			pathEl+"TnPZ_susyID_data2016_veryLoose.root",
			pathEl+"TnPZ_susyID_data2017_veryLoose.root",
			pathEl+"TnPZ_susyID_data2018_veryLoose.root",
			"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
	vl_Zel_Data->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01}, 2018);
    vl_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003}, 2018);
    vl_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.007, 0.03}, 2017);
    vl_Zel_Data->applySystematic_ptRange(20., 999., {0.003, 0.0001, 0.0007}, 2017);
    vl_Zel_Data->applySystematic_ptRange(0., 20., {0.02, 0.04, 0.02}, 2016);
    vl_Zel_Data->applySystematic_ptRange(20., 999., {0.003, 0.005, 0.003}, 2016);

	std::cout<<"Beginning Muon init"<<std::endl;
	/////MUON ID/////////
	std::cout<<"J id data .."<<std::endl;
	id_Jmu_Data = new E_bank(2016,
				path+"TnPJ_MuonID_data2016_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2017_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2018_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Jmu_Data->applySystematics(0.0151);	
	id_Jmu_Data->applySystematics(std::vector<double>{0.001,0.001}, 2018);
	id_Jmu_Data->applySystematics(std::vector<double>{0.001,0.001}, 2017);
	id_Jmu_Data->applySystematics(std::vector<double>{0.002,0.001}, 2016);

	std::cout<<"Z id data .."<<std::endl;
	id_Zmu_Data = new E_bank(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Zmu_Data->applySystematics(0.0054);
	id_Zmu_Data->applySystematics(std::vector<double>{0.0007,0.0003}, 2018);
	id_Zmu_Data->applySystematics(std::vector<double>{0.001,0.0003}, 2017);
    id_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002}, 2016);
    
	std::cout<<"J id mc .. "<<std::endl;
	id_Jmu_MC = new E_bank(2016,
				path+"TnPJ_MuonID_mc2016_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2017_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2018_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Jmu_MC->applySystematics(0.0151);	
	id_Jmu_MC->applySystematics(std::vector<double>{0.001,0.001}, 2018);
	id_Jmu_MC->applySystematics(std::vector<double>{0.001,0.001}, 2017);
	id_Jmu_MC->applySystematics(std::vector<double>{0.002,0.001}, 2016);

	std::cout<<"Z id mc .. "<<std::endl;
	id_Zmu_MC = new E_bank(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//id_Zmu_MC->applySystematics(0.0054);
	id_Zmu_MC->applySystematics(std::vector<double>{0.0007,0.0003}, 2018);
	id_Zmu_MC->applySystematics(std::vector<double>{0.001,0.0003}, 2017);
    id_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002}, 2016);

	/////MUON ISO///////////

	std::cout<<"Z iso data .. "<<std::endl;
	iso_med_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//iso_med_Zmu_Data->applySystematics(0.0057);	
	iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.004},2018);
	iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.002},2017);
	iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.004},2016);
	iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.007,0.002},2018);
	iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.007,0.004},2017);
	iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.005,0.0009},2016);
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
	iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.004},2018);
	iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.002},2017);
	iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.004},2016);
	iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.007,0.002},2018);
	iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.007,0.004},2017);
	iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.005,0.0009},2016);
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
	sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.0004,0.001},2018);
	sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002},2017);
	sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002},2016);	
	sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.003,0.002},2018);
	sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.005,0.003},2017);
	sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.002,0.002},2016);
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
	sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.0004,0.001},2018);
	sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002},2017);
	sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002},2016);	
	sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.003,0.002},2018);
	sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.005,0.003},2017);
	sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.002,0.002},2016);
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
	vl_Zmu_Data->applySystematics(std::vector<double>{0.0002, 0.0001}, 2018);
	vl_Zmu_Data->applySystematics(std::vector<double>{0.001, 0.001}, 2017);
	vl_Zmu_Data->applySystematics(std::vector<double>{0.0004, 0.002}, 2016);	
	vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.0008,0.003}, 2018);
	vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.001,0.0003}, 2017);
	vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.002,0.003}, 2016);
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
	vl_Zmu_MC->applySystematics(std::vector<double>{0.0002, 0.0001}, 2018);
	vl_Zmu_MC->applySystematics(std::vector<double>{0.001, 0.001}, 2017);
	vl_Zmu_MC->applySystematics(std::vector<double>{0.0004, 0.002}, 2016);	
	vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.0008,0.003}, 2018);
	vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.001,0.0003}, 2017);
	vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.002,0.003}, 2016);
	vl_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	//vl_Zmu_MC->applySystematics(std::vector<double>{0.001,0.001});

	/////fullsim fast sim SF
	std::cout<<" el simF"<<std::endl;

	id_SimF_el =  new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_el_sf_ID.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_el_sf_ID.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_ID.root",
                                "");


	iso_med_SimF_el = new E_bank(2016,
                		pathSimF+"TTJets_FastVsFull_Summer16_el_sf_iso.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_el_sf_iso.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_iso.root",
                                "");

	sip_isomed_SimF_el = new E_bank(2016,
				pathSimF+"TTJets_FastVsFull_Summer16_el_sf_sip3D.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_el_sf_sip3D.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_sip3D.root",
                                ""); 
        vl_SimF_el = new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_el_sf_vl.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_el_sf_vl.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_vl.root",
                                "");

	std::cout<<" mu simF"<<std::endl;
        id_SimF_mu =  new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_ID.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_ID.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_ID.root",
                                "");

        iso_med_SimF_mu = new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_iso.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_iso.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_iso.root",
                                "");

        sip_isomed_SimF_mu = new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_sip3D.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_sip3D.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_sip3D.root",
                                "");
        vl_SimF_mu = new E_bank(2016,
                                pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_vl.root",
                                pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_vl.root",
                                pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_vl.root",
                                "");



	
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

//sim sf definitions
double E_bank_manager::getIdSimSF(double pt, double eta, int pdg, int year){
	if(abs(pdg) == 11){//electron
                return id_SimF_el->getEfficiency(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return id_SimF_mu->getEfficiency(pt,eta,year);
        }
        return -1;

}
double  E_bank_manager::getIdSimSFerr(double pt, double eta, int pdg, int year){
	if(abs(pdg) == 11){//electron
                return id_SimF_el->getError(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return id_SimF_mu->getError(pt,eta,year);
        }
        return -1;
}
std::pair<double, double>  E_bank_manager::getIdSimSFpair( double pt, double eta, int pdg, int year){
	if(abs(pdg) == 11){//electron
                return id_SimF_el->getPair(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return id_SimF_mu->getPair(pt,eta,year);
        }
        return  std::make_pair(-1.,-1.);
}

double E_bank_manager::getIsoSimSF(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return iso_med_SimF_el->getEfficiency(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return iso_med_SimF_mu->getEfficiency(pt,eta,year);
        }
        return -1;

}
double  E_bank_manager::getIsoSimSFerr(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return iso_med_SimF_el->getError(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return iso_med_SimF_mu->getError(pt,eta,year);
        }
        return -1;
}
std::pair<double, double>  E_bank_manager::getIsoSimSFpair( double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return iso_med_SimF_el->getPair(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return iso_med_SimF_mu->getPair(pt,eta,year);
        }
        return  std::make_pair(-1.,-1.);
}

double E_bank_manager::getSipSimSF(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return sip_isomed_SimF_el->getEfficiency(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return sip_isomed_SimF_mu->getEfficiency(pt,eta,year);
        }
        return -1;

}
double  E_bank_manager::getSipSimSFerr(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return sip_isomed_SimF_el->getError(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return sip_isomed_SimF_mu->getError(pt,eta,year);
        }
        return -1;
}
std::pair<double, double>  E_bank_manager::getSipSimSFpair( double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return sip_isomed_SimF_el->getPair(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return sip_isomed_SimF_mu->getPair(pt,eta,year);
        }
        return  std::make_pair(-1.,-1.);
}

double E_bank_manager::getVlSimSF(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return vl_SimF_el->getEfficiency(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return vl_SimF_mu->getEfficiency(pt,eta,year);
        }
        return -1;

}
double  E_bank_manager::getVlSimSFerr(double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return vl_SimF_el->getError(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return vl_SimF_mu->getError(pt,eta,year);
        }
        return -1;
}
std::pair<double, double>  E_bank_manager::getVlSimSFpair( double pt, double eta, int pdg, int year){
        if(abs(pdg) == 11){//electron
                return vl_SimF_el->getPair(pt,eta,year);
        }
        if(abs(pdg) == 13){//muon
                return vl_SimF_mu->getPair(pt,eta,year);
        }
        return  std::make_pair(-1.,-1.);
}



//void make2d( E_bank_manager* bankman, E_bank* dat, E_bank* mc, int year, int rank, std::string name){//involves no splitting
void make2d ( E_bank_manager* bankman,int year, int rank, std::string name, //E_bank_fit* bankid_data, E_bank_fit* bankid_mc, 
E_bank* bankid_data, E_bank* bankid_mc, 
E_bank* bankiso_data=0,E_bank* bankiso_mc=0, 
E_bank* banksip_data=0,E_bank* banksip_mc=0,
E_bank* bankJ_data=0, E_bank* bankJ_mc=0 ){


	//RANK 0 VL RANK 1 GOLD RANK 2 SILVER RANK 3 BRONZE

	std::map<std::pair<double,double>, std::pair<double,double> > dataid;
	std::map<std::pair<double,double>, std::pair<double,double> > dataidlow;
	

	std::map<std::pair<double,double>, std::pair<double,double> > dataiso;
	std::map<std::pair<double,double>, std::pair<double,double> > dataisolow;
/*
	if(rank==0){
		if(year==2016){
			dataid = bankid_data->_map16;
			dataidlow = bankid_data->_fitmap16;
		}
		if(year==2017){
			dataid = bankid_data->_map17;
			dataidlow = bankid_data->_fitmap17;
		}
		if(year==2018){
			dataid = bankid_data->_map18;
			dataidlow = bankid_data->_fitmap18;

		}
	}
*/
//	if(rank>0){
		if(year==2016){
			dataiso = bankiso_data->_map16;
			dataisolow = bankiso_data->_map16;

		}
		if(year==2017){

			dataiso = bankiso_data->_map17;
			dataisolow = bankiso_data->_map17;
			
		}
		if(year==2018){
			dataiso = bankiso_data->_map18;
			dataisolow = bankiso_data->_map18;	
		}
//	}

	TH2D* hmc;
	TH2D* hdata;
	TH2D* hsf;
	//count bins, get keys as bin edges
//pick from iso because i can control these bins
	std::set<double> ptset;
	std::set<double> etaset;

//	if(rank>0){//loop on iso
/*		for(auto& itr : dataisolow){
			ptset.insert(itr.first.first);
			etaset.insert(itr.first.second);
		}
*/ //remove datisolow for electrons it doesnt have fit at low pt edges
		for(auto& itr : dataiso){
			//f(itr.first.first >=20){
				ptset.insert(itr.first.first);
			//}
			etaset.insert(itr.first.second);
		}

		ptset.insert( bankiso_data->_ptUpperEdge );
		etaset.insert( bankiso_data->_etaUpperEdge);
		
		ptset.erase(0);//electrons start at 5
std::cout<<"edges \n";
for( auto itr : ptset){
	std::cout<<itr<<" ";
}
std::cout<<std::endl;
for(auto itr : etaset){
	std::cout<<itr<<" ";
}
std::cout<<std::endl;
//	}	
/*	if(rank==0){//loop on id
		for(auto& itr : dataidlow){
			if(itr.first.first <20){
				ptset.insert(itr.first.first);
			}
			etaset.insert(itr.first.second);
		}
		for(auto& itr : dataid){
			ptset.insert(itr.first.first);
			etaset.insert(itr.first.second);
		}	
		ptset.insert( bankid_data->_ptUpperEdge );
		etaset.insert( bankid_data->_etaUpperEdge);
	}
*/
	int nbinseta = etaset.size()-1;
	int nbinspt = ptset.size()-1;
	int nedgeeta = etaset.size();
	int nedgept = ptset.size();
	//create bin edges from set
	double* ptedge = new double[nedgept];
	double* etaedge = new double[nedgeeta];

	
	int i=0;
	for(auto& itr: etaset ){
		etaedge[i] = itr;	
		//std::cout<<itr<<" ";
		i++; 
	}
	std::cout<<std::endl;
	i=0;
	for(auto& itr: ptset){
		ptedge[i] = itr;
		//std::cout<<itr<<" ";
		i++;
	}
	/* debugging prints
	std::cout<<std::endl;
	//quick test
	std::cout<<nbinspt<<" "<<nbinseta<<std::endl;
	for(int i=0; i<nedgept; i++){
		std::cout<<ptedge[i]<<" ";
	}	
	std::cout<<std::endl;
	for(int i=0; i<nedgeeta; i++){
		std::cout<<etaedge[i]<<" ";
	}
	std::cout<<std::endl;
	
	*/

	gStyle->SetPalette(kAquamarine);
	hmc = new TH2D(("h"+std::to_string(year)+std::to_string(rank)+"m").c_str(),";#eta;p_{T} [GeV];Data/MC",nbinseta ,etaedge, nbinspt, ptedge);
	hdata = new TH2D(("h"+std::to_string(year)+std::to_string(rank)+"d").c_str(),";#eta;p_{T} [GeV];Data/MC",nbinseta ,etaedge, nbinspt, ptedge);
	hsf = new TH2D(("h"+std::to_string(year)+std::to_string(rank)+"s").c_str(),";#eta;p_{T} [GeV];Data/MC",nbinseta ,etaedge, nbinspt, ptedge);
	
	//loop bin by bin get values+errors from map and set bin content
	int nbinsX,nbinsY;
	nbinsX= hmc->GetNbinsX();
	nbinsY= hmc->GetNbinsY();
	double binContent1,binContent2,binContent3, binError1,binError2,binError3;
	double effid_data,effiso_data,effsip_data, effid_mc,effiso_mc,effsip_mc;
	double errid_data,erriso_data,errsip_data, errid_mc,erriso_mc,errsip_mc;
	for(int i=1; i<=nbinsX; i++){
		for(int j=1; j<=nbinsY; j++){
		
			//in the case of VL
			if(rank == 0 ){
				//hax treat vl like iso because extrap, this gets us past needing to slicej), 
				effid_data = bankiso_data->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
				errid_data = bankiso_data->getError( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
				effid_mc = bankiso_mc->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
				errid_mc = bankiso_mc->getError(hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					
			
			
				
			//set VL bincontent
				binContent2 = effid_data;
				binError2 = errid_data;
				binContent1 = effid_mc;
				binError1 = errid_mc;
				binContent3 = binContent2/binContent1;
				binError3 = bankman->getSFerr(binContent2,binContent1, binError2, binError1);
	
			}//end rank0 check
			if(rank > 0 ){
					
				
//				if( hdata->GetYaxis()->GetBinLowEdge(j) >= 20. ){
					//std::cout<<"high pt eval id "<< hdata->GetYaxis()->GetBinLowEdge(j)<<" "<<hdata->GetXaxis()->GetBinLowEdge(i)<<std::endl;
					effid_data = bankid_data->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errid_data = bankid_data->getError( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					effid_mc = bankid_mc->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errid_mc = bankid_mc->getError(hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
//				}
/*				if(hdata->GetYaxis()->GetBinLowEdge(j)<20.){//pick jpsi in id
					//std::cout<<"low pt eval id "<<hdata->GetYaxis()->GetBinLowEdge(j)<<" "<<hdata->GetXaxis()->GetBinLowEdge(i)<<std::endl;
					effid_data = bankJ_data->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errid_data = bankJ_data->getError( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					effid_mc = bankJ_mc->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errid_mc = bankJ_mc->getError(hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
				}
*///remove jpsi selection
					//everything else sets normally
					//std::cout<<"eval iso "<<hdata->GetYaxis()->GetBinLowEdge(j)<<" "<<hdata->GetXaxis()->GetBinLowEdge(i)<<std::endl;
					effiso_data = bankiso_data->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					erriso_data = bankiso_data->getError( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					effiso_mc = bankiso_mc->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					erriso_mc = bankiso_mc->getError(hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);

					//std::cout<<"eval sip "<<hdata->GetYaxis()->GetBinLowEdge(j)<<" "<<hdata->GetXaxis()->GetBinLowEdge(i)<<std::endl;
					effsip_data = banksip_data->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errsip_data = banksip_data->getError( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					effsip_mc = banksip_mc->getEfficiency( hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					errsip_mc = banksip_mc->getError(hdata->GetYaxis()->GetBinLowEdge(j), hdata->GetXaxis()->GetBinLowEdge(i), year);
					
				
			// gsb content
			if(rank == 1){//gold
				binContent2 = effid_data*effiso_data*effsip_data;
				binError2 = bankman->getGSerr(effid_data, effiso_data, effsip_data, errid_data, erriso_data, errsip_data );

				binContent1 = effid_mc*effiso_mc*effsip_mc;
				binError1 = bankman->getGSerr(effid_mc, effiso_mc, effsip_mc, errid_mc, erriso_mc, errsip_mc );

				binContent3 = binContent2/binContent1;
				binError3 = bankman->getSFerr(binContent2,binContent1, binError2, binError1);
			}
			if(rank == 2){//silver
				binContent2 = effid_data*effiso_data*(1.-effsip_data);
				binError2 = bankman->getGSerr(effid_data, effiso_data, (1.-effsip_data), errid_data, erriso_data, errsip_data );

				binContent1 = effid_mc*effiso_mc*(1.-effsip_mc);
				binError1 = bankman->getGSerr(effid_mc, effiso_mc, (1.-effsip_mc), errid_mc, erriso_mc, errsip_mc );

				binContent3 = binContent2/binContent1;
				binError3 = bankman->getSFerr(binContent2,binContent1, binError2, binError1);
			}
			if(rank == 3){//bronze
				binContent2 = 1.-(effid_data*effiso_data);
				binError2 = bankman->getBerr(effid_data, effiso_data, errid_data, erriso_data );

				binContent1 = 1.-(effid_mc*effiso_mc);
				binError1 = bankman->getBerr(effid_mc, effiso_mc, errid_mc, erriso_mc );

				binContent3 = binContent2/binContent1;
				binError3 = bankman->getSFerr(binContent2,binContent1, binError2, binError1);
			}
			
			}//end rank>0 check
			///checks
			/*debugging prints
			double X = hdata->GetYaxis()->GetBinLowEdge(j);
			double Y = hdata->GetXaxis()->GetBinLowEdge(i);
			std::cout<<"bin "<<i<<" "<<j<<" X "<<X<<" Y "<<Y<<std::endl;
			std::cout<<"data "<<effid_data<<" "<<effiso_data<<" "<<effsip_data<<std::endl;
			//std::cout<<"err  "<<errid_data<<" "<<erriso_data<<" "<<errsip_data<<std::endl;
			std::cout<<"mc "<<effid_mc<<" "<<effiso_mc<<" "<<effsip_mc<<std::endl;
			//std::cout<<"err  "<<errid_mc<<" "<<erriso_mc<<" "<<errsip_mc<<std::endl;
			std::cout<<binContent1<<" "<<binContent2<<" "<<binContent3<<std::endl;
			std::cout<<std::endl;
			*/
			//fill hists
			hmc->SetBinContent(i,j,binContent1);
			hmc->SetBinError(i,j,binError1);
			hdata->SetBinContent(i,j,binContent2);
			hdata->SetBinError(i,j,binError2);
			hsf->SetBinContent(i,j,binContent3);
			hsf->SetBinError(i,j,binError3);
		}//end j loop
	}//end i loop
	//from copied histograms create new super hist with 



	TCanvas* c = new TCanvas();
	c->SetRightMargin(0.23);
	
	hmc->SetMarkerColor(kRed);
	hsf->Draw("COLZ");
	for (int i=1; i<=nbinsX; i++) {
      		for (int j=1; j<=nbinsY; j++) {
        	 auto t = new TText(hmc->GetXaxis()->GetBinCenter(i),  //was +0.25
                            hmc->GetYaxis()->GetBinCenter(j)+1.95,
                            Form(" %4.3f +/- %4.3f",hmc->GetBinContent(i,j), hmc->GetBinError(i,j)));
		 t->SetTextSize(0.02);
		 t->SetTextColor(kRed+2);
        	 t->SetTextAlign(22);
        	 t->Draw();
      		}
  	}
	for (int i=1; i<=nbinsX; i++) {
      		for (int j=1; j<=nbinsY; j++) {
        	 auto t = new TText(hdata->GetXaxis()->GetBinCenter(i),  //was -0.25
                            hdata->GetYaxis()->GetBinCenter(j),  //was +0.5
                            Form(" %4.3f +/- %4.3f",hdata->GetBinContent(i,j), hdata->GetBinError(i,j)));
		 t->SetTextSize(0.02);
        	 t->SetTextAlign(22);
        	 t->Draw();
      		}
  	}
	for (int i=1; i<=nbinsX; i++) {
      		for (int j=1; j<=nbinsY; j++) {
        	 auto t = new TText(hsf->GetXaxis()->GetBinCenter(i),
                            hsf->GetYaxis()->GetBinCenter(j)-1.95,
                            Form(" %4.3f +/- %4.3f",hsf->GetBinContent(i,j), hsf->GetBinError(i,j)));
		 t->SetTextSize(0.02);
		 t->SetTextColor(kMagenta+2);
        	 t->SetTextAlign(22);
        	 t->Draw();
      		}
  	}
	auto tlabel= new TText(0.1,94, Form(name.c_str()) );
	tlabel->SetTextSize(0.04);	
	tlabel->Draw();
	auto tleg1= new TText(0.1,91, Form("Data") );
	tleg1->SetTextSize(0.03);	
	tleg1->Draw();
	auto tleg2= new TText(0.45,91, Form("MC") );
	tleg2->SetTextSize(0.03);
	tleg2->SetTextColor(kRed+2);	
	tleg2->Draw();
	auto tleg3= new TText(0.7,91, Form("Data/MC") );
	tleg3->SetTextColor(kMagenta+2);
	tleg3->SetTextSize(0.03);
	tleg3->Draw(); 

	c->SaveAs(("h_"+std::to_string(year)+"_"+std::to_string(rank)+".pdf").c_str());



}
void make2d( E_bank_manager* bankman, E_bank_fit* dat, E_bank_fit* mc, int year){//has a split map
	TH2D* hmc;
	TH2D* hdata;
	TH2D* hsf;

}
int main(){
setTDRStyle();
	E_bank_manager* e1 = new E_bank_manager();
	
	
	//hax, can do this to do electrons! DONT DO THIS
	//it has a random 50/50 chance of not being about to access it's map from parent class ...... error is not reproducible, it is random .. WHAT?!
	//E_bank_fit* id_Zmu_Data = (E_bank_fit*) e1->id_Zmu_Data;
	//E_bank_fit* id_Zmu_MC = (E_bank_fit*) e1->id_Zmu_MC; 
	



	e1->id_Jmu_Data->printMap(e1->id_Jmu_Data->_map16);
	e1->id_Jmu_Data->printMap(e1->id_Jmu_MC->_map16);
	//id_Zmu_Data->printMap(id_Zmu_Data->_map16);
	//id_Zmu_MC->printMap(id_Zmu_MC->_map16);
	//e1->id_Zmu_Data->printMap(e1->id_Zmu_Data->_map16);

/*
	//2016 GSB VL
	make2d(e1,2016,1,"2016 Gold Eff.",
	//id_Zmu_Data, id_Zmu_MC,
	e1->id_Zmu_Data, e1->id_Zmu_MC,	
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);


	make2d(e1,2016,2,"2016 Silver Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);
	
	make2d(e1,2016,3,"2016 Bronze Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);

	//first set of objects are just placeholders that will get sliced
	make2d(e1,2016,0,"2016 V.L. Eff.",
	e1->vl_Zmu_Data, e1->vl_Zmu_MC,
	e1->vl_Zmu_Data, e1->vl_Zmu_MC);

	//2017 GSB VL
	make2d(e1,2017,1,"2017 Gold Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);

	make2d(e1,2017,2,"2017 Silver Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);
	
	make2d(e1,2017,3,"2017 Bronze Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);

	make2d(e1,2017,0,"2017 V.L. Eff.",
	e1->vl_Zmu_Data, e1->vl_Zmu_MC,
	 e1->vl_Zmu_Data, e1->vl_Zmu_MC);


	//2018 GSB VL
	make2d(e1,2018,1,"2018 Gold Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);

	make2d(e1,2018,2,"2018 Silver Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);
	
	make2d(e1,2018,3,"2018 Bronze Eff.",
	e1->id_Zmu_Data, e1->id_Zmu_MC,
	e1->iso_med_Zmu_Data, e1->iso_med_Zmu_MC,
	e1->sip_isomed_Zmu_Data, e1->sip_isomed_Zmu_MC,
	e1->id_Jmu_Data, e1->id_Jmu_MC);

	make2d(e1,2018,0,"2018 V.L. Eff.",
	e1->vl_Zmu_Data, e1->vl_Zmu_MC,
	 e1->vl_Zmu_Data, e1->vl_Zmu_MC);
	 
	*/ 
	   	make2d(e1,2016,1,"2016 Electron Gold Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);

	make2d(e1,2016,2,"2016 Electron Silver Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);
	
	make2d(e1,2016,3,"2016 Electron Bronze Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC);

	make2d(e1,2016,0,"2016 Electron V.L. Eff.",
	e1->vl_Zel_Data, e1->vl_Zel_MC,
	 e1->vl_Zel_Data, e1->vl_Zel_MC);
	
	
	
	 	//2017 GSB VL
	make2d(e1,2017,1,"2017 Electron Gold Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);

	make2d(e1,2017,2,"2017 Electron Silver Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);
	
	make2d(e1,2017,3,"2017 Electron Bronze Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);

	make2d(e1,2017,0,"2017 Electron V.L. Eff.",
	e1->vl_Zel_Data, e1->vl_Zel_MC,
	 e1->vl_Zel_Data, e1->vl_Zel_MC);

    
    make2d(e1,2018,1,"2018 Electron Gold Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);

	make2d(e1,2018,2,"2018 Electron Silver Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);
	
	make2d(e1,2018,3,"2018 Electron Bronze Eff.",
	e1->id_Zel_Data, e1->id_Zel_MC,
	e1->iso_med_Zel_Data, e1->iso_med_Zel_MC,
	e1->sip_isomed_Zel_Data, e1->sip_isomed_Zel_MC
	);

	make2d(e1,2018,0,"2018 Electron V.L. Eff.",
	e1->vl_Zel_Data, e1->vl_Zel_MC,
	 e1->vl_Zel_Data, e1->vl_Zel_MC);




}




