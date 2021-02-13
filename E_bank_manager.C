
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
};


E_bank_manager::E_bank_manager(){
	//hardcoded systematics
	//	

	//big collection of paths
	//muon path
	std::string path = "../";
	//muon filename

	/////MUON ID/////////
	id_Jmu_Data = new E_bank(2016,
				path+"TnPJ_MuonID_data2016_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2017_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_data2018_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	id_Jmu_Data->applySystematics(0.0151);	

	id_Zmu_Data = new E_bank(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	id_Zmu_Data->applySystematics(0.0054);

	id_Jmu_MC = new E_bank(2016,
				path+"TnPJ_MuonID_mc2016_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2017_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				path+"TnPJ_MuonID_mc2018_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	id_Jmu_MC->applySystematics(0.0151);	

	id_Zmu_MC = new E_bank(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	id_Zmu_MC->applySystematics(0.0054);

	/////MUON ISO///////////
	iso_med_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	iso_med_Zmu_Data->applySystematics(0.0057);	
	iso_med_Zmu_Data->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
				
	iso_med_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	iso_med_Zmu_MC->applySystematics(0.0057);
	iso_med_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");

	//////MUON SIP//////////
	sip_isomed_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	sip_isomed_Zmu_Data->applySystematics(0.0057);
	sip_isomed_Zmu_Data->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
				
	sip_isomed_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	sip_isomed_Zmu_MC->applySystematics(0.0057);
	sip_isomed_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");


	/////////MUON VERY LOOSE//////////
	vl_Zmu_Data = new E_bank_fit(2016,
				path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	vl_Zmu_Data->applySystematics(0.0054);
	vl_Zmu_Data->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");

	vl_Zmu_MC = new E_bank_fit(2016,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
	vl_Zmu_MC->applySystematics(0.0054);
	vl_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
		
}
///////////eff def
double E_bank_manager::getMCIdEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
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
		//todo
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
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataIsoEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getMCSipEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getEfficiency(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataSipEfficiency(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getEfficiency(pt,eta,year);
	}
	return -1;
}
/////////////error def
double E_bank_manager::getMCIdError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
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
		//todo
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
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getError(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataIsoError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getError(pt,eta,year);
	}
	return -1;
}	
double E_bank_manager::getMCSipError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getError(pt,eta,year);
	}
	return -1;
}
double E_bank_manager::getDataSipError(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getError(pt,eta,year);
	}
	return -1;
}
////////////////pair def
std::pair<double,double> E_bank_manager::getMCIdPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
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
		//todo
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
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_MC->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getDataIsoPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return iso_med_Zmu_Data->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getMCSipPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_MC->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}
std::pair<double,double> E_bank_manager::getDataSipPair(double pt, double eta, int pdg, int year=0){
	if(abs(pdg) == 11){//electron
		//todo
	}
	if(abs(pdg) == 13){//muon
		return sip_isomed_Zmu_Data->getPair(pt,eta,year);
	}
	return std::make_pair(-1.,-1.);
}

//testing
int main(){
	E_bank_manager* e1 = new E_bank_manager();
	e1->id_Jmu_Data->printMap(e1->id_Jmu_Data->_map16);
	e1->id_Zmu_Data->printMap(e1->id_Zmu_Data->_map16);
	//e1->id_Zmu_Data->applySystematics(0.0054);
	e1->id_Zmu_Data->printMap(e1->id_Zmu_Data->_map16);
	std::cout<<e1->getDataIdEfficiency(6,1.3, 13)<<std::endl;
	std::cout<<e1->getDataIdEfficiency(99,0.0001,13)<<std::endl;
	std::cout<<e1->getDataIdEfficiency(13.99,2.11,13)<<std::endl;
	e1->iso_med_Zmu_Data->printMap(e1->iso_med_Zmu_Data->_map16);
	e1->iso_med_Zmu_Data->printMap(e1->iso_med_Zmu_Data->_fitmap16);
	e1->vl_Zmu_Data->printMap(e1->vl_Zmu_Data->_fitmap16);
	/*std::cout<<"iso testing ... "<<std::endl;
	std::cout<<e1->getMCIsoPair(3.1,0.1,13,2017).first<<" "<<e1->getMCIsoPair(3.1,0.1,13,2017).second<<std::endl;
	std::cout<<e1->getDataSipPair(44.,1.36,13,2018).first<<" "<<e1->getDataSipPair(44.,1.36,13,2018).second<<std::endl;
	e1->iso_med_Zmu_MC->printMap(e1->iso_med_Zmu_MC->_map17);
	e1->iso_med_Zmu_MC->printMap(e1->iso_med_Zmu_MC->_fitmap17);
	e1->sip_isomed_Zmu_Data->printMap(e1->sip_isomed_Zmu_Data->_map18);
	e1->sip_isomed_Zmu_Data->printMap(e1->sip_isomed_Zmu_Data->_fitmap18);
	*/
}


