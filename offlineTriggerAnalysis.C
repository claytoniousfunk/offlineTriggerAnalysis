
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <dirent.h>  
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h>


#include "/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/myProcesses/hiforest/plugin/eventMap_hiForest.h"
#include "/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/JetEnergyCorrections/JetCorrector.h"
#include "/afs/cern.ch/user/c/cbennett/condorSkim/headers/AnalysisSetupV2p0.h"



bool isQualityMuon_soft(double muPt,
		double muEta, 
		int muChi2NDF, 
		double muInnerD0, 
		double muInnerDz, 
		int muPixelHits, 
		int muIsTracker, 
		int muIsGlobal, 
		int muTrkLayers){

	bool result = true;
	if(muPt < 0.0 ||
			TMath::Abs(muEta) > 2.0 ||
			muChi2NDF == -99 ||
			TMath::Abs(muInnerD0) > 0.3 ||
			TMath::Abs(muInnerDz) > 20 ||
			muPixelHits <= 0 ||
			muIsTracker == 0 ||
			muIsGlobal == 0 ||
			muTrkLayers <= 5)
		result = false;

	return result;
}


bool isQualityMuon_hybridSoft(double muPt,
		double muEta, 
		int muChi2NDF, 
		double muInnerD0, 
		double muInnerDz, 
		int muMuonHits, 
		int muPixelHits, 
		int muIsTracker, 
		int muStations, 
		int muTrkLayers){

	bool result = true;
	if(muPt < 0.0 ||
			TMath::Abs(muEta) > 2.0 ||
			muChi2NDF == -99 ||
			muChi2NDF > 10 ||
			TMath::Abs(muInnerD0) > 0.2 ||
			TMath::Abs(muInnerDz) > 0.5 ||
			muMuonHits <= 0 ||
			muPixelHits <= 0 ||
			muIsTracker == 0 ||
			muStations  <= 1 ||
			muTrkLayers <= 5)
		result = false;

	return result;
}

void offlineTriggerAnalysis(
	TString input = "root://cmsxrootd.fnal.gov//store/user/cbennett/PbPbMB_Run2MiniAOD_slimmed_2Aug22/HIMinimumBias0/crab_PbPbMB_Run2MiniAOD_slimmed_2Aug22/220802_155050/0000/HiForestMiniAOD_100.root",
	TString output = "out.root"){


	TH1D *muPt_trigOn[5];
	TH1D *muPt_all[5];

	muPt_trigOn[0] = new TH1D("muPt_trigOn_C0","muPt_trigOn, cent. 0-90%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_trigOn[1] = new TH1D("muPt_trigOn_C1","muPt_trigOn, cent. 0-10%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_trigOn[2] = new TH1D("muPt_trigOn_C2","muPt_trigOn, cent. 10-30%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_trigOn[3] = new TH1D("muPt_trigOn_C3","muPt_trigOn, cent. 30-50%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_trigOn[4] = new TH1D("muPt_trigOn_C4","muPt_trigOn, cent. 50-90%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	
	muPt_all[0] = new TH1D("muPt_all_C0","muPt_all, cent. 0-90%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_all[1] = new TH1D("muPt_all_C1","muPt_all, cent. 0-10%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_all[2] = new TH1D("muPt_all_C2","muPt_all, cent. 10-30%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_all[3] = new TH1D("muPt_all_C3","muPt_all, cent. 30-50%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	muPt_all[4] = new TH1D("muPt_all_C4","muPt_all, cent. 50-90%; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	

	
	cout << "Opening file..." << endl;
	TFile *file = TFile::Open(input);
	cout << "File opened!" << endl;
	auto em = new eventMap(file);
	em->isMC = 0;
	cout << "Initializing variables..." << endl;
	em->init();
	em->loadJet("akCs4PFJetAnalyzer");
	em->loadMuon("ggHiNtuplizer");
	em->loadMuonTrigger("hltanalysis");
	em->loadTrack();
	//em->loadMuonAnalyzer("muonAnalyzer");
	Long64_t NEvents = em->evtTree->GetEntries();
	cout << "Number of events = " << NEvents << endl;
	cout << "Variables initialized!" << endl;

	// event filters
	std::string filters[3] = {"pprimaryVertexFilter","pphfCoincFilter2Th4","pclusterCompatibilityFilter"};
	em->regEventFilter(3,filters);
	

	
	
	// event loop
	int evi_frac = 0;
	for(int evi = 0; evi < NEvents; evi++){
	//for(int evi = 1; evi < 10000; evi++){	

		em->getEvent(evi);

		if((100*evi/NEvents)%5==0 && 100*evi/NEvents > evi_frac){
			cout<<"evt frac: "<<evi_frac<<"%"<<endl;
		}
		evi_frac = 100*evi/NEvents;

		if(em->checkEventFilter()) continue; // event filters
		if(em->vz>15.0) continue;
		if(em->hiBin > 180) continue;



		double w = em->HLT_HIL3Mu5_NHitQ10_tagging_v1_Prescl * 1.0; // mu 5 prescale

		if(w<=0) continue;

		w = 1.0; // set prescale
	
		//cout << "w = " << w << endl;

		//cout << "nMu = " << em->nMu << endl;
		//for(int i = 0; i < em->nInner ; i++){
		//	//cout << "innerNTrkHits = " << em->innerNTrkHits->at(i) << endl;

		//}

		//cout << "nMu = " << em->nMu << " | nInner = " << em->nInner << endl;
		
		for(int m = 0; m < em->nMu; m++){

			
			if(!isQualityMuon_hybridSoft(em->muPt->at(m),
                                                em->muEta->at(m),
                                                em->muChi2NDF->at(m),
                                                em->muInnerD0->at(m),
                                                em->muInnerDz->at(m),
						em->muMuonHits->at(m),
                                                em->muPixelHits->at(m),
                                                em->muIsTracker->at(m),
                                                em->muStations->at(m),
                                                em->muTrkLayers->at(m))) continue; // skip if muon doesnt pass quality cuts
			


			
			
			if(em->HLT_HIL3Mu5_NHitQ10_tagging_v1 == 1){

				muPt_trigOn[0]->Fill(em->muPt->at(m),w);
				if(em->hiBin >= 0 && em->hiBin <=20) muPt_trigOn[1]->Fill(em->muPt->at(m),w);
				if(em->hiBin > 20 && em->hiBin <=60) muPt_trigOn[2]->Fill(em->muPt->at(m),w);
				if(em->hiBin > 60 && em->hiBin <=100) muPt_trigOn[3]->Fill(em->muPt->at(m),w);
				if(em->hiBin > 100 && em->hiBin <=180) muPt_trigOn[4]->Fill(em->muPt->at(m),w);

			}


			muPt_all[0]->Fill(em->muPt->at(m),w);
			if(em->hiBin >= 0 && em->hiBin <=20) muPt_all[1]->Fill(em->muPt->at(m),w);
			if(em->hiBin > 20 && em->hiBin <=60) muPt_all[2]->Fill(em->muPt->at(m),w);
			if(em->hiBin > 60 && em->hiBin <=100) muPt_all[3]->Fill(em->muPt->at(m),w);
			if(em->hiBin > 100 && em->hiBin <=180) muPt_all[4]->Fill(em->muPt->at(m),w);
			

		}

	} // end event loop


auto wf = TFile::Open(output,"recreate");

for(int i = 0; i < 5; i++){
muPt_trigOn[i]->Write();
muPt_all[i]->Write();
}
wf->Close();

return;




}// end program


