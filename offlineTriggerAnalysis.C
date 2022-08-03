
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

bool isQualityMuon(double muPt,
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


	TH1D *muPt_trigOn = new TH1D("muPt_trigOn","muPt_trigOn; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);
	TH1D *muPt_all = new TH1D("muPt_all","muPt_all; muPt [GeV]; Entries",NMuPtBins,muPtMin,muPtMax);

	
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
	em->loadMuonAnalyzer("muonAnalyzer");
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


		//if(em->HLT_HIL3Mu5_NHitQ10_v1==0) continue; // mu5 trigger

		double w = em->HLT_HIL3Mu5_NHitQ10_tagging_v1_Prescl * 1.0; // mu 5 prescale

		if(w<=0) continue;

		//double w = 1.0; // mu 5 prescale
	
		//cout << "w = " << w << endl;

		//cout << "nMu = " << em->nMu << endl;
		for(int i = 0; i < em->nInner ; i++){
			//cout << "innerNTrkHits = " << em->innerNTrkHits->at(i) << endl;

		}

		cout << "nMu = " << em->nMu << " | nInner = " << em->nInner << endl;
		
		for(int m = 0; m < em->nMu; m++){

			//cout << "hello" << endl;
				
			if(!isQualityMuon(em->muPt->at(m),
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

				muPt_trigOn->Fill(em->muPt->at(m),w);	

			}


			
			//cout << "muPt = " << em->muPt->at(m) << " | w = " << w << endl;
			muPt_all->Fill(em->muPt->at(m),w);

		}

	} // end event loop


auto wf = TFile::Open(output,"recreate");


muPt_trigOn->Write();
muPt_all->Write();

wf->Close();

return;




}// end program


