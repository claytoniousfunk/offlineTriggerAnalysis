void plotMe(){

	TFile *f = TFile::Open("rootFiles/PbPbMB_offlineTriggerAnalysis_3Aug22.root");
	
	TString triggerName = "HLT_HIL3Mu5_NHitQ10_tagging_v1";	


	TH1D *m1, *m2;

	f->GetObject("muPt_trigOn",m1);
	f->GetObject("muPt_all",m2);

	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	c1->cd();
	TPad *p1 = new TPad("p1","p1",0,0,1,1);
	p1->SetLogy();
	p1->Draw();
	p1->cd();
	
	m2->GetXaxis()->SetRangeUser(0,30);
	m2->SetTitle("Muon spectra");
	m2->SetLineStyle(7);
	m2->SetStats(0);
	m1->SetStats(0);
	m2->Draw();
	

	m1->Draw("same");
	
	TLegend *leg = new TLegend(0.5,0.5,0.8,0.7);
	leg->SetTextSize(0.04);
	leg->SetBorderSize(0);
	leg->AddEntry(m2,"muons in all events","l");
	leg->AddEntry(m1,"muons in triggered events","l");
	leg->Draw();


	TH1D *r = (TH1D*) m1->Clone("r");
	r->Divide(m1,m2,1,1,"B");

	
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	c2->cd();
	TPad *p2 = new TPad("p2","p2",0,0,1,1);
	p2->Draw();
	p2->cd();
	r->GetXaxis()->SetRangeUser(0,30);
	r->GetYaxis()->SetTitle("ratio");
	r->SetTitle("Trigger turn-on curve");
	r->Draw();
	TLatex *la = new TLatex();
	la->SetTextSize(0.04);
	la->SetTextFont(42);
	la->DrawLatexNDC(0.15,0.7,triggerName);


}
