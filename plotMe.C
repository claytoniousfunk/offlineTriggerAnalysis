void plotMe(){

	TFile *f = TFile::Open("out.root");

	TH1D *m1, *m2;

	f->GetObject("muPt_trigOn",m1);
	f->GetObject("muPt_all",m2);

	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	c1->cd();
	TPad *p1 = new TPad("p1","p1",0,0,1,1);
	p1->SetLogy();
	p1->Draw();
	p1->cd();
	


	m2->SetLineStyle(7);
	m2->Draw();


	m1->Draw("same");


	TH1D *r = (TH1D*) m1->Clone("r");
	r->Divide(m1,m2,1,1,"B");

	
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	c2->cd();
	TPad *p2 = new TPad("p2","p2",0,0,1,1);
	p2->Draw();
	p2->cd();
	
	r->GetYaxis()->SetTitle("ratio");
	r->Draw();



}
