{


//   TFile* myfile = new TFile ("minbias_template.root");
   TFile* myfile = new TFile ("minbias_template_uncorr_iter1.root");
   TH2D* CHProb_Vs_Path1;
   TH2D* Prob_Vs_Path1;
   myfile->cd();
   CHProb_Vs_Path1= (TH2D*)gROOT->FindObject("Charge_Vs_Path1");
   Prob_Vs_Path1 = (TH2D*) CHProb_Vs_Path1->Clone("Prob_Vs_Path1");
   Prob_Vs_Path1->Reset();
   Prob_Vs_Path1->SetDirectory(0);

   for(int j=0;j<=Prob_Vs_Path1->GetXaxis()->GetNbins()+1;j++){
     double Ni = 0;
     for(int k=0;k<=Prob_Vs_Path1->GetYaxis()->GetNbins()+1;k++){Ni+=CHProb_Vs_Path1->GetBinContent(j,k);}

     for(int k=0;k<=Prob_Vs_Path1->GetYaxis()->GetNbins()+1;k++){
            double tmp = 0;
            for(int l=0;l<=k;l++){ tmp+=CHProb_Vs_Path1->GetBinContent(j,l);}

            if(Ni>0){
               Prob_Vs_Path1->SetBinContent ( j, k, tmp/Ni);
            }else{
               Prob_Vs_Path1->SetBinContent ( j, k, 0);
            }
     }
  }
   TCanvas* c1 = new TCanvas("c2", "c2", 600,600);
   c1->cd();
   CHProb_Vs_Path1->Draw("colz");
   TCanvas* c2 = new TCanvas("c2", "c2", 600,600);
   c2->cd();
   Prob_Vs_Path1->Draw("colz");
}

