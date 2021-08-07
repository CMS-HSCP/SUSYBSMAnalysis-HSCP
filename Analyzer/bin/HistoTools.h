int Color [] = {1,4,2,8,6,9,3,7,5,46,42};
int Marker[] = {20,22,21,23,29,27,2,30,24,25};
int Style [] = {1,2,5,7,9,10,11,12, 13, 14};
int GraphStyle [] = {20, 21, 22, 23, 24, 25, 26, 27};

//=============================================================
//
//     Symmetrize TH1D Histogram
//
//=============================================================

void symmetrizeHisto (TH1D* histo, int mode){
   int limit, shift = 0;
   if (histo->GetNbinsX()%2==0) limit = histo->GetNbinsX()/2;
   else {
      limit = (histo->GetNbinsX()-1)/2; // in that case ignore the middle bin
      shift = 2;
   }
   if (mode==0){
      for (int x=0; x <= limit; x++){
         histo->SetBinContent (x, 0.5*histo->GetBinContent(x) + 0.5*histo->GetBinContent(histo->GetNbinsX() - x));
         histo->SetBinContent (histo->GetNbinsX() - x, histo->GetBinContent(x));
      }
   } else if (mode>0){
      int left  = histo->Integral(0, limit);
	   int right = histo->Integral(limit+shift, histo->GetNbinsX()+1);
      if (mode == 1){ // take the larger half only
         if (right > left){
            for (int x = 0; x <= limit; x++)
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	 } else {
            for (int x = histo->GetNbinsX(); x >= limit; x--)
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	 }
      } else if (mode == 2){ // take the larger of the two opposite bins
         for (int x = 0; x <= histo->GetNbinsX(); x++){
            if (histo->GetBinContent(histo->GetNbinsX() - x) >= histo->GetBinContent(x))
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	    else histo->SetBinContent(histo->GetNbinsX() - x, histo->GetBinContent(x)); 
	 }
      } else if (mode == 3){ // take the smaller half only
         if (right < left){
            for (int x = 0; x <= limit; x++)
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	 } else {
            for (int x = histo->GetNbinsX(); x >= limit; x--)
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	 }
      } else if (mode == 4){ // take the smaller of the two opposite bins
         for (int x = 0; x <= histo->GetNbinsX(); x++){
            if (histo->GetBinContent(histo->GetNbinsX() - x) <= histo->GetBinContent(x))
               histo->SetBinContent(x, histo->GetBinContent(histo->GetNbinsX() - x));
	    else histo->SetBinContent(histo->GetNbinsX() - x, histo->GetBinContent(x)); 
	 }
      }
   }
}

//=============================================================
//
//     Symmetrize TH2D Histogram
//
//=============================================================

void symmetrizeHisto (TH2D* histo, int mode){
   int limit, shift = 0;
   if (histo->GetNbinsX()%2==0) limit = histo->GetNbinsX()/2;
   else {
      limit = (histo->GetNbinsX()-1)/2; // in that case ignore the middle bin
      shift = 2;
   }
   if (mode == 0){ // take (histogram + mirrored histogram)/2
      for (int y=0; y <= histo->GetNbinsY()+1; y++){ // transform the overflow as well
         for (int x=0; x <= limit; x++){
            histo->SetBinContent(x, y, 0.5*histo->GetBinContent(x, y) + 0.5*histo->GetBinContent(histo->GetNbinsX() - x, y));
            histo->SetBinContent(histo->GetNbinsX() - x, y, histo->GetBinContent(x, y));
         }
      }
   } else if (mode > 0){
      int left  = histo->Integral(0, limit, 0, histo->GetNbinsY()+1);
      int right = histo->Integral(limit+shift, histo->GetNbinsX()+1, 0, histo->GetNbinsY()+1);
      if (mode == 1){ // take the larger half only
         if (right > left){
            for (int y=0; y <= histo->GetNbinsY()+1; y++){ // transform the overflow as well
               for (int x = 0; x <= limit; x++)
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
	    }
	 } else {
            for (int y=0; y <= histo->GetNbinsY()+1; y++){ // transform the overflow as well
               for (int x = histo->GetNbinsX(); x >= limit; x--)
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
            }
	 }
      } else if (mode == 2){ // take the larger of the two opposite bins
         for (int y = 0; y <= histo->GetNbinsY()+1; y++){
            for (int x = 0; x <= histo->GetNbinsX(); x++){
               if (histo->GetBinContent(histo->GetNbinsX() - x, y) >= histo->GetBinContent(x, y))
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
               else histo->SetBinContent(histo->GetNbinsX() - x, y, histo->GetBinContent(x, y)); 
	    }
	 }
      } else if (mode == 3){ // take the smaller half only
         if (right < left){
            for (int y=0; y <= histo->GetNbinsY()+1; y++){ // transform the overflow as well
               for (int x = 0; x <= limit; x++)
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
	    }
	 } else {
            for (int y=0; y <= histo->GetNbinsY()+1; y++){ // transform the overflow as well
               for (int x = histo->GetNbinsX(); x >= limit; x--)
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
            }
	 }
      } else if (mode == 4){ // take the smaller of the two opposite bins
         for (int y = 0; y <= histo->GetNbinsY()+1; y++){
            for (int x = 0; x <= histo->GetNbinsX(); x++){
               if (histo->GetBinContent(histo->GetNbinsX() - x, y) <= histo->GetBinContent(x, y))
                  histo->SetBinContent(x, y, histo->GetBinContent(histo->GetNbinsX() - x, y));
               else histo->SetBinContent(histo->GetNbinsX() - x, y, histo->GetBinContent(x, y)); 
	    }
	 }
      }
   }
}

//=============================================================
//
//     Compare weights in |Eta|<0 vs |ETA|>0
//
//=============================================================

void compareForwardToBackwardWeights (TH1D* EtaS, TH1D* EtaB){
   int limit = EtaS->GetNbinsX(),
       bins  = limit/*,
       shift = 1*/; 
   if (limit%2==0) limit = (limit+2) / 2;
   else {limit /= 2; /*shift = 2;*/}
   std::vector <double> eta;              // x value
   std::vector <double> etaSL;            // height of the bins in B region
   std::vector <double> etaSR;            // height of the bins in B region
   std::vector <double> etaBL;            // height of the bins in A region
   std::vector <double> etaBR;            // height of the bins in A region
   std::vector <double> weightsL;         // weights for |Eta| < 0
   std::vector <double> weightsR;         // weights for |Eta| > 0
   std::vector <double> ratioR2L;         // ratio of the weights

   for (int x = limit; x <= bins; x++){
      printf ("limit = %d, bins = %d, x = %d, bins-x+1 = %d etaL = %.2lf etaR = %.2lf\n",
         limit, bins, x, bins-x+1, EtaS->GetBinCenter(x), EtaS->GetBinCenter(bins-x+1));
      eta     .push_back(fabs(EtaS->GetBinCenter(x)));
      etaSL   .push_back(EtaS->GetBinContent(bins-x+1));
      etaBL   .push_back(EtaB->GetBinContent(bins-x+1));
      weightsL.push_back(EtaB->GetBinContent(bins-x+1)>0?((EtaS->GetBinContent(bins-x+1)*EtaB->Integral())/(EtaB->GetBinContent(bins-x+1)*EtaS->Integral())):0);
      etaSR   .push_back(EtaS->GetBinContent(x));
      etaBR   .push_back(EtaB->GetBinContent(x));
      weightsR.push_back(EtaB->GetBinContent(x)>0?((EtaS->GetBinContent(x)*EtaB->Integral())/(EtaB->GetBinContent(x)*EtaS->Integral())):0);
      ratioR2L.push_back(weightsL[weightsL.size()-1]>0?(weightsR[weightsR.size()-1]/weightsL[weightsL.size()-1]):0);
   }

   TGraph* EtaSL    = new TGraph ((int) eta.size(), &eta[0], &etaSL[0]);
   TGraph* EtaSR    = new TGraph ((int) eta.size(), &eta[0], &etaSR[0]);
   TGraph* EtaBL    = new TGraph ((int) eta.size(), &eta[0], &etaBL[0]);
   TGraph* EtaBR    = new TGraph ((int) eta.size(), &eta[0], &etaBR[0]);
   TGraph* WeightsL = new TGraph ((int) eta.size(), &eta[0], &weightsL[0]);
   TGraph* WeightsR = new TGraph ((int) eta.size(), &eta[0], &weightsR[0]);
   TGraph* RatioR2L = new TGraph ((int) eta.size(), &eta[0], &ratioR2L[0]);

   // save all these control graphs
   EtaSL->SaveAs ("EtaSL.root");
   EtaSR->SaveAs ("EtaSR.root");
   EtaBL->SaveAs ("EtaBL.root");
   EtaBR->SaveAs ("EtaBR.root");
   WeightsL->SaveAs ("WeightsL.root");
   WeightsR->SaveAs ("WeightsR.root");
   RatioR2L->SaveAs ("RatioR2L.root");
}


//=============================================================
//
//     Draw a list of TH1 and superimposed them
//
//=============================================================
void DrawSuperposedHistos(TH1** Histos, std::vector<std::string> legend, std::string Style_,  std::string Xlegend, std::string Ylegend, double xmin, double xmax, double ymin, double ymax, bool Normalize=false, bool same=false, bool lastBinOverflow=false, bool firstBinOverflow=false)
{
   int    N             = legend.size();

   double HistoMax      = -1;
   int    HistoHeighest = -1;

   for(int i=0;i<N;i++){
      if(!Histos[i])continue;
      if(Normalize && Histos[i]->Integral()!=0)Histos[i]->Scale(1.0/Histos[i]->Integral());
      Histos[i]->SetTitle("");
      Histos[i]->SetStats(kFALSE);
      Histos[i]->GetXaxis()->SetTitle(Xlegend.c_str());
      Histos[i]->GetYaxis()->SetTitle(Ylegend.c_str());
      Histos[i]->GetXaxis()->SetTitleOffset(1.1);
      Histos[i]->GetYaxis()->SetTitleOffset(1.40);
      Histos[i]->GetXaxis()->SetNdivisions(505);
      Histos[i]->GetYaxis()->SetNdivisions(505);
	   Histos[i]->GetXaxis()->SetTitleSize(0.05);
      if(xmin!=xmax)Histos[i]->SetAxisRange(xmin,xmax,"X");
      if(ymin!=ymax)Histos[i]->SetAxisRange(ymin,ymax,"Y");
      if(ymin==ymax && ymin<0)Histos[i]->SetMaximum(Histos[i]->GetMaximum()*fabs(ymax));
      Histos[i]->SetFillColor(0);
      Histos[i]->SetMarkerStyle(Marker[i]);
      Histos[i]->SetMarkerColor(Color[i]);
      Histos[i]->SetMarkerSize(1.5);
      Histos[i]->SetLineColor(Color[i]);
      Histos[i]->SetLineWidth(2);
      if(lastBinOverflow) {
         if(xmin!=xmax) {
            int lastBin=Histos[i]->GetXaxis()->FindBin(xmax);
            double sum=0;
            double error=0;
            for(int b=lastBin; b<Histos[i]->GetNbinsX()+2; b++) {sum+=Histos[i]->GetBinContent(b); error+=Histos[i]->GetBinError(b)*Histos[i]->GetBinError(b);}
            Histos[i]->SetBinContent(lastBin, sum);
            Histos[i]->SetBinError(lastBin, sqrt(error));
         }
         else {
            Histos[i]->SetBinContent(Histos[i]->GetNbinsX(), Histos[i]->GetBinContent(Histos[i]->GetNbinsX())+Histos[i]->GetBinContent(Histos[i]->GetNbinsX()+1));
            double error=sqrt(pow(Histos[i]->GetBinError(Histos[i]->GetNbinsX()),2)+pow(Histos[i]->GetBinError(Histos[i]->GetNbinsX()+1),2));
            Histos[i]->SetBinError(Histos[i]->GetNbinsX(), error);
         }
      }
      if(firstBinOverflow) {
         if(xmin!=xmax) {
            int firstBin=Histos[i]->GetXaxis()->FindBin(xmin);
            double sum=0;
            double error=0;
            for(int b=0; b<firstBin; b++) {sum+=Histos[i]->GetBinContent(b); error+=Histos[i]->GetBinError(b)*Histos[i]->GetBinError(b);}
            Histos[i]->SetBinContent(firstBin, sum);
            Histos[i]->SetBinError(firstBin, sqrt(error));
         }
         else {
            Histos[i]->SetBinContent(1, Histos[i]->GetBinContent(1)+Histos[i]->GetBinContent(0));
            double error=sqrt(pow(Histos[i]->GetBinError(1),2)+pow(Histos[i]->GetBinError(0),2));
            Histos[i]->SetBinError(1, error);
         }
      }
      if(Style_=="DataMC" && i==0){
         Histos[i]->SetFillColor(0);
         Histos[i]->SetMarkerStyle(20);
         Histos[i]->SetMarkerColor(1);
         Histos[i]->SetMarkerSize(1);
         Histos[i]->SetLineColor(1);
         Histos[i]->SetLineWidth(2);
      }

      if(Histos[i]->GetMaximum() >= HistoMax){
         HistoMax      = Histos[i]->GetMaximum();
         HistoHeighest = i;
      }
   }

   char Buffer[256];
   if(Style_=="DataMC"){
      if(HistoHeighest==0){
         Histos[HistoHeighest]->Draw("E1");
      }else{
         Histos[HistoHeighest]->Draw("HIST");
      }
      for(int i=0;i<N;i++){
           if(i==HistoHeighest)continue;
           if(i==0){
              Histos[i]->Draw("same E1");
           }else{
              Histos[i]->Draw("same");
           }
      }
   }else{
     if(same) {sprintf(Buffer,"same %s",Style_.c_str());
       Histos[HistoHeighest]->Draw(Buffer);}
     else Histos[HistoHeighest]->Draw(Style_.c_str());
      for(int i=0;i<N;i++){
           if(i==HistoHeighest)continue;
           if(Style_!=""){
	     sprintf(Buffer,"same %s",Style_.c_str());
           }else{
              sprintf(Buffer,"same");
           }
           Histos[i]->Draw(Buffer);
      }
   }
}