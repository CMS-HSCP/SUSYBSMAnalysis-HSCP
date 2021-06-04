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
       bins  = limit,
       shift = 1;
   if (limit%2==0) limit = (limit+2) / 2;
   else {limit /= 2; shift = 2;}
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
//     
//
//=============================================================
