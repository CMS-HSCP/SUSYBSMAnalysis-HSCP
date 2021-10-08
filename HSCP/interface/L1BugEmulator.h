#ifndef L1BugEmulator_h
#define L1BugEmulator_h

#include "TH2D.h"

//TODO rewrite from scratch and understand where the numbers come from
constexpr double preTrackingChangeL1IntLumi              = 29679.982; //because hardcoded is love, harcoded is life
constexpr double IntegratedLuminosity13TeV16             = 32170.50;  // pb -> not from brilcalc, but from the email, brilcalc numpy error

class L1BugEmulator{
   public:
      L1BugEmulator(){
         double Weight    = preTrackingChangeL1IntLumi/IntegratedLuminosity13TeV16;
         double pTBins[]  = {0.0, 52, 80, 120, 200, 500, 999999};
         double EtaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4, 999999};

         unsigned int NumPtBins  = static_cast<unsigned int> (sizeof(pTBins)/sizeof(double))  -1;
         unsigned int NumEtaBins = static_cast<unsigned int> (sizeof(EtaBins)/sizeof(double)) -1;

         // first period
         h1 = new TH2D ("Inefficiency_Run273158_274093", "Inefficiency_Run273158_274093", NumPtBins, pTBins, NumEtaBins, EtaBins);
         // eta<0.9 -> pass
         h1->SetBinContent(h1->FindBin( 25.0, 0.5), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin( 85.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin(125.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin(255.0, 0.5), 0.88);
         h1->SetBinContent(h1->FindBin(555.0, 0.5), 0.86);
         // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
         h1->SetBinContent(h1->FindBin( 25.0, 1.0), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 1.0), 0.70);
         h1->SetBinContent(h1->FindBin( 85.0, 1.0), 0.68);
         h1->SetBinContent(h1->FindBin(125.0, 1.0), 0.65);
         h1->SetBinContent(h1->FindBin(255.0, 1.0), 0.68);
         h1->SetBinContent(h1->FindBin(555.0, 1.0), 0.00);
         // 1.2<eta<2.1
         h1->SetBinContent(h1->FindBin( 25.0, 1.5), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 1.5), 0.86);
         h1->SetBinContent(h1->FindBin( 85.0, 1.5), 0.87);
         h1->SetBinContent(h1->FindBin(125.0, 1.5), 0.78);
         h1->SetBinContent(h1->FindBin(255.0, 1.5), 0.69);
         h1->SetBinContent(h1->FindBin(555.0, 1.5), 0.00);
         // 2.1<eta<2.4
         h1->SetBinContent(h1->FindBin( 25.0, 2.2), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 2.2), 0.75);
         h1->SetBinContent(h1->FindBin( 85.0, 2.2), 0.80);
         h1->SetBinContent(h1->FindBin(125.0, 2.2), 0.68);
         h1->SetBinContent(h1->FindBin(255.0, 2.2), 0.00);
         h1->SetBinContent(h1->FindBin(555.0, 2.2), 0.00);
         // 2.4<eta is impossible -> remove muon
         h1->SetBinContent(h1->FindBin( 25.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin( 55.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin( 85.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(125.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(255.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(555.0, 3.2), 0.00);

         // second period
         h2 = new TH2D ("Inefficiency_274094_", "Inefficiency_Run_274094_", NumPtBins, pTBins, NumEtaBins, EtaBins);
         // eta<0.9 -> pass
         h2->SetBinContent(h2->FindBin( 25.0, 0.5), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin( 85.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(125.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(255.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(555.0, 0.5), 0.92);
         // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
         h2->SetBinContent(h2->FindBin( 25.0, 1.0), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 1.0), 0.88);
         h2->SetBinContent(h2->FindBin( 85.0, 1.0), 0.88);
         h2->SetBinContent(h2->FindBin(125.0, 1.0), 0.89);
         h2->SetBinContent(h2->FindBin(255.0, 1.0), 0.84);
         h2->SetBinContent(h2->FindBin(555.0, 1.0), 0.00);
         // 1.2<eta<2.1
         h2->SetBinContent(h2->FindBin( 25.0, 1.5), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 1.5), 0.84);
         h2->SetBinContent(h2->FindBin( 85.0, 1.5), 0.84);
         h2->SetBinContent(h2->FindBin(125.0, 1.5), 0.76);
         h2->SetBinContent(h2->FindBin(255.0, 1.5), 0.63);
         h2->SetBinContent(h2->FindBin(555.0, 1.5), 0.00);
         // 2.1<eta<2.4
         h2->SetBinContent(h2->FindBin( 25.0, 2.2), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 2.2), 0.72);
         h2->SetBinContent(h2->FindBin( 85.0, 2.2), 0.76);
         h2->SetBinContent(h2->FindBin(125.0, 2.2), 0.78);
         h2->SetBinContent(h2->FindBin(255.0, 2.2), 0.66);
         h2->SetBinContent(h2->FindBin(555.0, 2.2), 0.00);
         // 2.4<eta is impossible -> remove muon
         h2->SetBinContent(h2->FindBin( 25.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin( 55.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin( 85.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(125.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(255.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(555.0, 3.2), 0.00);
      }

      ~L1BugEmulator(){
//         delete h1;
//         delete h2;
      }

      bool PassesL1Inefficiency(double pT, double Eta){
         if ((rand()%9999*1.0)/10000 < Weight)
            return (((rand()%9999)*1.0/10000) < h1->GetBinContent(h1->FindBin(pT, Eta)));
         else
            return (((rand()%9999)*1.0/10000) < h2->GetBinContent(h2->FindBin(pT, Eta)));
      }

   private:
      TH2D* h1;
      TH2D* h2;
      double Weight;
};


#endif
