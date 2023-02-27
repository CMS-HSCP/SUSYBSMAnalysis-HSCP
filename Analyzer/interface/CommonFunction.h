
#ifndef SUSYBSMAnalysis_Analyzer_CommonFunction_h
#define SUSYBSMAnalysis_Analyzer_CommonFunction_h

//=======================================================================================
//
// Global use of the SaturationCorrection class
// Need to have the class in Analysis_CommonFunction.h & in Analysis_Step1_EventLoop.C
// Step1 uses the functions in CommonFunction.h
// Need to load the correction parameters from a file
//
//=======================================================================================
#include "SaturationCorrection.h"  // New procedure for the correction of the saturation phenomena
SaturationCorrection sc;
void LoadCorrectionParameters() {
  char PathToParameters[2048];
  sprintf(PathToParameters, "%s/src/SUSYBSMAnalysis/HSCP/data/CorrectionParameters.root", getenv("CMSSW_BASE"));
  TFile* fileparameters = TFile::Open(PathToParameters);
  TTree* treeparameters = (TTree*)fileparameters->Get("tree");
  sc.SetTree(*treeparameters);
  sc.ReadParameters();
}

//=======================================================================================
//
//  Global use of some variables used to correct pixel charge
//
//=======================================================================================

#define calmax 50000

   int icalibL1;
   int icalibL1_2017;
   int icalibL1_2018;
   int icalibL2;
   int icalibL2_2017;
   int icalibL2_2018;
   int icalibL3;
   int icalibL3_2017;
   int icalibL3_2018;
   int icalibL4;
   int icalibL4_2017;
   int icalibL4_2018;
   int icalibR1;
   int icalibR1_2017;
   int icalibR1_2018;
   int icalibR2;
   int icalibR2_2017;
   int icalibR2_2018;
   int irunMinValL1[calmax];
   float scaleValL1[calmax];
   float errorScaleValL1[calmax];
   int irunMinValL2[calmax];
   float scaleValL2[calmax];
   float errorScaleValL2[calmax];
   int irunMinValL3[calmax];
   float scaleValL3[calmax];
   float errorScaleValL3[calmax];
   int irunMinValL4[calmax];
   float scaleValL4[calmax];
   float errorScaleValL4[calmax];
   int irunMinValR1[calmax];
   float scaleValR1[calmax];
   float errorScaleValR1[calmax];
   int irunMinValR2[calmax];
   float scaleValR2[calmax];
   float errorScaleValR2[calmax];


//====================================================================
//
//  Get pixel scale factors
//
//====================================================================

void loadSFPixel(int sampleType_) {
   if (sampleType_ > 0) return;
   std::ifstream file_calib1;
   file_calib1.open("CorrFact2017PixL1.txt");
   icalibL1=0;
   if (!file_calib1) { std::cerr << "cannot open file CorrFact2017PixL1.txt " << std::endl; }
   else
   {
      while (!file_calib1.eof () && icalibL1<calmax) {
       file_calib1 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
       icalibL1++;
      }
   }
   file_calib1.close ();

   icalibL1_2017=icalibL1;
   std::ifstream file_calib2;
   file_calib2.open("CorrFact2018PixL1.txt");
   if (!file_calib2) { std::cerr << "cannot open file CorrFact2018PixL1.txt " << std::endl; }
   else
   {
      while (!file_calib2.eof () && icalibL1<calmax) {
       file_calib2 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
       icalibL1++;
      }
   }
   if (irunMinValL1[icalibL1-1]==0) icalibL1--;
   icalibL1_2018=icalibL1-icalibL1_2017;
   file_calib2.close ();


   std::ifstream file_calib3;
   file_calib3.open("CorrFact2017PixL2.txt");
   icalibL2=0;
   if (!file_calib3) { std::cerr << "cannot open file CorrFact2017PixL2.txt " << std::endl; }
   else
   {
      while (!file_calib3.eof () && icalibL2<calmax) {
       file_calib3 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
       icalibL2++;
      }
   }
   file_calib3.close ();

   icalibL2_2017=icalibL2;
   std::ifstream file_calib4;
   file_calib4.open("CorrFact2018PixL2.txt");
   if (!file_calib4) { std::cerr << "cannot open file CorrFact2018PixL2.txt " << std::endl; }
   else
   {
      while (!file_calib4.eof () && icalibL2<calmax) {
       file_calib4 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
       icalibL2++;
      }
   }
   if (irunMinValL2[icalibL2-1]==0) icalibL2--;
   icalibL2_2018=icalibL2-icalibL2_2017;
   file_calib4.close ();

   std::ifstream file_calib5;
   file_calib5.open("CorrFact2017PixL3.txt");
   icalibL3=0;
   if (!file_calib5) { std::cerr << "cannot open file CorrFact2017PixL3.txt " << std::endl; }
   else
   {
      while (!file_calib5.eof () && icalibL3<calmax) {
       file_calib5 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
       icalibL3++;
      }
   }
   file_calib5.close ();

   icalibL3_2017=icalibL3;
   std::ifstream file_calib6;
   file_calib6.open("CorrFact2018PixL3.txt");
   if (!file_calib6) { std::cerr << "cannot open file CorrFact2018PixL3.txt " << std::endl; }
   else
   {
      while (!file_calib6.eof () && icalibL3<calmax) {
       file_calib6 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
       icalibL3++;
      }
   }
   if (irunMinValL3[icalibL3-1]==0) icalibL3--;
   icalibL3_2018=icalibL3-icalibL3_2017;
   file_calib6.close ();

   std::ifstream file_calib7;
   file_calib7.open("CorrFact2017PixL4.txt");
   icalibL4=0;
   if (!file_calib7) { std::cerr << "cannot open file CorrFact2017PixL4.txt " << std::endl; }
   else
   {
      while (!file_calib7.eof () && icalibL4<calmax) {
       file_calib7 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
       icalibL4++;
      }
   }
   file_calib7.close ();

   icalibL4_2017=icalibL4;
   std::ifstream file_calib8;
   file_calib8.open("CorrFact2018PixL4.txt");
   if (!file_calib8) { std::cerr << "cannot open file CorrFact2018PixL4.txt " << std::endl; }
   else
   {
      while (!file_calib8.eof () && icalibL4<calmax) {
       file_calib8 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
       icalibL4++;
      }
   }
   if (irunMinValL4[icalibL4-1]==0) icalibL4--;
   icalibL4_2018=icalibL4-icalibL4_2017;
   file_calib8.close ();

   std::ifstream file_calib9;
   file_calib9.open("CorrFact2017PixR1.txt");
   icalibR1=0;
   if (!file_calib9) { std::cerr << "cannot open file CorrFact2017PixR1.txt " << std::endl; }
   else
   {
      while (!file_calib9.eof () && icalibR1<calmax) {
       file_calib9 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
       icalibR1++;
      }
   }
   file_calib9.close ();

   icalibR1_2017=icalibR1;
   std::ifstream file_calib10;
   file_calib10.open("CorrFact2018PixR1.txt");
   if (!file_calib10) { std::cerr << "cannot open file CorrFact2018PixR1.txt " << std::endl; }
   else
   {
      while (!file_calib10.eof () && icalibR1<calmax) {
       file_calib10 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
       icalibR1++;
      }
   }
   if (irunMinValR1[icalibR1-1]==0) icalibR1--;
   icalibR1_2018=icalibR1-icalibR1_2017;
   file_calib10.close ();

   std::ifstream file_calib11;
   file_calib11.open("CorrFact2017PixR2.txt");
   icalibR2=0;
   if (!file_calib11) { std::cerr << "cannot open file CorrFact2017PixR2.txt " << std::endl; }
   else
   {
      while (!file_calib11.eof () && icalibR2<calmax) {
       file_calib11 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
       icalibR2++;
      }
   }
   file_calib11.close ();

   if (irunMinValR2[icalibR2-1]==0) icalibR2--;
   icalibR2_2017=icalibR2;
   std::ifstream file_calib12;
   file_calib12.open("CorrFact2018PixR2.txt");
   if (!file_calib12) { std::cerr << "cannot open file CorrFact2018PixR2.txt " << std::endl; }
   else
   {
      while (!file_calib12.eof () && icalibR2<calmax) {
       file_calib12 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
       icalibR2++;
      }
   }
   icalibR2_2018=icalibR2-icalibR2_2017;
   file_calib12.close ();
}

//====================================================================
//
//  Get pixel scale factors
//
//====================================================================
float GetSFPixel(int subdetid_, UInt_t detid_, string year, int run) {

    if(year=="") return 1.;
   int pix = 0;
   int layerorside = 0;

   if (subdetid_==1) {
        pix =1 ;
        if (year=="2016") {
                layerorside = int((detid_>>16)&0xF);
        }
        else {
                layerorside = int((detid_>>20)&0xF);
        }
   }
   if (subdetid_==2) {
        pix =2 ;
        layerorside = int((detid_>>23)&0x3); // 1=FPIX- 2=FPIX+
   }

   float scale =1 ;
   if (pix==1) {
     if (layerorside==1) {

       if (year=="2017") {
        for (int i=0; i<icalibL1_2017; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibL1_2017; i<icalibL1-1; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
        if (run >= irunMinValL1[icalibL1-1]) {
                scale = scaleValL1[icalibL1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year=="2017") {
        for (int i=0; i<icalibL2_2017; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibL2_2017; i<icalibL2-1; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
        if (run >= irunMinValL2[icalibL2-1]) {
                scale = scaleValL2[icalibL2-1];
                return scale;
        }
       }

     }
     else if (layerorside==3) {
     }

       if (year=="2017") {
        for (int i=0; i<icalibL3_2017; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibL3_2017; i<icalibL3-1; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
        if (run >= irunMinValL3[icalibL3-1]) {
                scale = scaleValL3[icalibL3-1];
                return scale;
        }
       }

     else if (layerorside==4) {

       if (year=="2017") {
        for (int i=0; i<icalibL4_2017; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibL4_2017; i<icalibL4-1; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
        if (run >= irunMinValL4[icalibL4-1]) {
                scale = scaleValL4[icalibL4-1];
                return scale;
        }
       }

     }
     else {
        std::cout << "unknowm layer number in Barrel Pixel " << layerorside << std::endl;
        return 0;
     }

   }
   else if (pix==2) {
     if (layerorside==1) {

       if (year=="2017") {
        for (int i=0; i<icalibR1_2017; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibR1_2017; i<icalibR1-1; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
        if (run >= irunMinValR1[icalibR1-1]) {
                scale = scaleValR1[icalibR1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year=="2017") {
        for (int i=0; i<icalibR2_2017; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
       }
       else if (year=="2018") {
        for (int i=icalibR2_2017; i<icalibR2-1; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
        if (run >= irunMinValR2[icalibR2-1]) {
                scale = scaleValR2[icalibR2-1];
                return scale;
        }
       }

     }
     else {
        std::cout << "unknowm ring number in EndCap Pixel " << layerorside << std::endl;
        return 0;
     }
   }

   return scale;

}


//====================================================================
//
//  General purpose code
//
//====================================================================

bool startsWith(std::string s, std::string sub) { return s.find(sub) == 0 ? 1 : 0; }

bool endsWith(std::string s, std::string sub) { return s.rfind(sub) == (s.length() - sub.length()) ? 1 : 0; }

void split(std::string s, char delim, std::vector<float>& v) {
  std::string t;
  for (uint i = 0; i < s.length(); i++) {
    if (s[i] == delim) {
      v.push_back(stof(t));
      t = "";
    } else {
      t += s[i];
    }
  }
  if (t != "")
    v.push_back(stof(t));  //include last word
}

void split(std::string s, char delim, std::vector<std::string>& v) {
  std::string t;
  for (uint i = 0; i < s.length(); i++) {
    if (s[i] == delim) {
      v.push_back(t);
      t = "";
    } else {
      t += s[i];
    }
  }
  if (t != "")
    v.push_back(t);  //include last word
}

// create a directory/subdirectory on disk
void MakeDirectories(std::string path) { system((std::string("mkdir -p ") + path).c_str()); }

//=============================================================
//
//    factorial function used to compute probQ
//
//=============================================================

int factorial(int n) { return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n; }

struct HitDeDx {
  float dedx;
  float dx;
  bool isSat;
  bool passClusterCleaning;
  bool isInside;
  unsigned char subDet;
};
typedef std::vector<HitDeDx> HitDeDxCollection;

class dedxHIPEmulator {
private:
  TH1D* ratePdfPixel;
  TH1D* ratePdfStrip;

  float eventRatePixel;
  float eventRateStrip;

  bool is2016;

public:
  void setEventRate(float ratePixel = -1, float rateStrip = -1) {
    if (ratePixel < 0) {
      eventRatePixel = ratePdfPixel->GetRandom();
      eventRatePixel -= 3.2;  //2.4; //subtract rate already present in the MC
      eventRatePixel *= 100;  //for random generator usage
    } else {
      eventRatePixel = ratePixel;
    }
    if (rateStrip < 0) {
      eventRateStrip = ratePdfStrip->GetRandom();
      eventRateStrip -= 1.1;  //0.8; //subtract rate already present in the MC
      eventRateStrip *= 100;  //for random generator usage
    } else {
      eventRateStrip = rateStrip;
    }

    eventRatePixel = std::max(eventRatePixel, 0.f);
    eventRateStrip = std::max(eventRateStrip, 0.f);
  }

  void setPeriodHIPRate(bool is2016_, const char* ratePdfPixelName = nullptr, const char* ratePdfStripName = nullptr) {
    if (is2016 == is2016_ && ratePdfPixel && ratePdfStrip)
      return;
    else {
      is2016 = is2016_;
      char tmp1[1024];
      char tmp2[1024];
      if (ratePdfPixel && !ratePdfPixelName)
        sprintf(tmp1, "%s_%d", ratePdfPixel->GetName(), is2016 ? 2016 : 2015);
      else if (ratePdfPixelName)
        sprintf(tmp1, "%s_%d", ratePdfPixelName, is2016 ? 2016 : 2015);
      else
        sprintf(tmp1, "ratePdfPixel_%d", is2016 ? 2016 : 2015);

      if (ratePdfStrip && !ratePdfStripName)
        sprintf(tmp2, "%s_%d", ratePdfStrip->GetName(), is2016 ? 2016 : 2015);
      else if (ratePdfStripName)
        sprintf(tmp2, "%s_%d", ratePdfStripName, is2016 ? 2016 : 2015);
      else
        sprintf(tmp2, "ratePdfStrip_%d", is2016 ? 2016 : 2015);

      std::cout << "Histo1 name is " << tmp1 << std::endl;
      std::cout << "Histo2 name is " << tmp2 << std::endl;
      if (ratePdfPixel) {
        ratePdfPixel->SetBins(is2016 ? 20 : 26, 0, is2016 ? 10 : 13);
        ratePdfPixel->SetNameTitle(tmp1, tmp1);
      } else
        ratePdfPixel = new TH1D(tmp1, tmp1, is2016_ ? 20 : 26, 0, is2016_ ? 10 : 13);
      if (ratePdfStrip) {
        ratePdfStrip->SetBins(is2016 ? 20 : 26, 0, is2016 ? 10 : 13);
        ratePdfStrip->SetNameTitle(tmp2, tmp2);
      } else
        ratePdfStrip = new TH1D(tmp2, tmp2, is2016_ ? 20 : 26, 0, is2016_ ? 10 : 13);
      if (!is2016_) {
        ratePdfPixel->SetBinContent(2, 116789);
        ratePdfPixel->SetBinContent(3, 2501600);
        ratePdfPixel->SetBinContent(4, 7088437);
        ratePdfPixel->SetBinContent(5, 3.771881e+07);
        ratePdfPixel->SetBinContent(6, 2.24745e+07);
        ratePdfPixel->SetBinContent(7, 7.324569e+07);
        ratePdfPixel->SetBinContent(8, 1.703727e+07);
        ratePdfPixel->SetBinContent(9, 1.72861e+07);
        ratePdfPixel->SetBinContent(10, 1.059352e+07);
        ratePdfPixel->SetBinContent(11, 2.179018e+07);
        ratePdfPixel->SetBinContent(12, 7.105593e+07);
        ratePdfPixel->SetBinContent(13, 1495028);
        ratePdfPixel->SetBinContent(14, 163321);
        ratePdfPixel->SetBinContent(15, 140044);
        ratePdfPixel->SetBinContent(16, 548);
        ratePdfPixel->Scale(1.0 / ratePdfPixel->Integral());

        ratePdfStrip->SetBinContent(4, 122);
        ratePdfStrip->SetBinContent(5, 312);
        ratePdfStrip->SetBinContent(6, 1848);
        ratePdfStrip->SetBinContent(7, 21443);
        ratePdfStrip->SetBinContent(8, 408057);
        ratePdfStrip->SetBinContent(9, 3.178217e+07);
        ratePdfStrip->SetBinContent(10, 1.560515e+08);
        ratePdfStrip->SetBinContent(11, 9.134156e+07);
        ratePdfStrip->SetBinContent(12, 1306459);
        ratePdfStrip->SetBinContent(13, 1759553);
        ratePdfStrip->SetBinContent(14, 29769);
        ratePdfStrip->SetBinContent(15, 3097);
        ratePdfStrip->SetBinContent(16, 1322);
        ratePdfStrip->SetBinContent(17, 177);
        ratePdfStrip->SetBinContent(18, 400);
        ratePdfStrip->Scale(1.0 / ratePdfStrip->Integral());
      } else {
        ratePdfPixel->SetBinContent(1, 72346);
        ratePdfPixel->SetBinContent(2, 55300);
        ratePdfPixel->SetBinContent(3, 703500);
        ratePdfPixel->SetBinContent(4, 1258900);
        ratePdfPixel->SetBinContent(5, 2542300);
        ratePdfPixel->SetBinContent(6, 9.472088e+07);
        ratePdfPixel->SetBinContent(7, 6.0498e+07);
        ratePdfPixel->SetBinContent(8, 8.668988e+08);
        ratePdfPixel->SetBinContent(9, 1.311922e+09);
        ratePdfPixel->SetBinContent(10, 1.682412e+09);
        ratePdfPixel->SetBinContent(11, 1.377868e+09);
        ratePdfPixel->SetBinContent(12, 1.369749e+09);
        ratePdfPixel->SetBinContent(13, 1.350122e+09);
        ratePdfPixel->SetBinContent(14, 1.153147e+09);
        ratePdfPixel->SetBinContent(15, 2.164707e+09);
        ratePdfPixel->SetBinContent(16, 7.892244e+08);
        ratePdfPixel->SetBinContent(17, 5.456613e+08);
        ratePdfPixel->SetBinContent(18, 8.04757e+07);
        ratePdfPixel->SetBinContent(19, 3642000);
        ratePdfPixel->SetBinContent(20, 595200);
        ratePdfPixel->SetBinContent(23, 34700);
        ratePdfPixel->SetBinContent(24, 218000);
        ratePdfPixel->Scale(1.0 / ratePdfPixel->Integral());

        ratePdfStrip->SetBinContent(1, 28146);
        ratePdfStrip->SetBinContent(2, 255300);
        ratePdfStrip->SetBinContent(3, 3341120);
        ratePdfStrip->SetBinContent(4, 1.383308e+08);
        ratePdfStrip->SetBinContent(5, 1.469652e+08);
        ratePdfStrip->SetBinContent(6, 8.804663e+07);
        ratePdfStrip->SetBinContent(7, 3.27495e+07);
        ratePdfStrip->SetBinContent(8, 2.400779e+08);
        ratePdfStrip->SetBinContent(9, 1.818388e+08);
        ratePdfStrip->SetBinContent(10, 4.019324e+08);
        ratePdfStrip->SetBinContent(11, 4.590071e+08);
        ratePdfStrip->SetBinContent(12, 6.454876e+08);
        ratePdfStrip->SetBinContent(13, 7.114447e+08);
        ratePdfStrip->SetBinContent(14, 9.504209e+08);
        ratePdfStrip->SetBinContent(15, 1.033392e+09);
        ratePdfStrip->SetBinContent(16, 6.776536e+08);
        ratePdfStrip->SetBinContent(17, 8.913748e+08);
        ratePdfStrip->SetBinContent(18, 1.203708e+09);
        ratePdfStrip->SetBinContent(19, 1.42717e+09);
        ratePdfStrip->SetBinContent(20, 8.380354e+08);
        ratePdfStrip->SetBinContent(21, 9.456285e+08);
        ratePdfStrip->SetBinContent(22, 8.246163e+08);
        ratePdfStrip->SetBinContent(23, 5.379617e+08);
        ratePdfStrip->SetBinContent(24, 2.384386e+08);
        ratePdfStrip->SetBinContent(25, 1.229594e+08);
        ratePdfStrip->SetBinContent(26, 1.087633e+08);
        ratePdfStrip->SetBinContent(27, 1.37972e+07);
        ratePdfStrip->Scale(1.0 / ratePdfStrip->Integral());
      }
    }
  }

  dedxHIPEmulator(bool is2016_ = false,
                  const char* ratePdfPixelName = nullptr,
                  const char* ratePdfStripName = nullptr) {
    setPeriodHIPRate(is2016_, ratePdfPixelName, ratePdfStripName);
  }
  ~dedxHIPEmulator() {}

  float getEventRatePixel() { return eventRatePixel; }
  float getEventRateStrip() { return eventRateStrip; }

  float fakeHIP(unsigned int subDet, float dedx) {
    if (subDet < 3 && rand() % 10000 < eventRatePixel)
      dedx = (0.6 + ((rand() % 15000) / 10000.0));
    if (subDet >= 3 && rand() % 10000 < eventRateStrip)
      dedx = (0.6 + ((rand() % 15000) / 10000.0));
    return dedx;
  }

  void fakeHIP(HitDeDxCollection& hitDeDx) {
    for (unsigned int h = 0; h < hitDeDx.size(); h++) {
      hitDeDx[h].dedx = fakeHIP(hitDeDx[h].subDet, hitDeDx[h].dedx);
    }
  }
};

class L1BugEmulator {
private:
  TH2D* h1;
  TH2D* h2;
  float Weight;

public:
  L1BugEmulator(float preTrackingChangeL1IntLumi = 1., float IntegratedLuminosity = 1.) {
    //UNUSED//float Weight    = preTrackingChangeL1IntLumi/IntegratedLuminosity;
    /*?maybe*/ Weight = preTrackingChangeL1IntLumi / IntegratedLuminosity;
    float pTBins[] = {0.0, 52, 80, 120, 200, 500, 999999};
    float EtaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4, 999999};

    unsigned int NumPtBins = static_cast<unsigned int>(sizeof(pTBins) / sizeof(float)) - 1;
    unsigned int NumEtaBins = static_cast<unsigned int>(sizeof(EtaBins) / sizeof(float)) - 1;

    // first period
    h1 = new TH2D(
        "Inefficiency_Run273158_274093", "Inefficiency_Run273158_274093", NumPtBins, pTBins, NumEtaBins, EtaBins);
    // eta<0.9 -> pass
    h1->SetBinContent(h1->FindBin(25.0, 0.5), 0.00);  // we do not trigger on muons with pT<50
    h1->SetBinContent(h1->FindBin(55.0, 0.5), 0.92);
    h1->SetBinContent(h1->FindBin(85.0, 0.5), 0.92);
    h1->SetBinContent(h1->FindBin(125.0, 0.5), 0.92);
    h1->SetBinContent(h1->FindBin(255.0, 0.5), 0.88);
    h1->SetBinContent(h1->FindBin(555.0, 0.5), 0.86);
    // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
    h1->SetBinContent(h1->FindBin(25.0, 1.0), 0.00);  // we do not trigger on muons with pT<50
    h1->SetBinContent(h1->FindBin(55.0, 1.0), 0.70);
    h1->SetBinContent(h1->FindBin(85.0, 1.0), 0.68);
    h1->SetBinContent(h1->FindBin(125.0, 1.0), 0.65);
    h1->SetBinContent(h1->FindBin(255.0, 1.0), 0.68);
    h1->SetBinContent(h1->FindBin(555.0, 1.0), 0.00);
    // 1.2<eta<2.1
    h1->SetBinContent(h1->FindBin(25.0, 1.5), 0.00);  // we do not trigger on muons with pT<50
    h1->SetBinContent(h1->FindBin(55.0, 1.5), 0.86);
    h1->SetBinContent(h1->FindBin(85.0, 1.5), 0.87);
    h1->SetBinContent(h1->FindBin(125.0, 1.5), 0.78);
    h1->SetBinContent(h1->FindBin(255.0, 1.5), 0.69);
    h1->SetBinContent(h1->FindBin(555.0, 1.5), 0.00);
    // 2.1<eta<2.4
    h1->SetBinContent(h1->FindBin(25.0, 2.2), 0.00);  // we do not trigger on muons with pT<50
    h1->SetBinContent(h1->FindBin(55.0, 2.2), 0.75);
    h1->SetBinContent(h1->FindBin(85.0, 2.2), 0.80);
    h1->SetBinContent(h1->FindBin(125.0, 2.2), 0.68);
    h1->SetBinContent(h1->FindBin(255.0, 2.2), 0.00);
    h1->SetBinContent(h1->FindBin(555.0, 2.2), 0.00);
    // 2.4<eta is impossible -> remove muon
    h1->SetBinContent(h1->FindBin(25.0, 3.2), 0.00);
    h1->SetBinContent(h1->FindBin(55.0, 3.2), 0.00);
    h1->SetBinContent(h1->FindBin(85.0, 3.2), 0.00);
    h1->SetBinContent(h1->FindBin(125.0, 3.2), 0.00);
    h1->SetBinContent(h1->FindBin(255.0, 3.2), 0.00);
    h1->SetBinContent(h1->FindBin(555.0, 3.2), 0.00);

    // second period
    h2 = new TH2D("Inefficiency_274094_", "Inefficiency_Run_274094_", NumPtBins, pTBins, NumEtaBins, EtaBins);
    // eta<0.9 -> pass
    h2->SetBinContent(h2->FindBin(25.0, 0.5), 0.00);  // we do not trigger on muons with pT<50
    h2->SetBinContent(h2->FindBin(55.0, 0.5), 0.92);
    h2->SetBinContent(h2->FindBin(85.0, 0.5), 0.92);
    h2->SetBinContent(h2->FindBin(125.0, 0.5), 0.92);
    h2->SetBinContent(h2->FindBin(255.0, 0.5), 0.92);
    h2->SetBinContent(h2->FindBin(555.0, 0.5), 0.92);
    // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
    h2->SetBinContent(h2->FindBin(25.0, 1.0), 0.00);  // we do not trigger on muons with pT<50
    h2->SetBinContent(h2->FindBin(55.0, 1.0), 0.88);
    h2->SetBinContent(h2->FindBin(85.0, 1.0), 0.88);
    h2->SetBinContent(h2->FindBin(125.0, 1.0), 0.89);
    h2->SetBinContent(h2->FindBin(255.0, 1.0), 0.84);
    h2->SetBinContent(h2->FindBin(555.0, 1.0), 0.00);
    // 1.2<eta<2.1
    h2->SetBinContent(h2->FindBin(25.0, 1.5), 0.00);  // we do not trigger on muons with pT<50
    h2->SetBinContent(h2->FindBin(55.0, 1.5), 0.84);
    h2->SetBinContent(h2->FindBin(85.0, 1.5), 0.84);
    h2->SetBinContent(h2->FindBin(125.0, 1.5), 0.76);
    h2->SetBinContent(h2->FindBin(255.0, 1.5), 0.63);
    h2->SetBinContent(h2->FindBin(555.0, 1.5), 0.00);
    // 2.1<eta<2.4
    h2->SetBinContent(h2->FindBin(25.0, 2.2), 0.00);  // we do not trigger on muons with pT<50
    h2->SetBinContent(h2->FindBin(55.0, 2.2), 0.72);
    h2->SetBinContent(h2->FindBin(85.0, 2.2), 0.76);
    h2->SetBinContent(h2->FindBin(125.0, 2.2), 0.78);
    h2->SetBinContent(h2->FindBin(255.0, 2.2), 0.66);
    h2->SetBinContent(h2->FindBin(555.0, 2.2), 0.00);
    // 2.4<eta is impossible -> remove muon
    h2->SetBinContent(h2->FindBin(25.0, 3.2), 0.00);
    h2->SetBinContent(h2->FindBin(55.0, 3.2), 0.00);
    h2->SetBinContent(h2->FindBin(85.0, 3.2), 0.00);
    h2->SetBinContent(h2->FindBin(125.0, 3.2), 0.00);
    h2->SetBinContent(h2->FindBin(255.0, 3.2), 0.00);
    h2->SetBinContent(h2->FindBin(555.0, 3.2), 0.00);
  }

  ~L1BugEmulator() {
    // delete h1;
    // delete h2;
  }

  bool PassesL1Inefficiency(float pT, float Eta) {
    if ((rand() % 9999 * 1.0) / 10000 < Weight)
      return (((rand() % 9999) * 1.0 / 10000) < h1->GetBinContent(h1->FindBin(pT, Eta)));
    else
      return (((rand() % 9999) * 1.0 / 10000) < h2->GetBinContent(h2->FindBin(pT, Eta)));
  }
};

#ifdef FWCORE
class HIPTrackLossEmulator {
private:
  TH1D* h;
  float lossRate;

public:
  HIPTrackLossEmulator() {
    float NVertexBins[] = {0, 5, 10, 15, 22, 28, 35, 9999};
    h = new TH1D("TrackLoss_vs_NVertex",
                 "TrackLoss_vs_NVertex",
                 static_cast<unsigned int>(sizeof(NVertexBins) / sizeof(float)) - 1,
                 NVertexBins);
    h->SetBinContent(h->FindBin(2.5), 1.000);
    h->SetBinContent(h->FindBin(7.5), 0.998);
    h->SetBinContent(h->FindBin(12.5), 0.995);
    h->SetBinContent(h->FindBin(17.5), 0.992);
    h->SetBinContent(h->FindBin(24.5), 0.985);
    h->SetBinContent(h->FindBin(30.0), 0.980);
    h->SetBinContent(h->FindBin(37.0), 0.970);

    lossRate = 0.0;
  }

  ~HIPTrackLossEmulator() {
    //        delete h;
  }

  void SetHIPTrackLossRate(const edm::Event& iEvent,
                           edm::EDGetTokenT<std::vector<reco::Vertex>> offlinePrimaryVerticesToken) {
    std::vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken);
    if (vertexColl.size() < 1) {
      lossRate = 1.0;
      return;
    } else
      lossRate = 1.0 - h->GetBinContent(h->FindBin(vertexColl.size()));
  }

  bool TrackSurvivesHIPInefficiency() { return (((rand() % 999999) * 1.0 / 1000000) > lossRate); }
};
#endif  //FWCORE

//=============================================================
//      Common functions
//=============================================================

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy = false) {
  size_t pos = Path.find("/");
  if (pos < 256) {
    std::string firstPart = Path.substr(0, pos);
    std::string endPart = Path.substr(pos + 1, Path.length());
    TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
    if (TMP != nullptr)
      return GetObjectFromPath(TMP, endPart, GetACopy);

    printf("ObjectNotFound: %s::%s\n", File->GetName(), Path.c_str());
    return nullptr;
  } else {
    if (GetACopy) {
      return (File->Get(Path.c_str()))->Clone();
    } else {
      return File->Get(Path.c_str());
    }
  }
}

// similar to the above code
TObject* GetObjectFromPath(TDirectory* Container, TDirectory* File, std::string Path, bool GetACopy = false) {
  TObject* toreturn = GetObjectFromPath(File, Path, GetACopy);
  if (TH1* th1 = dynamic_cast<TH1*>(toreturn))
    th1->SetDirectory(Container);
  return toreturn;
}

TH1D* GetProjectionFromPath(TDirectory* File, std::string Path, int CutIndex, std::string Name) {
  TH2D* tmp = (TH2D*)GetObjectFromPath(File, Path, false);
  if (!tmp)
    return nullptr;
  return tmp->ProjectionY(Name.c_str(), CutIndex + 1, CutIndex + 1, "o");
}

// return the TypeMode from a string inputPattern
int TypeFromPattern(const std::string& InputPattern) {
  if (InputPattern.find("Type0", 0) < std::string::npos) {
    return 0;
  } else if (InputPattern.find("Type1", 0) < std::string::npos) {
    return 1;
  } else if (InputPattern.find("Type2", 0) < std::string::npos) {
    return 2;
  } else if (InputPattern.find("Type3", 0) < std::string::npos) {
    return 3;
  } else if (InputPattern.find("Type4", 0) < std::string::npos) {
    return 4;
  } else if (InputPattern.find("Type5", 0) < std::string::npos) {
    return 5;
  } else {
    return 6;
  }
}

// define the legend corresponding to a Type
std::string LegendFromType(const std::string& InputPattern) {
  switch (TypeFromPattern(InputPattern)) {
    case 0:
      return std::string("Tracker - Only");
      break;
    case 1:
      return std::string("Tracker + Muon");
      break;
    case 2:
      return std::string("Tracker + TOF");
      break;
    case 3:
      return std::string("Muon - Only");
      break;
    case 4:
      return std::string("|Q|>1e");
      break;
    case 5:
      return std::string("|Q|<1e");
      break;
    default:
      std::string("unknown");
  }
  return std::string("unknown");
}

std::vector<int> convert(const std::vector<unsigned char>& input) {
  std::vector<int> output;
  for (unsigned int i = 0; i < input.size(); i++) {
    output.push_back((int)input[i]);
  }
  return output;
}

// compute deltaR between two point (eta,phi) (eta,phi)
float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  while (dphi > M_PI)
    dphi -= 2 * M_PI;
  while (dphi <= -M_PI)
    dphi += 2 * M_PI;
  return sqrt(deta * deta + dphi * dphi);
}

#ifdef FWCORE
// compute the distance between a "reconstructed" HSCP candidate and the closest generated HSCP
float DistToHSCP(const susybsm::HSCParticle& hscp,
                  const std::vector<reco::GenParticle>& genColl,
                  int& IndexOfClosest,
                  int TypeMode) {
  reco::TrackRef track;
  if (TypeMode != 3)
    track = hscp.trackRef();
  else {
    reco::MuonRef muon = hscp.muonRef();
    if (muon.isNull())
      return false;
    track = muon->standAloneMuon();
  }
  if (track.isNull())
    return false;

  float RMin = 9999;
  IndexOfClosest = -1;
  for (unsigned int g = 0; g < genColl.size(); g++) {
    if (genColl[g].pt() < 5)
      continue;
    if (genColl[g].status() != 1)
      continue;
    int AbsPdg = abs(genColl[g].pdgId());
    if (AbsPdg < 1000000 && AbsPdg != 17)
      continue;
    float dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
    if (dR < RMin) {
      RMin = dR;
      IndexOfClosest = g;
    }
  }
  return RMin;
}
#endif  //FWCORE

std::vector<int> Correction(const std::vector<int>& Q,
                            const int& label,
                            const float& rsat,
                            const float& thresholdSat,
                            const float& thresholdDeltaQ,
                            const float& thresholdrsat) {
  const unsigned N = Q.size();
  if (N < 3)
    return Q;
  std::vector<int> Qcorr;
  int total_charge = 0;
  int nsat = 0;
  float Qmin = 0.;
  float Qplus = 0.;
  bool testQmin = true;
  for (unsigned int i = 0; i < N; i++) {
    total_charge += Q[i];
    if (Q[i] > 253)
      nsat++;
    if (Q[i] > 253 && testQmin && i > 1) {
      testQmin = false;
      Qmin = Q[i - 1];
    }
    if (Q[i] > 253) {
      Qplus = Q[i + 1];
    }
    Qcorr.push_back(Q[i]);
  }
  float DeltaQ = abs(Qmin - Qplus);
  float charge_corr = sc.ChargeCorrected(total_charge, label, N, nsat);
  float DeltaChargeCorr = charge_corr - total_charge;
  if (rsat >= thresholdrsat && Qmin >= thresholdSat && Qplus >= thresholdSat && DeltaQ <= thresholdDeltaQ) {
    if (N == 1) {
      vector<int>::iterator maxQ = max_element(Qcorr.begin(), Qcorr.end());
      Qcorr[std::distance(Qcorr.begin(), maxQ)] += DeltaChargeCorr;
    }
    if (N == 2 && label < 5) {
      vector<int>::iterator maxQ = max_element(Qcorr.begin(), Qcorr.end());
      if (Qcorr[std::distance(Qcorr.begin(), maxQ + 1)] < 254)
        return Qcorr;
      Qcorr[std::distance(Qcorr.begin(), maxQ)] += DeltaChargeCorr / 2;
      Qcorr[std::distance(Qcorr.begin(), maxQ + 1)] += DeltaChargeCorr / 2;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////
  //// if a no saturated strip is between two saturated strips the method return the cluster without correction
  ////
  //// if more than one saturated strip the correction method is only used in the barrel (cut on label<5)
  ////
  //// the difference between the charge & the corrected charge is shared with the saturated strips
  ////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  return Qcorr;
}



std::vector<int> SaturationCorrection(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);

//---  only for one max well-defined
 if(Q.size()<2 || Q.size()>8){
        for (unsigned int i=0;i<Q.size();i++){
                QII.push_back((int) Q[i]);
        }
        return QII;
  }
 if(way){
          vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
          if(*mQ>253){
                 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
                 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
                     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
          }
      else{
          return Q; // no saturation --> no x-talk inversion
      }
  }
//---
 // do nothing else
 return Q;
}

std::vector<int> CrossTalkInv(const std::vector<int>& Q,
                              const float x1 = 0.10,
                              const float x2 = 0.04,
                              bool way = true,
                              float threshold = 20,
                              float thresholdSat = 25,
                              bool isClusterCleaning = false) {
  const unsigned N = Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  //  bool debugbool=false;
  TMatrix A(N, N);

  //---
  if (Q.size() < 2 || Q.size() > 8) {
    for (unsigned int i = 0; i < Q.size(); i++) {
      QII.push_back((int)Q[i]);
    }
    return QII;
  }

  if(way){
      std::vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end()) ;
      if(*mQ>253){
         if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
         if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
             QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
      }
   }
  //---

  for (unsigned int i = 0; i < N; i++) {
    A(i, i) = a;
    if (i < N - 1) {
      A(i + 1, i) = x1;
      A(i, i + 1) = x1;
    } else
      continue;
    if (i < N - 2) {
      A(i + 2, i) = x2;
      A(i, i + 2) = x2;
    }
  }

  if (N == 1)
    A(0, 0) = 1 / a;
  else
    A.InvertFast();

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      QI[i] += A(i, j) * (float)Q[j];
    }
  }

  for (unsigned int i = 0; i < QI.size(); i++) {
    if (QI[i] < threshold)
      QI[i] = 0;
    QII.push_back((int)QI[i]);
  }

  return QII;
}

#ifdef FWCORE
bool clusterCleaning(std::vector<int> ampls, int crosstalkInv = 0, uint8_t* exitCode = nullptr) {


  // ---------------- Count the number of maximas    --------------------------
  //----------------------------------------------------------------------------
  Int_t NofMax = 0;
  Int_t recur255 = 1;
  Int_t recur254 = 1;
  bool MaxOnStart = false;
  bool MaxInMiddle = false, MaxOnEnd = false;
  Int_t MaxPos = 0;
  // Start with a max
  if (ampls.size() != 1 &&
      ((ampls[0] > ampls[1]) ||
       (ampls.size() > 2 && ampls[0] == ampls[1] && ampls[1] > ampls[2] && ampls[0] != 254 && ampls[0] != 255) ||
       (ampls.size() == 2 && ampls[0] == ampls[1] && ampls[0] != 254 && ampls[0] != 255))) {
    NofMax = NofMax + 1;
    MaxOnStart = true;
  }

  // Max in the middle (strip between other ones)
  if (ampls.size() > 2) {
    for (unsigned int i = 1; i < ampls.size() - 1; i++) {
      if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
          (ampls.size() > 3 && i > 0 && i < ampls.size() - 2 && ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
           ampls[i] > ampls[i + 2] && ampls[i] != 254 && ampls[i] != 255)) {
        NofMax = NofMax + 1;
        MaxInMiddle = true;
        MaxPos = i;
      }
      if (ampls[i] == 255 && ampls[i] == ampls[i - 1]) {
        recur255 = recur255 + 1;
        MaxPos = i - (recur255 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
      if (ampls[i] == 254 && ampls[i] == ampls[i - 1]) {
        recur254 = recur254 + 1;
        MaxPos = i - (recur254 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
    }
  }
  // Max at the end of the cluster
  if (ampls.size() > 1) {
    if (ampls[ampls.size() - 1] > ampls[ampls.size() - 2] ||
        (ampls.size() > 2 && ampls[ampls.size() - 1] == ampls[ampls.size() - 2] &&
         ampls[ampls.size() - 2] > ampls[ampls.size() - 3]) ||
        ampls[ampls.size() - 1] == 255) {
      NofMax = NofMax + 1;
      MaxOnEnd = true;
    }
  }
  // If only one strip is hit
  if (ampls.size() == 1) {
    NofMax = 1;
  }

  // --- SHAPE SELECTION
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //
  //   bool shapetest=true;
  bool shapecdtn = false;
  if (exitCode)
    *exitCode = 255;

  if (crosstalkInv == 1) {
    if (NofMax == 1) {
      shapecdtn = true;
      if (exitCode)
        *exitCode = 0;
    }
    return shapecdtn;
  }

  //      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
  Float_t C_M = 0.0;
  Float_t C_D = 0.0;
  Float_t C_Mn = 10000;
  Float_t C_Dn = 10000;
  Float_t C_Mnn = 10000;
  Float_t C_Dnn = 10000;
  Int_t CDPos;
  Float_t coeff1 = 1.7;
  Float_t coeff2 = 2.0;
  Float_t coeffn = 0.10;
  Float_t coeffnn = 0.02;
  Float_t noise = 4.0;

  if (NofMax == 1) {
    if (MaxOnStart == true) {
      C_M = (Float_t)ampls[0];
      C_D = (Float_t)ampls[1];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[2];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 2;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[2];
        C_Dnn = (Float_t)ampls[3];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 3;
      }
    }

    if (MaxOnEnd == true) {
      C_M = (Float_t)ampls[ampls.size() - 1];
      C_D = (Float_t)ampls[ampls.size() - 2];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[0];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 4;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[ampls.size() - 3];
        C_Dnn = (Float_t)ampls[ampls.size() - 4];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 5;
      }
    }

    if (MaxInMiddle == true) {
      C_M = (Float_t)ampls[MaxPos];
      int LeftOfMaxPos = MaxPos - 1;
      if (LeftOfMaxPos <= 0)
        LeftOfMaxPos = 0;
      int RightOfMaxPos = MaxPos + 1;
      if (RightOfMaxPos >= (int)ampls.size())
        RightOfMaxPos = ampls.size() - 1;
      //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
      if (ampls[LeftOfMaxPos] < ampls[RightOfMaxPos]) {
        C_D = (Float_t)ampls[RightOfMaxPos];
        C_Mn = (Float_t)ampls[LeftOfMaxPos];
        CDPos = RightOfMaxPos;
      } else {
        C_D = (Float_t)ampls[LeftOfMaxPos];
        C_Mn = (Float_t)ampls[RightOfMaxPos];
        CDPos = LeftOfMaxPos;
      }
      if (C_Mn < coeff1 * coeffn * C_M + coeff2 * coeffnn * C_D + 2 * noise || C_M == 255) {
        if (ampls.size() == 3)
          shapecdtn = true;
        else if (ampls.size() > 3) {
          if (CDPos > MaxPos) {
            if (ampls.size() - CDPos - 1 == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 == 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 > 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = (Float_t)ampls[CDPos + 2];
            }
            if (MaxPos >= 2) {
              C_Mnn = (Float_t)ampls[MaxPos - 2];
            } else if (MaxPos < 2)
              C_Mnn = 0;
          }
          if (CDPos < MaxPos) {
            if (CDPos == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (CDPos == 1) {
              C_Dn = (Float_t)ampls[0];
              C_Dnn = 0;
            }
            if (CDPos > 1) {
              C_Dn = (Float_t)ampls[CDPos - 1];
              C_Dnn = (Float_t)ampls[CDPos - 2];
            }
            if (ampls.size() - LeftOfMaxPos > 1 && MaxPos + 2 < (int)(ampls.size()) - 1) {
              C_Mnn = (Float_t)ampls[MaxPos + 2];
            } else
              C_Mnn = 0;
          }
          if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
              C_Mnn <= coeff1 * coeffn * C_Mn + coeff2 * coeffnn * C_M + 2 * noise &&
              C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
            shapecdtn = true;
          }
        }
      } else if (exitCode)
        *exitCode = 6;
    }
  } else if (NofMax > 1 && exitCode)
    *exitCode = 1;  // more than one maximum
  if (ampls.size() == 1) {
    shapecdtn = true;
  }
  if (shapecdtn && exitCode)
    *exitCode = 0;

  return shapecdtn;
}

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

const float TkModGeomThickness[] = {1,
                                     0.029000,
                                     0.029000,
                                     0.047000,
                                     0.047000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.047000,
                                     0.047000,
                                     0.047000};
const float TkModGeomLength[] = {1,
                                  5.844250,
                                  5.844250,
                                  9.306700,
                                  9.306700,
                                  5.542900,
                                  4.408000,
                                  5.533000,
                                  4.258000,
                                  4.408000,
                                  5.533000,
                                  5.758000,
                                  7.363125,
                                  9.204400,
                                  10.235775};
const float TkModGeomWidthB[] = {1,
                                  3.072000,
                                  3.072000,
                                  4.684800,
                                  4.684800,
                                  4.579445,
                                  5.502407,
                                  3.158509,
                                  4.286362,
                                  5.502407,
                                  4.049915,
                                  3.561019,
                                  6.002559,
                                  5.235483,
                                  3.574395};
const float TkModGeomWidthT[] = {1,
                                  3.072000,
                                  3.072000,
                                  4.684800,
                                  4.684800,
                                  3.095721,
                                  4.322593,
                                  4.049915,
                                  3.146580,
                                  4.322593,
                                  3.158509,
                                  2.898798,
                                  4.824683,
                                  4.177638,
                                  4.398049};

bool isHitInsideTkModule(const LocalPoint hitPos, const DetId& detid, const SiStripCluster* cluster = nullptr) {
  if (detid.subdetId() < 3) {
    return true;
  }  //do nothing for pixel modules
  SiStripDetId SSdetId(detid);
  int moduleGeometry = SSdetId.moduleGeometry();

  //clean along the apv lines
  if (cluster &&
      (cluster->firstStrip() % 128 == 0 || (cluster->firstStrip() + cluster->amplitudes().size() % 128 == 127)))
    return false;

  float nx, ny;
  if (moduleGeometry <= 4) {
    ny = hitPos.y() / TkModGeomLength[moduleGeometry];
    nx = hitPos.x() / TkModGeomWidthT[moduleGeometry];
  } else {
    float offset =
        TkModGeomLength[moduleGeometry] * (TkModGeomWidthT[moduleGeometry] + TkModGeomWidthB[moduleGeometry]) /
        (TkModGeomWidthT[moduleGeometry] -
         TkModGeomWidthB[moduleGeometry]);  // check sign if GeomWidthT[moduleGeometry] < TkModGeomWidthB[moduleGeometry] !!!
    float tan_a = TkModGeomWidthT[moduleGeometry] / std::abs(offset + TkModGeomLength[moduleGeometry]);
    ny = hitPos.y() / TkModGeomLength[moduleGeometry];
    nx = hitPos.x() / (tan_a * std::abs(hitPos.y() + offset));
  }

  // "blacklists" for the gaps and edges
  switch (moduleGeometry) {
    case 0:
      return true;
    case 1:
      if (fabs(ny) > 0.96 || fabs(nx) > 0.98)
        return false;
      break;
    case 2:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 3:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) < 0.04)
        return false;
      break;
    case 4:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) < 0.04)
        return false;
      break;
    case 5:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98)
        return false;
      break;
    case 6:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 7:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.98)
        return false;
      break;
    case 8:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 9:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 10:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 11:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 12:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99 || (-0.17 < ny && ny < -0.07))
        return false;
      break;
    case 13:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99 || (-0.10 < ny && ny < -0.01))
        return false;
      break;
    case 14:
      if (fabs(ny) > 0.95 || fabs(nx) > 0.98 || (0.01 < ny && ny < 0.12))
        return false;
      break;
    default:
      std::cerr << "Unknown module geometry! Exiting!" << std::endl;
      exit(EXIT_FAILURE);
  }

  return true;
}

reco::DeDxData computedEdx(const float& track_eta,
                           const int& run_number,
                           string year,
                           const reco::DeDxHitInfo* dedxHits,
                           float* scaleFactors,
                           TH3* templateHisto = nullptr,
                           bool usePixel = false,
                           bool useStrip = true,
                           bool useClusterCleaning = true,
                           bool useTruncated = false,
                           bool mustBeInside = false,
                           size_t MaxStripNOM = 999,
                           bool correctFEDSat = false,
                           int crossTalkInvAlgo = 0,
                           float dropLowerDeDxValue = 0.0,
                           float* dEdxErr = nullptr,
                           bool useTemplateLayer = false,
                           bool skipPixelL1 = false,
                           int  skip_templates_ias = 0,
                           bool symmetricSmirnov = false,
                           bool useMorrisMethod = false) {

  if (!dedxHits)
    return reco::DeDxData(-1, -1, -1);

  std::vector<float> vect;
  std::vector<float> vectStrip;
  std::vector<float> vectPixel;

  // loop in order to have the number of saturated clusters in a track
  unsigned int nsatclust = 0;
  for (unsigned int t = 0; t < dedxHits->size(); t++) {
    DetId detid(dedxHits->detId(t));
    bool test_sat = false;
    if (detid.subdetId() < 3)
      continue;
    const SiStripCluster* cluster = dedxHits->stripCluster(t);
    std::vector<int> amplitudes = convert(cluster->amplitudes());
    for (unsigned int s = 0; s < amplitudes.size(); s++) {
      if (amplitudes[s] > 253)
        test_sat = true;
    }
    if (test_sat)
      nsatclust++;
  }
  float rsat = (float)nsatclust / (float)dedxHits->size();

  unsigned int NSat = 0;
  unsigned int SiStripNOM = 0;
  //float lowerStripDeDx=1000; UNUSED
  //int lowerStripDeDxIndex=-1; UNUSED
  for (unsigned int h = 0; h < dedxHits->size(); h++) {
    DetId detid(dedxHits->detId(h));
    if (!usePixel && detid.subdetId() < 3)
      continue;  // skip pixels
    if (!useStrip && detid.subdetId() >= 3)
      continue;  // skip strips

    if (mustBeInside &&
        !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId() >= 3 ? dedxHits->stripCluster(h) : nullptr))
      continue;
    if (detid.subdetId() >= 3 && ++SiStripNOM > MaxStripNOM)
      continue;  // skip remaining strips, but not pixel

    int ClusterCharge = dedxHits->charge(h);

    // Tav: This really should be done through the topology
    if (skipPixelL1 && detid.subdetId() == 1 && ((detid >> 20) & 0xF) == 1) //decoding mask for 2017-2018
      continue;

    if (detid.subdetId() >= 3) {  //for strip only
      
      if(fabs(track_eta) < 1.0 && !(detid.subdetId() == 3 || detid.subdetId() == 5)) continue; // eta < 1.0 -> only TIB+TOB hits
      if(fabs(track_eta) > 1.0 && track_eta < 1.7 && !(detid.subdetId() == 3 || detid.subdetId() == 4 || detid.subdetId() == 6)) continue; // 1.0 < eta < 1.7 -> only TIB+TID+TEC hits
      if(fabs(track_eta) > 1.7 && !(detid.subdetId() == 4 || detid.subdetId() == 6)) continue; // eta > 1.7 -> only TID+TEC hits

      SiStripDetId Sdetid(dedxHits->detId(h));
      const SiStripCluster* cluster = dedxHits->stripCluster(h);
      std::vector<int> amplitudes = convert(cluster->amplitudes());
      std::vector<int> amplitudesPrim = CrossTalkInv(amplitudes,0.10,0.04,true);

      // why is this hardcoded now?
      //if (useClusterCleaning && !clusterCleaning(amplitudes, crossTalkInvAlgo))
      if (useClusterCleaning && !clusterCleaning(amplitudesPrim, 1))
        continue;

      //////////////////////////////////////////////////////////////
      //
      //     0: no correction & no cross-talk inversion
      //     1: standard correction & use of cross-talk inversion
      //     2: correction from fits & no cross-talk inversion
      //     3: correction from fits & use of cross-talk inversion (no recorrection -- see bool=false)
      //
      //////////////////////////////////////////////////////////////

      //if (crossTalkInvAlgo>0) amplitudes = CrossTalkInv(amplitudes, 0.10, 0.04, true);

crossTalkInvAlgo=1;

      if (crossTalkInvAlgo == 1)
        //amplitudes = CrossTalkInv(amplitudes, 0.10, 0.04, true);
        amplitudes = SaturationCorrection(amplitudes,0.10,0.04,true,20,25);
      if (crossTalkInvAlgo == 2)
        amplitudes = Correction(amplitudes, Sdetid.moduleGeometry(), rsat, 25, 40, 0.6);
      if (crossTalkInvAlgo == 3)
        amplitudes =
            CrossTalkInv(Correction(amplitudes, Sdetid.moduleGeometry(), rsat, 25, 40, 0.6), 0.10, 0.04, false);

      float gain = 1.0;
      bool isSatCluster = false;
      ClusterCharge = 0;
      for (unsigned int s = 0; s < amplitudes.size(); s++) {
        int StripCharge = amplitudes[s];
        if (StripCharge < 254) {
          // TAV: it's kinda funny to divide with 1.0
          StripCharge = (int)(StripCharge / gain);
          if (StripCharge >= 1024) {
            StripCharge = 255;
          } else if (StripCharge >= 254) {
            StripCharge = 254;
          }
        }

        if (StripCharge >= 254) {
          isSatCluster = true;
        }
        if (StripCharge >= 255 && correctFEDSat) {
          StripCharge = 512;
        }
        ClusterCharge += StripCharge;
      }
      if (isSatCluster)
        NSat++;
    }

    float scaleFactor = scaleFactors[0];
    if (detid.subdetId() < 3){
      scaleFactor *= scaleFactors[1];  // add pixel scaling
      float pixelScaling = GetSFPixel(detid.subdetId(), detid, year, run_number);
      scaleFactor *= pixelScaling;
    }

    if (templateHisto) {  //save discriminator probability
      float ChargeOverPathlength =
          scaleFactor * ClusterCharge / (dedxHits->pathlength(h) * 10.0 * (detid.subdetId() < 3 ? 265 : 1));

      int moduleGeometry = 0;  // underflow for debug
      int layer = 0;
      if (detid.subdetId() < 3) {
        moduleGeometry = 15;
      }  // 15 == pixel
      else {
        SiStripDetId SSdetId(detid);
        moduleGeometry = SSdetId.moduleGeometry();
      }
      if (detid.subdetId() == 3) {
        layer = ((detid >> 14) & 0x7);
      }  //TIB
      if (detid.subdetId() == 4) {
        layer = ((detid >> 9) & 0x3) + 10;
      }  //TID
      if (detid.subdetId() == 5) {
        layer = ((detid >> 14) & 0x7) + 4;
      }  //TOB
      if (detid.subdetId() == 6) {
        layer = ((detid >> 5) & 0x7) + 13;
      }  //TEC

      //skip templates ias = 1 --> skip pixel, TIB, TID, 3 first TEC layers
      if (skip_templates_ias == 1 && (
                     detid.subdetId()<5 ||
                     layer == 14 ||
                     layer == 15 ||
                     layer == 16
                  )
        ){continue;}

      //skip templates ias = 2 --> pixel only, with pixL1 or not
      if (skip_templates_ias == 2 && (
                  detid.subdetId()>2 ||
                  //(skipPixelL1 && detid.subdetId() == 1 && ((detid >> 16) & 0xF) == 1) //decoding mask for 2016 !
                  (skipPixelL1 && detid.subdetId() == 1 && ((detid >> 20) & 0xF) == 1) //decoding mask for 2017-2018
                  )
        ){continue;}

      int BinX = templateHisto->GetXaxis()->FindBin(moduleGeometry);
      if (useTemplateLayer)
        BinX = templateHisto->GetXaxis()->FindBin(layer);
      int BinY = templateHisto->GetYaxis()->FindBin(dedxHits->pathlength(h) * 10.0);  //*10 because of cm-->mm
      int BinZ = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
      float Prob = templateHisto->GetBinContent(BinX, BinY, BinZ);
      //printf("%i %i %i  %f\n", BinX, BinY, BinZ, Prob);
      vect.push_back(Prob);  //save probability
      //printf("%i %i %i %i  %f\n", layer, BinX, BinY, BinZ, Prob);
    } else {
      float Norm = (detid.subdetId() < 3) ? 3.61e-06 : 3.61e-06 * 265;
      float ChargeOverPathlength = scaleFactor * Norm * ClusterCharge / dedxHits->pathlength(h);

      vect.push_back(ChargeOverPathlength);  //save charge
      if (detid.subdetId() < 3)
        vectPixel.push_back(ChargeOverPathlength);
        //vectPixel.push_back(ClusterCharge/dedxHits->pathlength(h));
      if (detid.subdetId() >= 3)
        vectStrip.push_back(ChargeOverPathlength);
        //vectStrip.push_back(ClusterCharge/dedxHits->pathlength(h));
      //           printf("%i - %f / %f = %f\n", h, scaleFactor*Norm*dedxHits->charge(h), dedxHits->pathlength(h), ChargeOverPathlength);
    }
  }

  if (dropLowerDeDxValue > 0) {
    std::vector<float> tmp(vect.size());
    std::copy(vect.begin(), vect.end(), tmp.begin());
    std::sort(tmp.begin(), tmp.end(), std::greater<float>());
    int nTrunc = tmp.size() * dropLowerDeDxValue;
    vect.clear();
    for (unsigned int t = 0; t + nTrunc < tmp.size(); t++) {
      vect.push_back(tmp[t]);
    }
  }

  float result;
  int size = vect.size();

  if (size > 0) {
    if (templateHisto) {
      if (useMorrisMethod) {
        // FiStrips discriminator
        float alpha = 1;
        for (int i = 0; i < size; i++) {
          alpha *= vect[i];
        }
        float logAlpha = log(alpha);
        float probQm = 0;
        for (int i = 0; i < size; i++) {
          probQm += ((pow(-logAlpha, i)) / (factorial(i)));
        }
        result = alpha * probQm;
      } else {
        //Ias discriminator
        result = 1.0 / (12 * size);
        std::sort(vect.begin(), vect.end(), std::less<float>());
        for (int i = 1; i <= size; i++) {
          if(!symmetricSmirnov) result += vect[i - 1] * pow(vect[i - 1] - ((2.0 * i - 1.0) / (2.0 * size)), 2); //Ias 
          else result += pow(vect[i - 1] - ((2.0 * i - 1.0) / (2.0 * size)), 2); //Is 
        }
        result *= (3.0 / size);
       }
    } else {  //dEdx estimator
      if (useTruncated) {
        //truncated40 estimator
        std::sort(vect.begin(), vect.end(), std::less<float>());
        result = 0;
        int nTrunc = size * 0.40;
        for (int i = 0; i + nTrunc < size; i++) {
          result += vect[i];
        }
        result /= (size - nTrunc);
      } else {
        //harmonic2 estimator (Ih)
        result = 0;
        float expo = -2;
        if (dEdxErr) *dEdxErr = 0;

        for (int i = 0; i < size; i++) {
          result += pow(vect[i], expo);
          if (dEdxErr) *dEdxErr += pow(vect[i], 2 * (expo - 1)) * pow(0.01, 2);
        }
        result = pow(result / size, 1. / expo);
        if (dEdxErr) *dEdxErr = result * result * result * sqrt(*dEdxErr) / size;
      }
      //           printf("Ih = %f\n------------------\n",result);
    }
  } else {
    result = -1;
  }
  return reco::DeDxData(result, NSat, size);
}
#endif  //FWCORE

// compute mass out of a momentum and dEdx value
float GetMass(float P, float I, float dEdxK, float dEdxC) {
  float& K = dEdxK;
  float& C = dEdxC;

  if (I - C < 0)
    return -1;
  return sqrt((I - C) / K) * P;
}

float GetMassErr(float P, float PErr, float dEdx, float dEdxErr, float M, float dEdxK, float dEdxC) {
  if (M < 0)
    return -1;
  float KErr = 0.2;
  //float CErr     = 0.4; UNUSED
  //float cErr     = 0.01; UNUSED
  float Criteria = dEdx - dEdxC;
  float Fac1 = P * P / (2 * M * dEdxK);
  float Fac2 = pow(2 * M * M * dEdxK / (P * P), 2);
  float MassErr = Fac1 * sqrt(Fac2 * pow(PErr / P, 2) + Criteria * Criteria * pow(KErr / dEdxK, 2) +
                               dEdxErr * dEdxErr + dEdxC * dEdxC);

  if (std::isnan(MassErr) || std::isinf(MassErr))
    MassErr = -1;

  return MassErr / M;
}

// pz compute Ick out of dEdx value
float GetIck(float I, float dEdxK, float dEdxC) {
  float& K = dEdxK;
  float& C = dEdxC;

  return (I - C) / K;
}

// compute mass out of a beta and momentum value
float GetMassFromBeta(float P, float beta) {
  float gamma = 1 / sqrt(1 - beta * beta);
  return P / (beta * gamma);
}

// compute mass out of a momentum and TOF value
float GetTOFMass(float P, float TOF) { return GetMassFromBeta(P, 1 / TOF); }

// estimate beta from a dEdx value, if dEdx is below the allowed threshold returns a negative beta value
float GetIBeta(float I, float dEdxK, float dEdxC) {
  float& K = dEdxK;
  float& C = dEdxC;

  float a = K / (I - C);
  float b2 = a / (a + 1);
  if (b2 < 0)
    return -1 * sqrt(b2);
  return sqrt(b2);
}

#ifdef FWCORE
float deltaROpositeTrack(const susybsm::HSCParticleCollection& hscpColl, const susybsm::HSCParticle& hscp) {
  reco::TrackRef track1 = hscp.trackRef();

  float maxDr = -0.1;
  for (unsigned int c = 0; c < hscpColl.size(); c++) {
    reco::TrackRef track2;
    if (!hscpColl[c].trackRef().isNull()) {
      track2 = hscpColl[c].trackRef();
    } else if (!hscpColl[c].muonRef().isNull() && hscpColl[c].muonRef()->combinedQuality().updatedSta) {
      track2 = hscpColl[c].muonRef()->standAloneMuon();
    } else {
      continue;
    }

    if (fabs(track1->pt() - track2->pt()) < 1 &&
        deltaR(track1->eta(), track1->phi(), track2->eta(), track2->phi()) < 0.1)
      continue;  //Skip same tracks
    //  float dR = deltaR(-1*track1->eta(), M_PI+track1->phi(), track2->eta(), track2->phi());
    TVector3 v1 = TVector3(track1->momentum().x(), track1->momentum().y(), track1->momentum().z());
    TVector3 v2 = TVector3(track2->momentum().x(), track2->momentum().y(), track2->momentum().z());
    float dR = v1.Angle(v2);
    if (dR > maxDr)
      maxDr = dR;
  }
  return maxDr;
}

//=============================================================
//
//     Method for Counting the number of muon stations used in track fit only counting DT and CSC stations.
//
//=============================================================
unsigned int muonStations(const reco::HitPattern& hitPattern) {
  unsigned int stations[4] = {0, 0, 0, 0};
  for (int i = 0; i < hitPattern.numberOfAllHits(reco::HitPattern::HitCategory::TRACK_HITS); i++) {
    uint32_t pattern = hitPattern.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i);
    if (pattern == 0)
      break;
    if (hitPattern.muonHitFilter(pattern) &&
        (int(hitPattern.getSubStructure(pattern)) == 1 || int(hitPattern.getSubStructure(pattern)) == 2) &&
        hitPattern.getHitType(pattern) == 0) {
      stations[hitPattern.getMuonStation(pattern) - 1] = 1;
    }
  }
  return stations[0] + stations[1] + stations[2] + stations[3];
}

//=============================================================
//
//     Trigger patterns
//
//=============================================================
bool passTriggerPatterns(edm::Handle<edm::TriggerResults> trigger,
                         const edm::TriggerNames triggerNames,
                         std::vector<std::string> patterns) {
  for (unsigned int i = 0; i < triggerNames.triggerNames().size(); i++) {
    TString name = triggerNames.triggerNames()[i];
    for (TString const& pattern : patterns) {
        if (pattern.Length()==0) {
          return false;
        }
      if (name.Contains(pattern) && trigger->accept(i))
        return true;
    }
  }
  return false;
}
bool passTriggerPatterns(edm::Handle<edm::TriggerResults> trigger,
                         const edm::TriggerNames triggerNames,
                         std::string pattern) {
  if (pattern.length()==0) {
    return false;
  }
  for (uint i = 0; i < triggerNames.triggerNames().size(); i++) {
    TString name = triggerNames.triggerNames()[i];
    if (name.Contains(pattern) && trigger->accept(i))
      return true;
  }
  return false;
}

//=============================================================
//
//     selection of gen HSCP
//
//=============================================================
bool isGoodGenHSCP(const reco::GenParticle& gen, bool onlyCharged) {
  if (gen.status() != 1)
    return false;
  if (gen.pt() < 5)
    return false;
  int AbsPdg = abs(gen.pdgId());
  if (AbsPdg < 1000000 && AbsPdg != 17)
    return false;
  if (onlyCharged && (AbsPdg == 1000993 || AbsPdg == 1009313 || AbsPdg == 1009113 || AbsPdg == 1009223 ||
                      AbsPdg == 1009333 || AbsPdg == 1092114 || AbsPdg == 1093214 || AbsPdg == 1093324))
    return false;  //Skip neutral gluino RHadrons
  if (onlyCharged && (AbsPdg == 1000622 || AbsPdg == 1000642 || AbsPdg == 1006113 || AbsPdg == 1006311 ||
                      AbsPdg == 1006313 || AbsPdg == 1006333))
    return false;  //skip neutral stop RHadrons
  if (onlyCharged && AbsPdg == 1000022)
    return false;  //skip neutralino
  return true;
}

//=============================================================
//
//     Count Number of Charged HSCP
//
//=============================================================
int HowManyChargedHSCP(const std::vector<reco::GenParticle>& genColl) {
  int toReturn = 0;
  for (unsigned int g = 0; g < genColl.size(); g++) {
    if (isGoodGenHSCP(genColl[g], true))
      toReturn++;
  }
  return toReturn;
}

//=============================================================
//
//     Compute Decay length of generated HSCP
//
//=============================================================
void GetGenHSCPDecayLength(const std::vector<reco::GenParticle>& genColl,
                           float& length1,
                           float& length2,
                           bool onlyCharged) {
  length1 = -1;
  length2 = -1;
  for (auto const& mcParticle : genColl) {
    if (!isGoodGenHSCP(mcParticle, onlyCharged))
      continue;

    //code below copied from disapearing track source code.  FIXME we should make sure that the logic is ok.
    //https://raw.githubusercontent.com/OSU-CMS/DisappTrks/master/SignalMC/plugins/DecayAnalyzer.cc
    const reco::Candidate* mother = &mcParticle;
    const reco::Candidate* daughter = mother;  // FIXME

    // Descend the decay chain until no daughters have the same PDG ID as mcParticle.
    while (true) {
      bool foundDauSamePdgId = false;
      for (uint i = 0; i < daughter->numberOfDaughters(); i++) {
        if (daughter->daughter(i)->pdgId() == mcParticle.pdgId()) {
          foundDauSamePdgId = true;
          mother = daughter;
          daughter = daughter->daughter(i);
          break;
        }
      }
      if (!foundDauSamePdgId)
        break;
    }

    // Now daughter has no daughters with the same PDG ID as mcParticle.
    // Next choose the daughter with the outermost production vertex, in case there are multiple vertices
    // (e.g., an electron delta ray can produce a vertex before the decay vertex)
    float radiusLastVtx = -99;
    int idxDauLastVtx = -99;
    for (uint i = 0; i < daughter->numberOfDaughters(); i++) {
      float radius = daughter->daughter(i)->vertex().R();
      if (radius > radiusLastVtx) {
        radiusLastVtx = radius;
        idxDauLastVtx = i;
      }
    }
    if (idxDauLastVtx < 0)
      continue;

    mother = daughter;
    daughter = daughter->daughter(idxDauLastVtx);

    TVector3 source(mcParticle.vx(), mcParticle.vy(), mcParticle.vz());
    TVector3 decay(daughter->vx(), daughter->vy(), daughter->vz());
    float ctau = (decay - source).Mag() / (mcParticle.p4().Beta() * mcParticle.p4().Gamma());

    if (length1 < 0) {
      length1 = ctau;
    } else if (length2 < 0) {
      length2 = ctau;
      return;
    }
  }
  if (length1 < 0) {
    length1 = 9999;
  }
  if (length2 < 0) {
    length2 = 9999;
  }
}

// Same as above but not forcing it to be HSCP
void GetGenBcgDecayLength(const std::vector<reco::GenParticle>& genColl,
                           float& length1,
                           float& length2,
                           bool onlyCharged) {
  length1 = -1;
  length2 = -1;
  for (auto const& mcParticle : genColl) {

    const reco::Candidate* mother = &mcParticle;
    const reco::Candidate* daughter = mother;  // FIXME

    // Descend the decay chain until no daughters have the same PDG ID as mcParticle.
    while (true) {
      bool foundDauSamePdgId = false;
      for (uint i = 0; i < daughter->numberOfDaughters(); i++) {
        if (daughter->daughter(i)->pdgId() == mcParticle.pdgId()) {
          foundDauSamePdgId = true;
          mother = daughter;
          daughter = daughter->daughter(i);
          break;
        }
      }
      if (!foundDauSamePdgId)
        break;
    }

    // Now daughter has no daughters with the same PDG ID as mcParticle.
    // Next choose the daughter with the outermost production vertex, in case there are multiple vertices
    // (e.g., an electron delta ray can produce a vertex before the decay vertex)
    float radiusLastVtx = -99;
    int idxDauLastVtx = -99;
    for (uint i = 0; i < daughter->numberOfDaughters(); i++) {
      float radius = daughter->daughter(i)->vertex().R();
      if (radius > radiusLastVtx) {
        radiusLastVtx = radius;
        idxDauLastVtx = i;
      }
    }
    if (idxDauLastVtx < 0)
      continue;

    mother = daughter;
    daughter = daughter->daughter(idxDauLastVtx);

    TVector3 source(mcParticle.vx(), mcParticle.vy(), mcParticle.vz());
    TVector3 decay(daughter->vx(), daughter->vy(), daughter->vz());
    float ctau = (decay - source).Mag() / (mcParticle.p4().Beta() * mcParticle.p4().Gamma());

    if (length1 < 0) {
      length1 = ctau;
    } else if (length2 < 0) {
      length2 = ctau;
      return;
    }
  }
  if (length1 < 0) {
    length1 = 9999;
  }
  if (length2 < 0) {
    length2 = 9999;
  }
}

//=============================================================
//
//     Returns the generated beta of the two firsts HSCP in the events
//
//=============================================================
void GetGenHSCPBeta(const std::vector<reco::GenParticle>& genColl, float& beta1, float& beta2, bool onlyCharged) {
  beta1 = -1;
  beta2 = -1;
  for (auto const& mcParticle : genColl) {
    if (!isGoodGenHSCP(mcParticle, onlyCharged))
      continue;
    if (beta1 < 0) {
      beta1 = mcParticle.p() / mcParticle.energy();
    } else if (beta2 < 0) {
      beta2 = mcParticle.p() / mcParticle.energy();
      return;
    }
  }
}

//============================================================
//
//      Return tk-based isolation with a PV constrain
//
//============================================================
float TkBasedIsolationWithPVConstrain(const std::vector<reco::Track>& tkTracks,
                                       const reco::Vertex& PV,
                                       reco::TrackRef TkRef,
                                       const float& dRmax,
                                       const float& dxyMax,
                                       const float& dzMax,
                                       const float& trackSelection = 0) {
  if (TkRef->pt() < 10)
    return -1;
  float SumPt = 0.;
  for (unsigned int i = 0; i < tkTracks.size(); i++) {
    reco::Track itTrack = tkTracks[i];
    if (itTrack.pt() < trackSelection)
      continue;
    if (fabs(itTrack.pt() - TkRef->pt()) < 0.1 && fabs(itTrack.eta() - TkRef->eta()) < 0.05)
      continue;
    float dR = deltaR(itTrack.eta(), itTrack.phi(), TkRef->eta(), TkRef->phi());
    if (dR > dRmax)
      continue;
    if (itTrack.dxy(PV.position()) > dxyMax || itTrack.dz(PV.position()) > dzMax)
      continue;
    SumPt += itTrack.pt();
  }
  return SumPt;
}

//============================================================
//
//      Return cut on sigma_pT / pT for as a function of
//      3 eta bins and for a given sigma_pT quantile
//
//============================================================

TF1 f_eta_0p8_q95("f_eta_0p8_q95","0.006725+0.000122*x",0,10000);
TF1 f_0p8_eta_1p7_q95("f_0p8_eta_1p7_q95","0.012919+0.000168*x",0,10000);
TF1 f_1p7_eta_2p1_q95("f_1p7_eta_2p1_q95","0.015408+0.000320*x",0,10000);

TF1 f_eta_0p8_q99("f_eta_0p8_q99","0.006428+0.000174*x",0,10000);
TF1 f_0p8_eta_1p7_q99("f_0p8_eta_1p7_q99","0.012825+0.000222*x",0,10000);
TF1 f_1p7_eta_2p1_q99("f_1p7_eta_2p1_q99","0.015597+0.000403*x",0,10000);

TF1 f_eta_0p8_q999("f_eta_0p8_q999","0.006260+0.000258*x",0,10000);
TF1 f_0p8_eta_1p7_q999("f_0p8_eta_1p7_q999","0.018252+0.000372*x",0,10000);
TF1 f_1p7_eta_2p1_q999("f_1p7_eta_2p1_q999","-0.008672+0.001016*x",0,10000);


float pTerr_over_pT_etaBin(const float& pt, const float& eta, float q=99){
    if(q==999){
        if(abs(eta)<=0.8) return std::max(f_eta_0p8_q999(pt),0.01);
        else if(abs(eta)>0.8 && abs(eta)<=1.7) return std::max(f_0p8_eta_1p7_q999(pt),0.01);
        else if(abs(eta)>1.7 && abs(eta)<=2.1) return std::max(f_1p7_eta_2p1_q999(pt),0.01);
    }else if(q==99){
        if(abs(eta)<=0.8) return std::max(f_eta_0p8_q99(pt),0.01);
        else if(abs(eta)>0.8 && abs(eta)<=1.7) return std::max(f_0p8_eta_1p7_q99(pt),0.01);
        else if(abs(eta)>1.7 && abs(eta)<=2.1) return std::max(f_1p7_eta_2p1_q99(pt),0.01);
    }else if(q==95){
        if(abs(eta)<=0.8) return std::max(f_eta_0p8_q95(pt),0.01);
        else if(abs(eta)>0.8 && abs(eta)<=1.7) return std::max(f_0p8_eta_1p7_q95(pt),0.01);
        else if(abs(eta)>1.7 && abs(eta)<=2.1) return std::max(f_1p7_eta_2p1_q95(pt),0.01);
    }
    return -1.;
}


#endif  //FWCORE

#endif
