#ifndef DEdxHIPEmulator_h
#define DEdxHIPEmulator_h

struct HitDeDx {
  double dedx;
  double dx;
  bool isSat;
  bool passClusterCleaning;
  bool isInside;
  unsigned char subDet;
};

class DEdxHIPEmulator {
private:
  TH1D* ratePdfPixel;
  TH1D* ratePdfStrip;

  double eventRatePixel;
  double eventRateStrip;

  bool is2016;

public:
  typedef std::vector<HitDeDx> HitDeDxCollection;
  void setEventRate(double ratePixel = -1, double rateStrip = -1) {
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

    eventRatePixel = std::max(eventRatePixel, 0.0);
    eventRateStrip = std::max(eventRateStrip, 0.0);
  }

  void setPeriodHIPRate(bool is2016_, const char* ratePdfPixelName = NULL, const char* ratePdfStripName = NULL) {
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

  DEdxHIPEmulator(bool is2016_ = false, const char* ratePdfPixelName = NULL, const char* ratePdfStripName = NULL) {
    setPeriodHIPRate(is2016_, ratePdfPixelName, ratePdfStripName);
  }
  ~DEdxHIPEmulator() {}

  double getEventRatePixel() { return eventRatePixel; }
  double getEventRateStrip() { return eventRateStrip; }

  double fakeHIP(unsigned int subDet, double dedx) {
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

#endif
