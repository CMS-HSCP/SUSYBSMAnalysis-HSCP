import ROOT

# tdrGrid: Turns the grid lines on (true) or off (false)

def setTDRStyle():
    tdrStyle = ROOT.TStyle("tdrStyle","Style for P-TDR");

    tdrStyle.SetPalette(1)

    #    For Legend

    tdrStyle.SetLegendBorderSize(0);
    
    # For the canvas:
    tdrStyle.SetCanvasBorderMode(0);
    tdrStyle.SetCanvasColor(ROOT.kWhite);
    tdrStyle.SetCanvasDefH(600); #Height of canvas
    tdrStyle.SetCanvasDefW(600); #Width of canvas
    tdrStyle.SetCanvasDefX(0);   #Position on screen
    tdrStyle.SetCanvasDefY(0);

    
    # For the Pad:
    tdrStyle.SetPadBorderMode(0);
    # tdrStyle.SetPadBorderSize(Width_t size = 1);
    tdrStyle.SetPadColor(ROOT.kWhite);
    tdrStyle.SetPadGridX(False);
    tdrStyle.SetPadGridY(False);
    #tdrStyle.SetPadGridX(True);
    #tdrStyle.SetPadGridY(True);
    tdrStyle.SetGridColor(0);
    tdrStyle.SetGridStyle(3);
    tdrStyle.SetGridWidth(1);

    #    For the frame:
    tdrStyle.SetFrameBorderMode(0);
    tdrStyle.SetFrameBorderSize(1);
    tdrStyle.SetFrameFillColor(0);
    tdrStyle.SetFrameFillStyle(0);
    tdrStyle.SetFrameLineColor(1);
    tdrStyle.SetFrameLineStyle(1);
    tdrStyle.SetFrameLineWidth(2);

    # For the histo:
    # tdrStyle->SetHistFillColor(1);
    # tdrStyle->SetHistFillStyle(0);
    tdrStyle.SetHistLineColor(1);
    tdrStyle.SetHistLineStyle(0);
    tdrStyle.SetHistLineWidth(1);
    # tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
    # tdrStyle->SetNumberContours(Int_t number = 20);

    #tdrStyle.SetEndErrorSize(1);
    #GHM  tdrStyle->SetEndErrorSize(2);#
    #  tdrStyle->SetErrorMarker(20);
    #tdrStyle.SetErrorX(1.);
    
    tdrStyle.SetMarkerStyle(20);

    #For the fit/function:
    tdrStyle.SetOptFit(0);
    tdrStyle.SetFitFormat("5.4g");
    tdrStyle.SetFuncColor(2);
    tdrStyle.SetFuncStyle(1);
    tdrStyle.SetFuncWidth(2);

    #For the date:
    tdrStyle.SetOptDate(0);
    # tdrStyle.SetDateX(Float_t x = 0.01);
    # tdrStyle.SetDateY(Float_t y = 0.01);
    
    #For the statistics box:
    tdrStyle.SetOptFile(0);
    tdrStyle.SetOptStat(0)
    #tdrStyle.SetOptStat('mr'); # To display the mean and RMS:   SetOptStat("mr");
    tdrStyle.SetStatColor(ROOT.kWhite);
    tdrStyle.SetStatFont(42);
    tdrStyle.SetStatFontSize(0.025);
    tdrStyle.SetStatTextColor(1);
    tdrStyle.SetStatFormat("6.4g");
    tdrStyle.SetStatBorderSize(1);
    tdrStyle.SetStatH(0.1);
    tdrStyle.SetStatW(0.15);
    # tdrStyle->SetStatStyle(Style_t style = 1001);
    
    # tdrStyle->SetStatX(Float_t x = 0);
    # tdrStyle->SetStatY(Float_t y = 0);
    # Margins:
    tdrStyle.SetPadTopMargin(0.065);
    tdrStyle.SetPadBottomMargin(0.12);
    tdrStyle.SetPadLeftMargin(0.15);
    tdrStyle.SetPadRightMargin(0.05);
    
    # For the Global title:
    
    tdrStyle.SetOptTitle(1);
    tdrStyle.SetTitleFont(42);
    tdrStyle.SetTitleColor(1);
    tdrStyle.SetTitleTextColor(1);
    tdrStyle.SetTitleFillColor(10);
    tdrStyle.SetTitleFontSize(0.0525);
    tdrStyle.SetTitleH(0); # Set the height of the title box
    tdrStyle.SetTitleW(0); # Set the width of the title box
    tdrStyle.SetTitleX(0.2); # Set the position of the title box
    tdrStyle.SetTitleY(1.0055); # Set the position of the title box
    tdrStyle.SetTitleStyle(1001);
    tdrStyle.SetTitleBorderSize(0);
    #tdrStyle.SetTitleAlign(23);
    tdrStyle.SetTitleAlign(13)
    
    # For the axis titles:

    tdrStyle.SetTitleColor(1, "XYZ");
    tdrStyle.SetTitleFont(42, "XY");
    tdrStyle.SetTitleFont(42, "Z");
    tdrStyle.SetTitleSize(0.05, "XY");
    tdrStyle.SetTitleSize(0.035, "Z");
    #tdrStyle.SetTitleSize(0.06, "XYZ");
    # tdrStyle.SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
    # tdrStyle.SetTitleYSize(Float_t size = 0.02);
    tdrStyle.SetTitleXOffset(1.0);
    #  tdrStyle.SetTitleYOffset(1.25);
    tdrStyle.SetTitleYOffset(1.0);
    # tdrStyle.SetTitleOffset(1.1, "Y"); # Another way to set the Offset

    # For the axis labels:

    tdrStyle.SetLabelColor(1, "XYZ");
    tdrStyle.SetLabelFont(42, "XYZ");
    tdrStyle.SetLabelOffset(0.005, "XYZ");
    #tdrStyle.SetLabelSize(0.04, "XYZ");
    tdrStyle.SetLabelSize(0.03, "XYZ");

    # For the axis:

    tdrStyle.SetAxisColor(1, "XYZ");
    tdrStyle.SetStripDecimals(ROOT.kTRUE);
    tdrStyle.SetTickLength(0.01, "X");
    tdrStyle.SetTickLength(0.01, "Y");
    
    tdrStyle.SetNdivisions(510, "XYZ");
    tdrStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
    tdrStyle.SetPadTickY(1);

    # Change for log plots:
    tdrStyle.SetOptLogx(0);
    tdrStyle.SetOptLogy(0);
    tdrStyle.SetOptLogz(0);
    
    # Postscript options:
    # tdrStyle.SetPaperSize(20.,20.);
    # tdrStyle.SetLineScalePS(Float_t scale = 3);
    # tdrStyle.SetLineStyleString(Int_t i, const char* text);
    # tdrStyle.SetHeaderPS(const char* header);
    # tdrStyle.SetTitlePS(const char* pstitle);
    
    # tdrStyle.SetBarOffset(Float_t baroff = 0.5);
    # tdrStyle.SetBarWidth(Float_t barwidth = 0.5);
    # tdrStyle.SetPaintTextFormat(const char* format = "g");
    # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    # tdrStyle.SetTimeOffset(Double_t toffset);
    # tdrStyle.SetHistMinimumZero(kTRUE);

    tdrStyle.cd();
    
    ROOT.gROOT.ForceStyle() 
  


#void fixOverlay() {
#  gPad.RedrawAxis();
#}

#endif
