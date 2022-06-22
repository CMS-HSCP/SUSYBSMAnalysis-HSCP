import argparse
import ROOT
from array import array
import ROOT, sys, os, time, re
import tdrstyle

def draw(file_map, h_name, h_name2, layer,  outfile, fileName):

#    Gain_40E14Ne_T253K = ROOT.TGraphErrors("./SignalBackgroundEff_HSCP.txt","%lg %lg")
    Gain_40E14Ne_T253K = ROOT.TGraphErrors(fileName,"%lg %lg")
    
    Gain_40E14Ne_T253K.SetMarkerStyle(20)
    Gain_40E14Ne_T253K.SetMarkerColor(1)
    Gain_40E14Ne_T253K.SetLineColor(1)
    Gain_40E14Ne_T253K.SetLineStyle(1)

    # ------------------------------------------------------
        
    mgGainResult =  ROOT.TMultiGraph()
    mgGainResult.Add(Gain_40E14Ne_T253K,"p");
    mgGainResult.SetTitle(";HSCP mass [GeV];Signal efficiency [%]")
    mgGainResult.GetYaxis().SetTitleSize(.05)
    mgGainResult.GetYaxis().SetTitleOffset(1.5)
    mgGainResult.GetYaxis().SetLabelSize(.04)
    mgGainResult.GetXaxis().SetTitleSize(.05)
    mgGainResult.GetXaxis().SetLabelSize(.04)
    mgGainResult.GetYaxis().SetRangeUser(0.,100.)
    
    legGainResult = ROOT.TLegend(.5,.30,.8,.545)
    legGainResult.SetTextFont(42)
    legGainResult.SetTextSize(0.035)
    legGainResult.AddEntry(Gain_40E14Ne_T253K,"R-hadron with gluino","LP")
        
    tex2 = ROOT.TLatex(0.18,0.96,"CMS");
    #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
    tex2.SetNDC();
    tex2.SetTextFont(61);
    tex2.SetTextSize(0.0375);
    tex2.SetLineWidth(2);

    tex3 = ROOT.TLatex(0.27,0.96,"Simulation")
    #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
    tex3.SetNDC();
    tex3.SetTextFont(52);
    tex3.SetTextSize(0.0285);
    tex3.SetLineWidth(2);
    
##    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = #frac{1}{1-p0#upointexp#left(-#frac{p1}{E}#right)}") #Morris
#    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = exp[A#upointexp(-b/E)]") #Morris2.0
##    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = #frac{1}{1-p0#upointexp#left(-#left(#frac{p1}{E}#right)^{p2}#right)}")
##    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = #frac{1}{1-(p0*x/p1)^{p2}}") #Muller-Moll
#    texFunction.SetNDC();
#    texFunction.SetTextFont(42)
#    texFunction.SetTextSize(0.035)
#    texFunction.SetTextColor(2)
    
    #    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = #frac{1}{1-p0#upointexp#left(-#left(#frac{p1}{E}#right)^{p2}#right)}")
#    texFunction = ROOT.TLatex(0.27,0.76,"g(E) = #frac{1}{1-(p0*x/p1)^{p2}}") #Muller-Moll
    
    
#    texFunctionPar0 = ROOT.TLatex(0.27,0.7,"#chi^{2}/ndf=0.3096/5")
#    texFunctionPar0.SetNDC();
#    texFunctionPar0.SetTextFont(42)
#    texFunctionPar0.SetTextSize(0.035)
##    texFunctionPar.SetTextColor(2)
#
#    texFunctionPar1 = ROOT.TLatex(0.27,0.65,"A=1.526x10^{9}")
#    texFunctionPar1.SetNDC();
#    texFunctionPar1.SetTextFont(42)
#    texFunctionPar1.SetTextSize(0.035)
#    texFunctionPar.SetTextColor(2)
#
#    texFunctionPar2 = ROOT.TLatex(0.27,0.6,"b=2.898x10^{7} V/cm")
#    texFunctionPar2.SetNDC();
#    texFunctionPar2.SetTextFont(42)
#    texFunctionPar2.SetTextSize(0.035)
    
    cGainResult = ROOT.TCanvas('cGainResult', 'cGainResult',800,800)
    cGainResult.SetRightMargin(0.08)
    mgGainResult.Draw('AP')
#    mgGainResult.Fit("expFit")
#    expFit.Draw("SAME")
    legGainResult.Draw("SAME")
    tex2.Draw("SAME");
    tex3.Draw("SAME");
    
#    texFunction.Draw("SAME")
#    texFunctionPar0.Draw("SAME")
#    texFunctionPar1.Draw("SAME")
#    texFunctionPar2.Draw("SAME")
    
    cGainResult.SaveAs("SignalEffHSCPGluino_"+fileName[:-4]+".png")
    
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()
#    ROOT.gStyle.SetOptFit(111)

 #-------------------------------------------------------------------------
 #  PARSE COMMAND LINE ARGUMENT
 #-------------------------------------------------------------------------
    
    parser = argparse.ArgumentParser(description='Charge MPV extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    parser.add_argument('-n', action="store", dest = 'h_name', default = 'analysis/h701_n1')
    parser.add_argument('-n2', action="store", dest = 'h_name2', default = 'analysis/h701_n1c')
    parser.add_argument('-lay', action="store", dest = 'layer', default = 'Layer 1')
    parser.add_argument('-o', action="store", dest = 'outfile', default = 'tmp.pdf')
    parser.add_argument('-f', action="store", dest = 'fileName', default = 'tmp.pdf')

    p = parser.parse_args()

    draw(p.inlist, p.h_name, p.h_name2, p.layer, p.outfile, p.fileName)
