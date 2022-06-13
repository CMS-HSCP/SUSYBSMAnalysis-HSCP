{

//  .L run2study.C

root -l
.L run2study.C
 TChain *f_2016b = new TChain("stage/ttree");
 f_2016b->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunBv2/200909_151602/0000/nt_aod_ul2016_*.root");
 f_2016b->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunBv2/200909_151602/0001/nt_aod_ul2016_*.root");
 run2study t0(f_2016b);
 t0.Loop(2016,"B");
.q

root -l
.L run2study.C
 TChain *f_2016c = new TChain("stage/ttree");
 f_2016c->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunC/200909_151617/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016c);
 t1.Loop(2016,"C");                                                                                                                                       
.q                                                                                                                                                        


root -l
.L run2study.C
 TChain *f_2016d = new TChain("stage/ttree");
 f_2016d->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunD/200909_151630/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016d);
 t1.Loop(2016,"D");
.q

root -l
.L run2study.C
 TChain *f_2016e = new TChain("stage/ttree");
 f_2016e->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunE/200909_151658/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016e);
 t1.Loop(2016,"E");
.q

root -l
.L run2study.C
 TChain *f_2016f = new TChain("stage/ttree");
 f_2016f->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_HIPM_RunF/200909_151751/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016f);
 t1.Loop(2016,"F");
.q


root -l
.L run2study.C
 TChain *f_2016f = new TChain("stage/ttree");
 f_2016f->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_RunF/200909_151805/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016f);
 t1.Loop(2016,"F2");
.q

root -l
.L run2study.C
 TChain *f_2016g = new TChain("stage/ttree");
 f_2016g->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_RunG/200909_151821/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016g);
 t1.Loop(2016,"G");
.q

root -l
.L run2study.C
 TChain *f_2016h = new TChain("stage/ttree");
 f_2016h->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodSept2020_CMSSW_10_6_12/SingleMuon/UL2016_RunH/200909_151833/0000/nt_aod_ul2016_*.root");
 run2study t1(f_2016h);
 t1.Loop(2016,"H");
.q




root -l
.L run2study.C
 TChain *f_2017b = new TChain("stage/ttree");
 f_2017b->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunB/200831_091142/0000/nt_aod_ul2017*.root");
 run2study t1(f_2017b);
 t1.Loop(2017,"B");
.q

root -l
.L run2study.C
 TChain *f_2017c = new TChain("stage/ttree");
 f_2017c->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunC/200831_091638/0000/nt_aod_ul2017*.root");
 f_2017c->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunC/200831_091638/0001/nt_aod_ul2017*.root");
 run2study t1(f_2017c);
 t1.Loop(2017,"C");
.q
root -l
.L run2study.C
 TChain *f_2017d = new TChain("stage/ttree");
 f_2017d->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunD/200831_091735/0000/nt_aod_ul2017*.root");
 run2study t1(f_2017d);
 t1.Loop(2017,"D");
.q
root -l
.L run2study.C
 TChain *f_2017e = new TChain("stage/ttree");
 f_2017e->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunE/200831_091816/0000/nt_aod_ul2017*.root");
 f_2017e->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunE/200831_091816/0001/nt_aod_ul2017*.root");
 run2study t1(f_2017e);
 t1.Loop(2017,"E");
.q
root -l
.L run2study.C
 TChain *f_2017f = new TChain("stage/ttree");
 f_2017f->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunF/200831_091901/0000/nt_aod_ul2017*.root");
 f_2017f->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodAug2020_CMSSW_10_6_12/SingleMuon/UL2017_RunF/200831_091901/0001/nt_aod_ul2017*.root");
 run2study t1(f_2017f);
 t1.Loop(2017,"F");
.q

/*

root -l
.L run2study.C
 TChain *f_2018a = new TChain("stage/ttree");
 f_2018a->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2018A/191024_133317/0000/nt_data_aod_*.root");
 run2study t1(f_2018a);
 t1.Loop(2018,"A");
.q

root -l
.L run2study.C
 TChain *f_2018b = new TChain("stage/ttree");
 f_2018b->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2018B/191024_133631/0000/nt_data_aod_*.root");
 run2study t1(f_2018b);
 t1.Loop(2018,"B");
.q

root -l
.L run2study.C
 TChain *f_2018c = new TChain("stage/ttree");
 f_2018c->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2018C/191024_133646/0000/nt_data_aod_*.root");
 run2study t1(f_2018c);
 t1.Loop(2018,"C");
.q

root -l
.L run2study.C
 TChain *f_2018d = new TChain("stage/ttree");
 f_2018d->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2018D/191024_133704/0000/nt_data_aod_*.root");
 run2study t1(f_2018d);
 t1.Loop(2018,"D");
.q
*/


root -l
.L run2study.C
 TChain *f_gluino = new TChain("stage/ttree");
 f_gluino->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/HSCPgluino_M_1400_TuneCP5_13TeV_pythia8/MC_Gl1400NoPuKazana/211008_132040/0000/nt_mc_aod_*.root");
  run2study t1(f_gluino);
  t1.Loop(1400,"MC_1400",false);
.q

root -l
.L run2study.C
 TChain *f_gluino = new TChain("stage/ttree");
 f_gluino->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/HSCPgluino_M_1000_TuneCP5_13TeV_pythia8/MC_Gl1000PuKazana/211008_133144/0000/nt_mc_aod_*.root");
  run2study t1(f_gluino);
  t1.Loop(1000,"MC_1000",false);
.q

//root -l
//.L run2study.C
// TChain *f_w17 = new TChain("stage/ttree");
// f_w17->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MC_W17/211008_130800/0000/nt_mc_aod_*.root");
//  run2study t1(f_w17);
//  t1.Loop(1000,"MC_w17",false);
//.q

root -l
.L run2study2.C
 TChain *f_w17 = new TChain("stage/ttree");
 f_w17->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/Prod2/211015_100302/0000/n*.root");
  run2study t1(f_w17);
  t1.Loop(1000,"MC_w17_0",false);
.q

root -l
.L run2study2.C
 TChain *f_w17_1 = new TChain("stage/ttree");
 f_w17_1->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/Prod2/211015_100302/0001/n*.root");
  run2study t1(f_w17_1);
  t1.Loop(1000,"MC_w17_1",false);
.q

root -l
.L run2study2.C
 TChain *f_w17_2 = new TChain("stage/ttree");
 f_w17_2->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2021_CMSSW_10_6_27/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/Prod2/211015_100302/0002/n*.root");
  run2study t2(f_w17_2);
  t2.Loop(1000,"MC_w17_2",false);
.q




