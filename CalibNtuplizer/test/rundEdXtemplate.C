{
 .L dEdXtemplate.C

 TChain *superTree = new TChain("stage/ttree");
 superTree->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP/minbias/test_minbias_aod1.root");
 superTree->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP/minbias/test_minbias_aod2.root");
 superTree->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP/minbias/test_minbias_aod3.root");
 superTree->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP/minbias/test_minbias_aod4.root");




 dEdXtemplate t(superTree);
 t.Loop("minbias",1);

  /*
  superTree->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP/ttbar/test_ttbar_aod2.root");
 dEdXtemplate t(superTree);
 t.Loop("ttbar",1);
  */
}
