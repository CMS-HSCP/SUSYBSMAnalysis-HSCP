root -l -b << EOF
   TString makeshared(gSystem->GetMakeSharedLib());
   makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated -Wno-unused-local-typedefs -Wno-attributes ");
   makeshared.ReplaceAll("-Woverloaded-virtual ", " ");
   makeshared.ReplaceAll("-Wshadow ", " -std=c++0x -D__USE_XOPEN2K8 ");
   cout << "Compilling with the following arguments: " << makeshared << endl;
   gSystem->SetMakeSharedLib(makeshared);
   gSystem->SetIncludePath("-I$ROOFITSYS/include");
   //.x Analysis_Step4_LimitComputation.C++("Final7TeV", "", "");
   //.x Analysis_Step4_LimitComputation.C++("Final8TeV", "", "");
   //.x Analysis_Step4_LimitComputation.C++("FinalCOMBRun1", "", "");
   //.x Analysis_Step4_LimitComputation.C++("Final13TeV15" , "", "");
   .x Analysis_Step4_LimitComputation.C++("Final13TeV16" , "", "");
//mk   .x Analysis_Step4_LimitComputation.C++("Final13TeV16G" , "", "");
//mk   .x Analysis_Step4_LimitComputation.C++("FinalCOMB2016", "", "");
EOF

