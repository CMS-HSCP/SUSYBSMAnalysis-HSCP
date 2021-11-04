#ifndef ModuleGeom_H
#define ModuleGeom_H

#include <exception>
#include <vector>
#include <unordered_map>
#include "TVector3.h"
#include <string>
#include "TChain.h"


class ModuleGeom{
  //A singleton containing a different class to provide the geometry info would be better than hiding global vars behind static
  //Even better would be using the DB to get the current geometry!
   public:
      float TP; TVector3 pos; TVector3 width; TVector3 length; TVector3 thick;
      static std::unordered_map<unsigned int, ModuleGeom*> static_geomMap;
   public :
        ModuleGeom(){}
        ~ModuleGeom(){}

        TVector3 toGlobal(TVector3 local ){ return (pos + local.x()*width.Unit() + local.y()*length.Unit() + local.z()*thick.Unit());}
        TVector3 toLocal (TVector3 global){ TVector3 o = global-pos;  return TVector3((o*width.Unit()), (o*length.Unit()), (o*thick.Unit()));}

        bool propagateParametersOnPlane(TVector3& pos, TVector3& momentum, TVector3& localPosOnPlane, bool debug=false){
           if(debug){printf("PPOP pos=(%f,%f,%f) width=(%f,%f,%f) length=(%f,%f,%f) thick=(%f,%f,%f)\n", pos.x(), pos.y(), pos.z(), width.x(), width.y(), width.z(), length.x(), length.y(), length.z(), thick.x(), thick.y(), thick.z()); }

           TVector3 x = toLocal(pos);
           TVector3 p = TVector3((momentum*width.Unit()), (momentum*length.Unit()), (momentum*thick.Unit())).Unit();

           if(debug)printf("PPOP pos global=(%f,%f,%f) --> local (%f,%f,%f)\n", pos.x(), pos.y(), pos.z(), x.x(), x.y(), x.z());           
           if(debug)printf("PPOP mom global=(%f,%f,%f) --> local (%f,%f,%f)\n", momentum.x(), momentum.y(), momentum.z(), p.x(), p.y(), p.z());           


           double s = -x.z()/p.z(); // sp.z() - x.z(); local z of plane always 0
           if(debug)printf("PPOP s=%f  --> can compute = %i\n", s, (p.z()==0 || (((p.x() != 0 || p.y() != 0) && p.z() == 0 && s!= 0)))==true?0:1);


           if (p.z()==0 || ((p.x() != 0 || p.y() != 0) && p.z() == 0 && s!= 0)) return false;
 
           localPosOnPlane = TVector3(x.x() + p.x()*s,
                                      x.y() + p.y()*s,
                                      x.z() + p.z()*s);    

           if(debug)printf("PPOS LPOP = (%f,%f,%f) --> mag=%f\n", x.x() + p.x()*s, x.y() + p.y()*s, x.z() + p.z()*s, localPosOnPlane.Mag());

           return true;
        }

        static ModuleGeom* get(unsigned int detId){return static_geomMap[detId]; }

        static void loadGeometry(std::string path){
            ModuleGeom::static_geomMap.clear(); //reset the geometry map

            TChain* t = new TChain("GeomDumper/geom");
            t->Add(path.c_str());

            unsigned int rawId;                 t->SetBranchAddress("rawId", &rawId);
            float trapezeParam;                 t->SetBranchAddress("trapezeParam", &trapezeParam);
            TVector3* posV    = new TVector3(); t->SetBranchAddress("pos",    &posV);
            TVector3* widthV  = new TVector3(); t->SetBranchAddress("width",  &widthV);
            TVector3* lengthV = new TVector3(); t->SetBranchAddress("length", &lengthV);
            TVector3* thickV  = new TVector3(); t->SetBranchAddress("thick",  &thickV);

            for (unsigned int ientry = 0; ientry < t->GetEntries(); ientry++) {
                t->GetEntry(ientry);
                
                ModuleGeom* mod = new ModuleGeom();
                mod->TP     = trapezeParam;
                mod->pos    = *posV;
                mod->width  = *widthV;
                mod->length = *lengthV;
                mod->thick  = *thickV;
                static_geomMap[rawId] = mod;
            }
            delete t;
        }
};
std::unordered_map<unsigned int, ModuleGeom*> ModuleGeom::static_geomMap; //need to define this here to reference the object


#endif
