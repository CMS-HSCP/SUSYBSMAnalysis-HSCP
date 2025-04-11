// Original Author:  Gael Coulon

#ifndef SUSYBSMAnalysis_Analyzer_SaturationCorrectionInStrip_h
#define SUSYBSMAnalysis_Analyzer_SaturationCorrectionInStrip_h

#include <numeric> // pour std::accumulate
#include <vector>
#include <iostream>
using namespace std;

// Layer:
//       1-4: TIB 
//      5-10: TOB
//     11-13: TEC ring
//     14-20: TEC ring

int Correction_FL_FR(const std::vector <int>&  Q, int layer, std::string TemplateFile)
{
    int N_sat = 0;
    int ThresholdSat = -1;

    if ((layer>=5 && layer<=10) || layer>=18) ThresholdSat = 25;
    else ThresholdSat = 16;

    // CORRECTION TEMPLATES
    std::string line;
    std::vector<double> template_FL;
    std::vector<double> template_FR;

    std::ifstream Template(TemplateFile);
    while (std::getline(Template, line))
    {
        std::istringstream iss(line);
        int i_layer=-1;
        double FL_a=-1, FR_a=-1;
        if (!(iss >> i_layer >> FL_a >> FR_a)) break; // error
            
        template_FL.push_back(FL_a);
        template_FR.push_back(FR_a);
    }
    Template.close();

    // NUMBER OF MAX
    for (unsigned int i=0; i<Q.size(); i++)
    {
        if (Q[i] >= 254) N_sat++;
    }

    int max_Q = *max_element(Q.begin(), Q.end());
    int i_max = find(Q.begin(), Q.end(), max_Q) - Q.begin();

    // NO CORRECTION IF TOO SMALL OR LARGE
    if (Q.size()<2 || Q.size()>8) return max_Q;

    if (N_sat == 1)
    {
        int MaxCorr = -1;
        // FULL LEFT
        if (max_Q==Q[0] && Q[i_max+1]>=ThresholdSat && Q[i_max+1]<254)
        {
            MaxCorr = template_FL[layer-1] * Q[i_max+1];

            if (max_Q == 255 && MaxCorr < 255) return 255;
            if (MaxCorr < 254) return 254;
            else return MaxCorr;
        }

        // FULL RIGHT
        if (max_Q==Q[Q.size()-1] && Q[i_max-1]>=ThresholdSat && Q[i_max-1]<254)
        {
            MaxCorr = template_FR[layer-1] * Q[i_max-1];

            if (max_Q == 255 && MaxCorr < 255) return 255;
            if (MaxCorr < 254) return 254;
            else return MaxCorr;
        }

        return max_Q;
    }

    // NO ONE SINGLE SATURATED STRIP --> no x-talk inversion
    return max_Q;
}

int Correction_LRC(const std::vector <int>&  Q, int layer, std::string TemplateFile, bool CENTER)
{
	// SETUP 
	int N_sat = 0;
	int Vmin = 0, Vmax = 0;

    for (unsigned int i=0; i<Q.size(); i++)
    {
        if (Q[i]>=254) N_sat++;
    }

	int max_Q = *max_element(Q.begin(), Q.end());
    int i_max = find(Q.begin(), Q.end(), max_Q) - Q.begin(); 

    if (Q[i_max-1] > Q[i_max+1])
    {
        Vmax = Q[i_max-1];
        Vmin = Q[i_max+1];
    }
    else
    {
        Vmax = Q[i_max+1];
        Vmin = Q[i_max-1];
    }

    // NO CORRECTION IF TOO SMALL OR LARGE
    if (Q.size()<2 || Q.size()>8) return max_Q;

    // CORRECTION TEMPLATES
    std::ifstream Template(TemplateFile);
    std::string line;
    std::vector<double> template_a1;
    std::vector<double> template_a2;
    
    if (CENTER)
    {
        while (std::getline(Template, line))
        {
            std::istringstream iss(line);
            int i_layer=-1;
            double a1=-1;
            if (!(iss >> i_layer >> a1)) break; // error
                
            template_a1.push_back(a1);
        }
        Template.close();

        int MaxCorr = 1.*(Vmax+Vmin)/2 * template_a1[layer-1];

        if (max_Q == 255 && MaxCorr < 255) return 255;
        if (MaxCorr < 254) return 254;
        else return MaxCorr;
    }

    while (std::getline(Template, line))
    {
        std::istringstream iss(line);
        int i_layer=-1;
        double a1=-1, a2=-1;
        if (!(iss >> i_layer >> a1 >> a2 )) break; // error
        
        template_a1.push_back(a1);
        template_a2.push_back(a2);
    }
    Template.close();

    int MaxCorr = Vmax * template_a1[layer-1] + Vmin * template_a2[layer-1];
    if (max_Q == 255 && MaxCorr < 255) return 255;
    if (MaxCorr < 254) return 254;
    else return MaxCorr;
}

int Correction_2strips(const std::vector <int>&  Q)
{
    int Qcorr = -1;
    int SumQ = accumulate(Q.begin(), Q.end(), 0);
    unsigned int index_maxLeft = -1;

        // NUMBER OF MAX + IF CONSECUTIVE OR NOT
    int N_sat = 0;
    bool two_cons_sat = false;

    if (Q[0]>=254) N_sat++;
    for (unsigned int i=1; i<Q.size(); i++)
    {
        if (Q[i]>=254) N_sat++;

        if (Q[i-1]>=254 && Q[i]>=254)
        {
            if (two_cons_sat) continue;
            two_cons_sat = true;
            index_maxLeft = i-1; 
        }
    }
    
        // NO CORRECTION if not 2 consecutive saturated strips, or if clusters too small or if max is at the edge
    if (N_sat!=2 || !two_cons_sat || Q.size() <= 3
        || index_maxLeft==0 || index_maxLeft==Q.size()-2) return SumQ;
    
        // CORRECTION
    Qcorr = 1./0.0613 * (Q[index_maxLeft-1] + Q[index_maxLeft+2])/2;

    return Qcorr;
}

void ClusterShape(const std::vector <int>& Q, bool &left, bool &right, bool &center, bool &FullLeft, bool &FullRight)
{
    // SETUP
    if (Q.empty()) return;

    auto it_max = std::max_element(Q.begin(), Q.end());
    if (it_max == Q.end()) return;

    int max_Q = *it_max;
    auto it_found = std::find(Q.begin(), Q.end(), max_Q);
    if (it_found == Q.end()) return;

    unsigned int i_max = it_found - Q.begin();

    int Nleft = -1, Nright = -1;
    if (Q.size()>=3 && i_max>0 && i_max<Q.size()-1) {Nleft = Q[i_max-1]; Nright = Q[i_max+1];}
    
    int N_sat = 0;
    for (unsigned int i=0; i<Q.size(); i++) { if (Q[i]>=254) N_sat++; }

    // SHAPE
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft>1.1*Nright) left=true;
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft<0.9*Nright) right=true;
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
    if (Q.size()>=2 && N_sat==1 && i_max==0) FullLeft=true;
    if (Q.size()>=2 && N_sat==1 && i_max==Q.size()-1) FullRight=true;

    return;
}

std::vector <int> ReturnCorrVec(const std::vector <int>& Q, const int layer, bool& AreSameCluster)
{
    // SETUP
    if (Q.empty()) return {};
    bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
    int MaxCorr = -1;


    // SHAPE
    ClusterShape(Q, left, right, center, FullLeft, FullRight);
    

    // CORRECTION
        // 2 consecutive saturated strips
    MaxCorr = Correction_2strips(Q);
    if (MaxCorr > accumulate(Q.begin(), Q.end(), 0))
    {
        std::vector <int> Qcorr2s;
        Qcorr2s.push_back(MaxCorr);

        AreSameCluster = false;
        return Qcorr2s;
    }
    
    MaxCorr = *max_element(Q.begin(), Q.end());
        // 1 saturated strip
    if (center) MaxCorr = Correction_LRC(Q, layer, "SUSYBSMAnalysis/HSCP/data/Template_CENTER.txt", true);
    else if (left || right) MaxCorr = Correction_LRC(Q, layer, "SUSYBSMAnalysis/HSCP/data/Template_LEFTRIGHT.txt", false);
    else if (FullLeft || FullRight) MaxCorr = Correction_FL_FR(Q, layer, "SUSYBSMAnalysis/HSCP/data/Template_FLFR.txt");
    

    // SUM CORRECTION
    std::vector <int> Qcorr;
    unsigned int i_max = find(Q.begin(), Q.end(), *max_element(Q.begin(), Q.end())) - Q.begin();
    for (unsigned int i=0; i<Q.size(); i++)
    {
        if (i != i_max) Qcorr.push_back(Q[i]);
        else Qcorr.push_back(MaxCorr);
    }

    if (accumulate(Qcorr.begin(), Qcorr.end(), 0) <= accumulate(Q.begin(), Q.end(), 0)) AreSameCluster = true;
    else AreSameCluster = false;
    
    return Qcorr;
}

std::vector <int> CrossTalkInvInStrip(const std::vector<int>& Q, const int layer, const std::string TemplateFile, bool IfCorrApplied = true, float threshold = 20) {

      // CROSS-TALK COEFF.      
  std::ifstream Template(TemplateFile);
  std::string line;
  std::vector <double> x1;
  std::vector <double> x2;
  
  while (std::getline(Template, line))
  {
      std::istringstream iss(line);
      int i_layer=-1;
      double x1_value=-1, x2_value=-1;
      if (!(iss >> i_layer >> x1_value >> x2_value)) break; // error
          
      x1.push_back(x1_value);
      x2.push_back(x2_value);
  }
  Template.close();


      // EXCLUSION PART: ANOMALOUS SHAPES or SATURATION
  if (Q.empty()) return {};
  std::vector<int> QII;
  if(Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    return QII;
  }

  if(IfCorrApplied)
  {
    int N_sat = 0;
    for (unsigned int i=0; i<Q.size(); i++) { if (Q[i]>=254) N_sat++;}
    
    if(N_sat>0)
    {
      bool AreSameCluster = true;
      QII = ReturnCorrVec(Q, layer, AreSameCluster);
      if (!AreSameCluster) return QII;    // saturation correction
    }
  }


      // CROSS-TALK INVERSION
  QII.clear();
  std::vector <double> a;
  for (unsigned int i = 0; i < x1.size(); i++) a.push_back(1 - 2*x1[i] - 2*x2[i]);
  
  const unsigned N = Q.size();
  std::vector<float> QI(N, 0);
  TMatrix A(N, N);

  for (unsigned int i = 0; i < N; i++) {
    A(i, i) = a[layer-1];
    if (i < N - 1) {
      A(i + 1, i) = x1[layer-1];
      A(i, i + 1) = x1[layer-1];
    }
    else continue;
    if (i < N - 2) {
      A(i + 2, i) = x2[layer-1];
      A(i, i + 2) = x2[layer-1];
    }
  }
  
  //    A =   a  x1 x2 0-------------0
  //          x1  a x1 x2 0          |
  //          x2 x1  a x1 x2 0       |
  //          0 x2 x1  a x1 x2 0     |
  //          |  0               0   |
  //          |     0              0 |
  //          |        0             0
  //          0----------0 x2 x1  a x1
  

  if (N == 1) A(0, 0) = 1 / a[layer-1];
  else A.InvertFast();

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      QI[i] += A(i, j) * (float)Q[j];
    }
  }

  for (unsigned int i = 0; i < QI.size(); i++) {
    if (QI[i] < threshold) QI[i] = 0;
    QII.push_back((int)QI[i]);
  }

  return QII;
}

#endif