#include "functions.h"
#include "TMath.h"

float deltaPhi(float phi1, float phi2) {
  float r = phi1 - phi2;
  float pi = TMath::Pi();

  if (r > 2 * pi)
    r -= 2 * pi;
  if (r < -2 * pi)
    r += 2 * pi;
  return r;
}

float deltaR(float phi1, float phi2, float eta1, float eta2) {
  float dphi = deltaPhi(phi1, phi2);

  return TMath::Sqrt(dphi * dphi + (eta1 - eta2) * (eta1 - eta2));
}
