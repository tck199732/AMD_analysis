#ifndef Physics_hh
#define Physics_hh

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <regex>

#include "TMath.h"
#include "TString.h"

namespace Physics
{

    double GetReactionBeta(const double &mass1, const double &mass2, const double &beam_energy_per_nucleon, const int &beam_nucleon);
    double GetBeamRapidity(const double &mass1, const double &mass2, const double &beam_energy_per_nucleon, const int &beam_nucleon);

    // general physics
    float GetPt(const float &px, const float &py);
    float GetP(const float &pt, const float &pz);
    float GetP(const float &px, const float &py, const float &pz);
    float GetEkin(const float &mass, const float &px, const float &py, const float &pz);
    float GetEkin(const float &mass, const float &p);
    float GetPhi(const float &px, const float &py);
    float GetTheta(const float &pt, const float &pz);
    float boostz(const float &mass, const float &pz, const float &ekin, const float &betacms);
    float GetRapidity(const float &ekin, const float &pz, const float &mass);

};

#endif