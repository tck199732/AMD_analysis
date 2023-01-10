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

    extern std::map<std::string, double> particle_mass;
    extern std::map<std::string, double> particle_mass_old;
    extern double amu_MeV;
    extern std::map<std::string, int> particle_A;
    extern std::map<std::string, int> particle_Z;

    std::string GetNucleiName(const int &Z, const int &A);
    double GetBeamRapidity(const std::string &reaction, const std::string &option = "new");
    double GetReactionBeta(const std::string &reaction, const std::string &option = "new");
    double GetBeamRapidity(const std::string &beam, const int &beamEperA, const std::string &option = "new");
    double GetReactionBeta(const std::string &beam, const std::string &target, const int &beamEperA, const std::string &option = "new");
    double GetNucleiMass(const std::string &symbol, const std::string &option = "new");
    double GetNucleiMass(const int &Z, const int &A, const std::string &option = "new");

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