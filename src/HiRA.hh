#ifndef HiRA_hh
#define HiRA_hh

#include <iostream>
#include <fstream>
#include <array>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

#include "TMath.h"
#include "TString.h"

class HiRA
{
protected:
    // kinergy cut for experiment e15190
    std::map<std::string, std::array<double, 2>> KINERGYCUT = {
        {"p", {20.0, 198.0}},
        {"d", {15.0, 263.0 / 2}},
        {"t", {12.0, 312 / 3.}},
        {"3He", {20.0, 200.0}},
        {"4He", {18.0, 200.0}},
    };

    std::map<std::pair<int, int>, std::string> MASS_CHARGE_TO_NAME = {
        {{1, 1}, "p"},
        {{2, 1}, "d"},
        {{3, 1}, "t"},
        {{3, 2}, "3He"},
        {{4, 2}, "4He"},
    };
    std::array<double, 2> THETACUT = {30., 75.};
    std::array<double, 2> PHICUT = {0, 360};

public:
    HiRA();
    ~HiRA(){};
    bool PassCharged(const int &Z) { return Z > 0; }
    bool PassKinergyCut(const std::string &particle_name, const double &kinergy_per_nucleon);
    bool PassKinergyCut(const int &A, const int &Z, const double &kinergy);
    bool PassAngularCut(const double &theta_deg, const double &phi_deg);

    void AddAnalParticle(const int &Z, const int &A, const std::string &name) { this->mMassChargeToName[{A, Z}] = name; }
    void SetKinergyCut(const std::string &, const double &, const double &);
    void SetThetaCut(const double &, const double &);
    void SetPhiCut(const double &, const double &);

    void ResetCounter() { mCounterPass = mCounterFail = 0; }
    long GetCountPass() { return mCounterPass; }
    long GetCountFail() { return mCounterFail; }
    void CountPass() { mCounterPass++; }
    void CountFail() { mCounterFail++; }

private:
    long mCounterPass, mCounterFail;
    std::map<std::string, std::array<double, 2>> mKinergyCut;
    std::array<double, 2> mThetaCut;
    std::array<double, 2> mPhiCut;
    std::map<std::pair<int, int>, std::string> mMassChargeToName;
};

inline void HiRA::SetKinergyCut(const std::string &name, const double &a, const double &b) { mKinergyCut[name] = {a, b}; }
inline void HiRA::SetThetaCut(const double &a, const double &b) { mThetaCut = {a, b}; }
inline void HiRA::SetPhiCut(const double &a, const double &b) { mPhiCut = {a, b}; }

#endif
