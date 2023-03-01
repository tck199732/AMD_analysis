#ifndef Microball_hh
#define Microball_hh

#include <iostream>
#include <fstream>
#include <array>
#include <tuple>
#include <sstream>
#include <map>
#include <float.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "TString.h"

class Microball
{
protected:
    static constexpr int NumRing = 10; // use ring 2, 3, 4, 5, 7, 8
    static constexpr int NumDet = 15;
    static constexpr int MaxA = 100;
    static constexpr int MaxZ = 100;

public:
    Microball();
    ~Microball(){};

    // option control
    void TurnOff_ChargedParticle() { this->Is_apply_cut_charged_particle = 0; }
    void TurnOff_MultipleHitCorrection() { this->Is_apply_cut_multiple_hit = 0; }
    void TurnOff_KinergyCorrection() { this->Is_apply_cut_kinergy = 0; }
    void TurnOff_CoverageCorrection() { this->Is_apply_cut_coverage = 0; }

    bool Get_Is_apply_cut_charged_particle() { return this->Is_apply_cut_charged_particle; }
    bool Get_Is_apply_cut_multiple_hit() { return this->Is_apply_cut_multiple_hit; }
    bool Get_Is_apply_cut_kinergy() { return this->Is_apply_cut_kinergy; }
    bool Get_Is_apply_cut_coverage() { return this->Is_apply_cut_coverage; }

    // read input information
    void ReadGeometryMap(const std::string &filename);
    void ReadThresholdKinergyMap(const std::string &filename);
    void ConfigurateSetup(const std::string &reaction, const std::string &filename);

    // Getters for angles
    double GetThetaMin(const int &ring, const int &det) { return this->ThetaMap[ring][0]; }
    double GetThetaMax(const int &ring, const int &det) { return this->ThetaMap[ring][1]; }
    double GetPhiMin(const int &ring, const int &det) { return this->PhiMap[{ring, det}][0]; }
    double GetPhiMax(const int &ring, const int &det) { return this->PhiMap[{ring, det}][1]; }

    // get phi limits for shifting phi calulated from atan2
    double GetPhiMinInRing(const int &ring);
    double GetPhiMaxInRing(const int &ring);

    // detector ID from angles
    std::pair<int, int> GetRingDetID(const double &thetalab, const double &phi);
    int GetRingID(const double &thetalab);
    int GetDetID(const double &thetalab, const double &phi);

    // Getters for threshold kinergy
    double GetThresholdKinergy(const int &ring_id, const int &aid, const int &zid);
    double GetThresholdKinergy(const double &thetalab, const int &aid, const int &zid);

    // CsI hit counting
    void ResetCsIHitMap();
    int GetCsIHits();
    int GetCsIHits(const int &ring, const int &det);
    void AddCsIHit(const double &thetalab, const double &phi);

    // filters
    bool IsChargedParticle(const int &Z);
    bool IsCovered(const double &thetalab, const double &phi);
    bool IsAccepted(const double &ekinlab, const double &thetalab, const int &aid, const int &zid);

    // check input information
    void ViewGeometryMap();
    void ViewThresholdKinergyMap();
    void ViewDetectorSetupMap();

private:
    std::map<std::array<int, 2>, int> CsIHitMap;
    std::map<int, std::array<double, 2>> ThetaMap;
    std::map<std::array<int, 2>, std::array<double, 2>> PhiMap;
    std::map<std::array<int, 3>, double> KinergyThresholdMap;
    std::map<int, std::vector<int>> DetectorSetupMap;

    bool Is_apply_cut_charged_particle; // if true, only count charged particle
    bool Is_apply_cut_multiple_hit;     // if true, only single hit on each CsI
    bool Is_apply_cut_kinergy;          // if true, remove particle with kinergy < threshold
    bool Is_apply_cut_coverage;         // if true, remove particle not covered by uBall
};

#endif
