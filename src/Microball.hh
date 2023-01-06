#ifndef Microball_hh
#define Microball_hh

#include <iostream>
#include <fstream>
#include <array>
#include <filesystem>
namespace fs = std::filesystem;

#include "TMath.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2D.h"

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

    void TurnOff_ChargedParticle() { this->OptionChargedParticle = 0; }
    void TurnOff_MultipleHitCorrection() { this->OptionMultipleHit = 0; }
    void TurnOff_KinergyCorrection() { this->OptionKinergy = 0; }
    void TurnOff_CoverageCorrection() { this->OptinoCoverage = 0; }

    void ReadGeometry(const std::string &filename);
    void ReadThreshold(const int &ring, const std::string &filename);
    void ReadConfiguration(const std::string &reaction, const std::string &filename);
    void RemoveCsI(const int &ring, const int &det);

    double GetThetaMin(const int &ring, const int &det) { return this->Theta[ring][det][0]; }
    double GetThetaMax(const int &ring, const int &det) { return this->Theta[ring][det][1]; }
    double GetPhiMin(const int &ring, const int &det) { return this->Phi[ring][det][0]; }
    double GetPhiMax(const int &ring, const int &det) { return this->Phi[ring][det][1]; }

    int GetRingID(const double &thetalab);
    int GetDetID(const double &thetalab, const double &phi);
    double GetThreshold(const double &thetalab, const int &aid, const int &zid);
    double GetThreshold(const int &ring_id, const int &aid, const int &zid) { return this->KinergyThreshold[ring_id][aid][zid]; }

    void ResetCsI();
    int GetCsIHits();
    int GetCsIHits(const int &ring, const int &det) { return this->CsIHits[ring][det]; }
    void AddCsIHit(const int &ring, const int &det) { this->CsIHits[ring][det]++; }
    void AddCsIHit(const double &thetalab, const double &phi);
    void HiRA_Coordinate();
    void ResetPhiRange();

    bool IsChargedParticle(const int &Z);
    bool IsCovered(const double &thetalab, const double &phi);
    bool IsAccepted(const double &ekinlab, const double &thetalab, const int &aid, const int &zid);
    bool IsReadyCsI(const double &thetalab, const double &phi);

private:
    int CsIHits[NumRing][NumDet];                 // number of hit on each CsI detector
    std::array<double, 2> Theta[NumRing][NumDet]; // range of theta of a CsI detector
    std::array<double, 2> Phi[NumRing][NumDet];   // range of phi
    double KinergyThreshold[NumRing][MaxA][MaxZ]; // threshold of kinergy-lab

    bool OptionChargedParticle; // if true, only count charged particle
    bool OptionMultipleHit;     // if true, only single hit on each CsI
    bool OptionKinergy;         // if true, remove particle with kinergy < threshold
    bool OptinoCoverage;        // if true, remove particle not covered by uBall
};

#endif
