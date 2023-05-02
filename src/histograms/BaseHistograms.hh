#ifndef BaseHistograms_hh
#define BaseHistograms_hh

#include <map>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "Particle.hh"

class BaseHistograms
{
public:
    BaseHistograms(const std::string &suffix);
    virtual ~BaseHistograms();
    virtual void Fill(const Particle &particle, const double &weight);
    virtual void Normalize(const double &scale);
    virtual void Write();

    std::string name, suffix;
    std::map<std::string, TH1D *> Histogram1D_Collection;
    std::map<std::string, TH2D *> Histogram2D_Collection;
    std::map<std::string, TH3D *> Histogram3D_Collection;

protected:
    std::vector<std::string> PARTICLENAMES = {
        "n",
        "p",
        "d",
        "t",
        "3He",
        "4He",
        "coal_n",
        "coal_p",
    };

    std::map<std::pair<int, int>, std::string> NUCLEINAMES = {
        {{0, 1}, "n"},
        {{1, 1}, "p"},
        {{1, 2}, "d"},
        {{1, 3}, "t"},
        {{2, 3}, "3He"},
        {{2, 4}, "4He"},
    };
};

class ImpactParameterMultiplicity : public BaseHistograms
{
public:
    ImpactParameterMultiplicity(const std::string &suffix);
    void Fill(const int &multi, const double &b, const double &weight);
};

class KinergyTheta : public BaseHistograms
{
public:
    KinergyTheta(const std::string &suffix, const std::string &frame = "cms");
    void Fill(const Particle &particle, const double &weight);

private:
    std::string frame;
};

class PtRapidity : public BaseHistograms
{
public:
    PtRapidity(const std::string &suffix);
    void Fill(const Particle &particle, const double &weight);
};

class PmagEmissionTime : public BaseHistograms
{
public:
    PmagEmissionTime(const std::string &suffix);
    void Fill(const Particle &particle, const double &weight);
};
#endif