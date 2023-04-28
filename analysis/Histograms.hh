#ifndef Histograms_hh
#define Histograms_hh

#include <map>
#include <string>
#include <vector>

#include "TH2D.h"
#include "Particle.hh"

class Histograms
{
public:
    Histograms(const std::string &suffix);
    ~Histograms() { ; }
    void Fill(const Particle &particle, const double &weight);
    void Normalize(const double &scale);
    void Write();
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;

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
#endif