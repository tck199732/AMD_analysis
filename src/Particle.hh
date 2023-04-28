#ifndef Particle_hh
#define Particle_hh

#include "TMath.h"
#include "Physics.hh"

class Particle
{
public:
    Particle(const int &N, const int &Z, const double &px_per_nucleon, const double &py_per_nucleon, const double &pz_per_nucleon, const double &m = 0., const std::string &frame = "cms");
    ~Particle() { ; }

    void Initialize(const double &betacms, const double &beam_rapidity = 1.);

    int N, Z, A;
    double mass;

    // same in lab and cms
    double px, py, phi, pmag_trans;

    // cms quantities
    double pz_cms;
    double theta_cms, kinergy_cms, pmag_cms, rapidity_cms;

    // lab quantities
    double pz_lab;
    double theta_lab, kinergy_lab, pmag_lab, rapidity_lab;

    // rapidity lab / beam rapidity
    double rapidity_lab_normed;

private:
    std::string _frame_at_construct;

protected:
    double NucleonMass = 938.272; // MeV/c^2
};
#endif