#include "Physics.hh"

double Physics::GetReactionBeta(const double &mass1, const double &mass2, const double &beam_energy_per_nucleon, const int &beam_nucleon)
{
    double beam_ke = beam_energy_per_nucleon * 1.0 * beam_nucleon;
    double beam_energy_tot = beam_ke + mass1;
    double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * mass1);
    double gamma = beam_energy_tot / mass1;
    return mom_beam / (gamma * mass1 + mass2);
}

double Physics::GetBeamRapidity(const double &mass1, const double &mass2, const double &beam_energy_per_nucleon, const int &beam_nucleon)
{
    double beam_ke = beam_energy_per_nucleon * 1.0 * beam_nucleon;
    double beam_energy_tot = beam_ke + mass1;
    double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * mass1);
    return 0.5 * TMath::Log((beam_energy_tot + mom_beam) / (beam_energy_tot - mom_beam));
}

float Physics::GetPt(const float &px, const float &py)
{
    return TMath::Sqrt(pow(px, 2.) + pow(py, 2.));
}

float Physics::GetP(const float &pt, const float &pz)
{
    return TMath::Sqrt(pow(pt, 2.) + pow(pz, 2.));
}

float Physics::GetP(const float &px, const float &py, const float &pz)
{
    return Physics::GetP(Physics::GetPt(px, py), pz);
}

float Physics::GetEkin(const float &mass, const float &p)
{
    return TMath::Sqrt(pow(p, 2.) + pow(mass, 2.)) - mass;
}

float Physics::GetEkin(const float &mass, const float &px, const float &py,
                       const float &pz)
{
    float p = Physics::GetP(px, py, pz);
    return Physics::GetEkin(mass, p);
}

float Physics::GetPhi(const float &px, const float &py)
{
    return TMath::ATan2(py, px);
}

float Physics::GetTheta(const float &pt, const float &pz)
{
    return TMath::ATan2(pt, pz);
}

/**
 * @brief Boost from lab framm to cms frame if betacms is positive
 *
 * @param mass
 * @param pz
 * @param ekin
 * @param betacms
 * @return float
 */
float Physics::boostz(const float &mass, const float &pz, const float &ekin,
                      const float &betacms)
{
    float gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));
    return gamma * (pz - betacms * (ekin + mass));
}

float Physics::GetRapidity(const float &ekin, const float &pz,
                           const float &mass)
{
    return 0.5 * TMath::Log((ekin + pz + mass) / (ekin - pz + mass));
}