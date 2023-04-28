#include "Particle.hh"

Particle::Particle(const int &N, const int &Z, const double &px_per_nucleon, const double &py_per_nucleon, const double &pz_cms_per_nucleon, const double &m)
{
    this->N = N;
    this->Z = Z;
    this->A = N + Z;

    this->px = px_per_nucleon * A;
    this->py = py_per_nucleon * A;
    this->pz_cms = pz_cms_per_nucleon * A;

    this->mass = m;

    if (m == 0.)
    {
        this->mass = this->A * this->NucleonMass;
    }

    // initialize frame-independent quantities
    this->phi = TMath::ATan2(this->py, this->px);
    this->pmag_trans = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
}

void Particle::Initialize(const double &betacms)
{
    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));

    // initialize cms quantities
    this->pmag_cms = TMath::Sqrt(pow(this->pmag_trans, 2.) + pow(this->pz_cms, 2.));
    this->kinergy_cms = TMath::Sqrt(pow(this->pmag_cms, 2.) + pow(this->mass, 2.)) - this->mass;
    this->theta_cms = TMath::ATan2(this->pmag_trans, this->pz_cms) * TMath::RadToDeg();
    this->rapidity_cms = 0.5 * TMath::Log((this->pmag_cms + this->pz_cms) / (this->pmag_cms - this->pz_cms));

    // initialize lab quantities
    this->pz_lab = gamma * (this->pz_cms + betacms * (this->kinergy_cms + this->mass));
    double pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz_lab, 2.));
    this->kinergy_lab = TMath::Sqrt(pow(pmag_lab, 2.) + pow(this->mass, 2.)) - this->mass;
    this->theta_lab = TMath::ATan2(pmag_trans, this->pz_lab) * TMath::RadToDeg();
    this->rapidity_lab = 0.5 * TMath::Log((pmag_lab + this->pz_lab) / (pmag_lab - this->pz_lab));
}