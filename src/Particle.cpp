#include "Particle.hh"

Particle::Particle(const int &N, const int &Z, const double &px_per_nucleon, const double &py_per_nucleon, const double &pz_per_nucleon, const double &m, const std::string &frame)
{
    this->N = N;
    this->Z = Z;
    this->A = N + Z;

    this->px = px_per_nucleon * A;
    this->py = py_per_nucleon * A;

    this->mass = m;

    if (m == 0.)
    {
        this->mass = this->A * this->NucleonMass;
    }
    this->_frame_at_construct = frame;

    // initialize frame-independent quantities
    this->phi = Physics::GetPhi(this->px, this->py);
    this->pmag_trans = Physics::GetPt(this->px, this->py);

    if (frame == "cms")
    {
        // initialize cms quantities
        this->pz_cms = pz_per_nucleon * A;
        this->pmag_cms = Physics::GetP(this->pmag_trans, this->pz_cms);
        this->kinergy_cms = Physics::GetEkin(this->mass, this->pmag_cms);
        this->theta_cms = Physics::GetTheta(this->pmag_trans, this->pz_cms);
        this->rapidity_cms = Physics::GetRapidity(this->kinergy_cms, this->pz_cms, this->mass);
    }

    else if (frame == "lab")
    {
        // initialize lab quantities
        this->pz_lab = pz_per_nucleon * A;
        this->pmag_lab = Physics::GetP(this->pmag_trans, this->pz_lab);
        this->kinergy_lab = Physics::GetEkin(this->mass, this->pmag_lab);
        this->theta_lab = Physics::GetTheta(this->pmag_trans, this->pz_lab);
        this->rapidity_lab = Physics::GetRapidity(this->kinergy_lab, this->pz_lab, this->mass);
    }
}

void Particle::Initialize(const double &betacms, const double &beam_rapidity)
{
    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));

    if (this->_frame_at_construct == "cms")
    {
        // construct lab quantities
        this->pz_lab = Physics::boostz(this->mass, this->pz_cms, this->kinergy_cms, -betacms);
        this->pmag_lab = Physics::GetP(this->pmag_trans, this->pz_lab);
        this->kinergy_lab = Physics::GetEkin(this->mass, this->pmag_lab);
        this->theta_lab = Physics::GetTheta(this->pmag_trans, this->pz_lab);
        this->rapidity_lab = Physics::GetRapidity(this->kinergy_lab, this->pz_lab, this->mass);
    }
    else if (this->_frame_at_construct == "lab")
    {
        // construct cms quantities
        this->pz_cms = Physics::boostz(this->mass, this->pz_lab, this->kinergy_lab, betacms);
        this->pmag_cms = Physics::GetP(this->pmag_trans, this->pz_cms);
        this->kinergy_cms = Physics::GetEkin(this->mass, this->pmag_cms);
        this->theta_cms = Physics::GetTheta(this->pmag_trans, this->pz_cms);
        this->rapidity_cms = Physics::GetRapidity(this->kinergy_cms, this->pz_cms, this->mass);
    }

    // normalize rapidity to beam rapidity
    this->rapidity_lab_normed = this->rapidity_lab / beam_rapidity;
}