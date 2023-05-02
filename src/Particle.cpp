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

    // by default, x, y, z, t are set to DBL_MIN
    this->x = DBL_MIN;
    this->y = DBL_MIN;
    this->z_cms = DBL_MIN;
    this->t_cms = DBL_MIN;
    this->z_lab = DBL_MIN;
    this->t_lab = DBL_MIN;
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

        //
        if (this->t_cms >= 0. && this->z_cms >= 0.)
        {
            this->t_lab = gamma * (this->t_cms + betacms * this->z_cms);
            this->z_lab = gamma * (this->z_cms + betacms * this->t_cms);
        }
    }
    else if (this->_frame_at_construct == "lab")
    {
        // construct cms quantities
        this->pz_cms = Physics::boostz(this->mass, this->pz_lab, this->kinergy_lab, betacms);
        this->pmag_cms = Physics::GetP(this->pmag_trans, this->pz_cms);
        this->kinergy_cms = Physics::GetEkin(this->mass, this->pmag_cms);
        this->theta_cms = Physics::GetTheta(this->pmag_trans, this->pz_cms);
        this->rapidity_cms = Physics::GetRapidity(this->kinergy_cms, this->pz_cms, this->mass);

        if (this->t_lab >= 0. && z_lab >= 0.)
        {
            this->t_cms = gamma * (this->t_lab - betacms * this->z_lab);
            this->z_cms = gamma * (this->z_lab - betacms * this->t_lab);
        }
    }

    // normalize rapidity to beam rapidity
    this->rapidity_lab_normed = this->rapidity_lab / beam_rapidity;
}

void Particle::SetXYZT(const double &x, const double &y, const double &z, const double &t, const std::string &frame = "cms")
{
    this->x = x;
    this->y = y;

    if (frame == "cms")
    {
        this->z_cms = z;
        this->t_cms = t;
    }
    else if (frame == "lab")
    {
        this->z_lab = z;
        this->t_lab = t;
    }
}