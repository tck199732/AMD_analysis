#include "HiRA.hh"

HiRA::HiRA()
{
    mKinergyCut = this->KINERGYCUT;
    mThetaCut = this->THETACUT;
    mPhiCut = this->PHICUT;
    mMassChargeToName = this->MASS_CHARGE_TO_NAME;
    mCounterPass = 0;
    mCounterFail = 0;
}

bool HiRA::PassAngularCut(const double &theta_deg, const double &phi_deg)
{
    return (
        theta_deg >= mThetaCut[0] && theta_deg <= mThetaCut[1] &&
        phi_deg >= mPhiCut[0] && phi_deg <= mPhiCut[1]);
}

bool HiRA::PassKinergyCut(const std::string &particle_name, const double &kinergy_per_nucleon)
{
    return this->mKinergyCut.count(particle_name) == 1 && kinergy_per_nucleon >= this->mKinergyCut[particle_name][0] && kinergy_per_nucleon <= this->mKinergyCut[particle_name][1];
}

bool HiRA::PassKinergyCut(const int &A, const int &Z, const double &kinergy)
{
    if (mMassChargeToName.count({A, Z}) == 0)
    {
        return false;
    }
    return PassKinergyCut(mMassChargeToName[{A, Z}], kinergy / A);
}