#include "../src/Physics.cpp"
#include "../src/Microball.cpp"
#include "../src/HiRA.cpp"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include <stdlib.h>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

/**
 * @brief Input the N, Z and momenta (not per nucleon !)
 *
 */
struct particle
{
    int N, Z;
    double px, py, pz_cms;
    double kinergy_lab;
    double theta_lab, phi;
    void initialize(const double &betacms);
};

void particle::initialize(const double &betacms)
{
    double mass = Physics::GetNucleiMass(this->Z, this->N + this->Z);
    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));
    double pmag_trans = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    double pmag_cms = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz_cms, 2.));

    double kinergy_cms = TMath::Sqrt(pow(pmag_cms, 2.) + pow(mass, 2.)) - mass;

    double pz_lab = gamma * (this->pz_cms + betacms * (kinergy_cms + mass));
    double pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(pz_lab, 2.));

    this->theta_lab = TMath::ATan2(pmag_trans, pz_lab) * TMath::RadToDeg();

    this->kinergy_lab = TMath::Sqrt(pow(pmag_lab, 2.) + pow(mass, 2.)) - mass;
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg();

    if (this->phi < 0)
    {
        this->phi += 360.;
    }
    if (this->phi > 360)
    {
        this->phi -= 360.;
    }
}

Microball *GetMicroBall(const std::string &reaction)
{
    Microball *uBall = new Microball();
    fs::path project_dir = std::getenv("PROJECT_DIR");
    fs::path database_dir = project_dir / "database/e15190/microball/acceptance";
    fs::path path_config = database_dir / "config.dat";
    fs::path path_geometry = database_dir / "geometry.dat";

    std::map<int, fs::path> path_thres;
    path_thres[2] = database_dir / "threshold_Sn_65mgcm2.dat";
    path_thres[3] = database_dir / "threshold_Sn_58mgcm2.dat";
    path_thres[4] = database_dir / "threshold_Sn_50mgcm2.dat";
    path_thres[5] = database_dir / "threshold_Sn_43mgcm2.dat";
    path_thres[7] = database_dir / "threshold_Sn_30mgcm2.dat";
    path_thres[8] = database_dir / "threshold_Sn_23mgcm2.dat";

    uBall->ReadGeometry(path_geometry.string());

    for (auto &[id, pth] : path_thres)
    {
        if (!fs ::exists(pth))
        {
            std::string msg = Form("%s does not exists.", pth.c_str());
            throw std::invalid_argument(msg.c_str());
        }
        uBall->ReadThreshold(id, pth.string());
    }

    uBall->ReadConfiguration(reaction, path_config.string());

    uBall->HiRA_Coordinate(); // add phi angle by 90, and move to phi=[0,360];

    return uBall;
}

bool ReadMicroballParticle(Microball *&uBall, const particle &particle)
{
    return uBall->IsChargedParticle(particle.Z) &&
           uBall->IsCovered(particle.theta_lab, particle.phi) &&
           uBall->IsAccepted(particle.kinergy_lab, particle.theta_lab, particle.N + particle.Z, particle.Z) && uBall->IsReadyCsI(particle.theta_lab, particle.phi);
}

bool ReadHiRAParticle(HiRA *&hira, const particle &particle)
{
    return hira->PassAngularCut(particle.theta_lab, particle.phi) && hira->PassCharged(particle.Z) && hira->PassKinergyCut(particle.Z + particle.N, particle.Z, particle.kinergy_lab);
}