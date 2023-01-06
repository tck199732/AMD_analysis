#include "../src/system_info.cpp"
#include "../src/particle.cpp"
#include "../src/Microball.cpp"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include <stdlib.h>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

// struct VetoWall
// {
//     std::map<std::string, std::array<double, 2>> Ekinlabcut;
//     std::map<std::string, std::array<double, 2>> thetalabcut;
//     std::map<std::string, std::array<double, 2>> philabcut;
//     void init();
//     bool pass_ekinlab(const particle &particle);
//     bool pass_angle(const particle &particle);
// };

// struct NeutronWall
// {
// };

struct HiRA
{
    std::map<std::string, std::array<double, 2>> Ekinlabcut = {
        {"p", {20.0, 198.0}},
        {"d", {15.0, 263.0 / 2}},
        {"t", {12.0, 312 / 3.}},
        {"3He", {20.0, 200.0}},
        {"4He", {18.0, 200.0}},
    };
    std::array<double, 2> thetalabcut = {30., 75.};
    std::array<double, 2> philabcut = {0, 360};
    bool pass_ekinlab(const particle &particle);
    bool pass_angle(const particle &particle);
};

bool HiRA::pass_ekinlab(const particle &particle)
{
    if (this->Ekinlabcut.count(particle.name) == 0)
    {
        return 0;
    }
    return (particle.ekinlab >= this->Ekinlabcut[particle.name][0] && particle.ekinlab <= this->Ekinlabcut[particle.name][1]);
}
bool HiRA::pass_angle(const particle &particle)
{
    return (particle.thetalab >= this->thetalabcut[0] && particle.thetalab <= this->thetalabcut[1] && particle.phi >= this->philabcut[0] && particle.phi <= this->philabcut[1]);
}

HiRA hira_detector;

Microball *GetMicroBall(const std::string &reaction)
{
    Microball *uBall = new Microball();
    fs::path project_dir = std::getenv("PROJECT_DIR");
    fs::path database_dir = project_dir / "database/e15190/microball/acceptance";
    fs::path path_config = path_base / "config.dat";
    fs::path path_geo = path_base / "geometry.dat";

    std::map<int, fs::path> path_thres;
    path_thres[2] = path_base / "threshold_Sn_65mgcm2.dat";
    path_thres[3] = path_base / "threshold_Sn_58mgcm2.dat";
    path_thres[4] = path_base / "threshold_Sn_50mgcm2.dat";
    path_thres[5] = path_base / "threshold_Sn_43mgcm2.dat";
    path_thres[7] = path_base / "threshold_Sn_30mgcm2.dat";
    path_thres[8] = path_base / "threshold_Sn_23mgcm2.dat";

    uBall->ReadGeometry(path_geo.string());

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

void ReadParticle(Microball *&uBall, const particle &particle)
{
    if (
        uBall->IsChargedParticle(particle.zid) &&
        uBall->IsCovered(particle.thetalab, particle.phi) &&
        uBall->IsAccepted(particle.ekinlab, particle.thetalab, particle.aid, particle.zid) && uBall->IsReadyCsI(particle.thetalab, particle.phi))
    {
        uBall->AddCsIHit(particle.thetalab, particle.phi);
    }
    return;
}
