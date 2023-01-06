#include "../src/Physics.cpp"
#include "../src/particle.cpp"
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
    return uBall->IsChargedParticle(particle.zid) &&
           uBall->IsCovered(particle.thetalab, particle.phi) &&
           uBall->IsAccepted(particle.ekinlab, particle.thetalab, particle.aid, particle.zid) && uBall->IsReadyCsI(particle.thetalab, particle.phi);
}

bool ReadHiRAParticle(HiRA *&hira, const particle &particle)
{
    return hira->PassAngularCut(particle.thetalab, particle.phi) && hira->PassCharged(particle.zid) && hira->PassKinergyCut(particle.name, particle.ekinlab);
}