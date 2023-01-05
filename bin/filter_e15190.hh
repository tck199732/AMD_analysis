#include "../src/system_info.cpp"
#include "../src/detector_uball.cpp"
#include "../src/particle.cpp"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
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

class MicroBall
{
public:
    MicroBall();
    ~MicroBall();
    void Configurate(const std::string &system);
    void Reset();
    void ReadParticle(const particle &particle);
    int GetNumberOfChargedParticles() { return this->counter; }

private:
    fs::path path_config, path_geo;
    std::map<int, fs::path> path_thres;
    detector_uball *uball;
    int counter = 0;
};

MicroBall::MicroBall()
{
    fs::path path_base = "../database/e15190/microball/acceptance";
    this->path_config = path_base / "config.dat";
    this->path_geo = path_base / "geometry.dat";
    this->path_thres[2] = path_base / "threshold_Sn_65mgcm2.dat";
    this->path_thres[3] = path_base / "threshold_Sn_58mgcm2.dat";
    this->path_thres[4] = path_base / "threshold_Sn_50mgcm2.dat";
    this->path_thres[5] = path_base / "threshold_Sn_43mgcm2.dat";
    this->path_thres[7] = path_base / "threshold_Sn_30mgcm2.dat";
    this->path_thres[8] = path_base / "threshold_Sn_23mgcm2.dat";

    this->uball = new detector_uball();
    this->uball->read_geometry(fs::absolute(this->path_geo));

    for (auto &[id, pth] : this->path_thres)
    {
        if (!fs ::exists(pth))
        {
            std::string msg = Form("%s does not exists.", pth.c_str());
            throw std::invalid_argument(msg.c_str());
        }
        this->uball->read_threshold(id, fs::absolute(pth));
    }

    this->uball->hira_coordinate(); // add phi angle by 90, and move to phi=[0,360];
}

MicroBall::~MicroBall() { delete this->uball; }

void MicroBall::Configurate(const std::string &system)
{
    std::ifstream stream(fs::absolute(this->path_config));
    stream.ignore(99, '\n');
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream iss(line);
        std::string sys;
        int ring_id, det_id;
        iss >> sys >> ring_id;
        if (sys != system)
        {
            continue;
        }

        std::vector<int> det_array;
        while (iss >> det_id)
        {
            det_array.push_back(det_id);
        }

        for (int i = 0; i < 15; i++)
        {
            if (std::find(det_array.begin(), det_array.end(), i) - det_array.begin() == det_array.size())
            {
                this->uball->remove_csi(ring_id, i);
            }
        }
    }
}

void MicroBall::Reset()
{
    this->uball->reset_csi();
    this->counter = 0;
}

void MicroBall::ReadParticle(const particle &particle)
{
    if (particle.zid == 0)
    {
        return;
    }
    if (this->uball->cover(particle.thetalab, particle.phi))
    {
        if (this->uball->punch_thr(particle.ekinlab, particle.thetalab, particle.aid, particle.zid))
        {
            if (this->uball->ready_csi(particle.thetalab, particle.phi))
            {
                this->counter++;
                this->uball->add_csi_hit(particle.thetalab, particle.phi);
            }
        }
    }
    return;
}
