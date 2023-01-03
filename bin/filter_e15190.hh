#include "../src/system_info.cpp"
#include "../src/detector_uball.cpp"
#include "../src/particle.cpp"
#include "../src/root_io.cpp"

#include <vector>
#include <string>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

fs::path PROJECT_DIR = "/home/kin/Desktop/amd_analysis";

struct event
{
    int Nc;
    double bimp;
    std::vector<particle> particles; // hira particles
};

struct vetowall
{
    std::map<std::string, std::array<double, 2>> Ekinlabcut;
    std::map<std::string, std::array<double, 2>> thetalabcut;
    std::map<std::string, std::array<double, 2>> philabcut;
    void init();
    bool pass_ekinlab(const particle &particle);
    bool pass_angle(const particle &particle);
};

struct lana
{
};

struct hira
{
    std::map<std::string, std::array<double, 2>> Ekinlabcut;
    std::array<double, 2> thetalabcut = {30., 75.};
    std::array<double, 2> philabcut = {0, 360};
    void init();
    bool pass_ekinlab(const particle &particle);
    bool pass_angle(const particle &particle);
};

struct uball
{
    fs::path e15190_uball_config = PROJECT_DIR / "database/microball/acceptance/e15190_config.dat";
    detector_uball *uball;
    int counter = 0;
    void init();
    void reset();
    void config(const std::string &system);
    void read_particle(const particle &particle);
};

void hira::init()
{
    this->Ekinlabcut["p"] = {20.0, 198.0};
    this->Ekinlabcut["d"] = {15.0, 263.0 / 2};
    this->Ekinlabcut["t"] = {12.0, 312 / 3.};
    this->Ekinlabcut["3He"] = {20.0, 200.0};
    this->Ekinlabcut["4He"] = {18.0, 200.0};
}

bool hira::pass_ekinlab(const particle &particle)
{
    if (this->Ekinlabcut.count(particle.name) == 0)
    {
        return 0;
    }
    return (particle.ekinlab >= this->Ekinlabcut[particle.name][0] && particle.ekinlab <= this->Ekinlabcut[particle.name][1]);
}
bool hira::pass_angle(const particle &particle)
{
    return (particle.thetalab >= this->thetalabcut[0] && particle.thetalab <= this->thetalabcut[1] && particle.phi >= this->philabcut[0] && particle.phi <= this->philabcut[1]);
}

void uball::init()
{
    fs::path path_base = PROJECT_DIR / "database/microball/acceptance";
    fs::path path_geo = path_base / "uball_geometry.dat";
    std::map<int, fs::path> path_thres;
    path_thres[2] = path_base / "threshold_Sn_65mgcm2.dat";
    path_thres[3] = path_base / "threshold_Sn_58mgcm2.dat";
    path_thres[4] = path_base / "threshold_Sn_50mgcm2.dat";
    path_thres[5] = path_base / "threshold_Sn_43mgcm2.dat";
    path_thres[7] = path_base / "threshold_Sn_30mgcm2.dat";
    path_thres[8] = path_base / "threshold_Sn_23mgcm2.dat";

    this->uball = new detector_uball();
    this->uball->read_geometry(fs::absolute(path_geo));

    for (auto &[id, pth] : path_thres)
    {
        if (fs ::exists(pth))
        {
            this->uball->read_threshold(id, fs::absolute(pth));
        }
        else
        {
            std::cout << fs::absolute(pth) << " does not exist." << std::endl;
        }
    }

    this->uball->hira_coordinate(); // add phi angle by 90, and move to phi=[0,360];
}

void uball::config(const std::string &system)
{
    std::ifstream stream(fs::absolute(this->e15190_uball_config));
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

void uball::reset()
{
    this->uball->reset_csi();
    this->counter = 0;
}

void uball::read_particle(const particle &particle)
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
