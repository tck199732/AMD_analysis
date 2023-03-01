#include "Physics.cpp"
#include "Microball.cpp"
#include "HiRA.cpp"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <getopt.h>
#include <filesystem>
namespace fs = std::filesystem;

struct particle
{
    int N, Z;
    double px, py, pz_cms;
    double kinergy_lab, pz_lab;
    double theta_lab, phi;
    void initialize(const double &betacms);
};

std::map<std::array<int, 2>, double> AME_MASS_TABLE;
void ReadAMETable(const std::string &filename, std::map<std::array<int, 2>, double> &ame_mass_table)
{
    std::ifstream infile(filename);

    infile.ignore(99, '\n');
    int Z, A;
    double mass;
    while (infile >> Z)
    {
        infile >> A >> mass;
        ame_mass_table[{Z, A}] = mass;
    }
}

void particle::initialize(const double &betacms)
{
    double mass = AME_MASS_TABLE[{this->Z, this->N + this->Z}];
    double A = this->N + this->Z;
    this->px = this->px * A;
    this->py = this->py * A;
    this->pz_cms = this->pz_cms * A;

    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));
    double pmag_trans = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    double pmag_cms = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz_cms, 2.));

    double kinergy_cms = TMath::Sqrt(pow(pmag_cms, 2.) + pow(mass, 2.)) - mass;

    double pz_lab = gamma * (this->pz_cms + betacms * (kinergy_cms + mass));
    double pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(pz_lab, 2.));

    this->theta_lab = TMath::ATan2(pmag_trans, pz_lab) * TMath::RadToDeg();

    this->kinergy_lab = TMath::Sqrt(pow(pmag_lab, 2.) + pow(mass, 2.)) - mass;
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg(); // from -180 to 180 here, need to convert to uball coordinate for each ring
}

void Initialize_MicroBall(Microball *&microball, const std::string &reaction)
{
    fs::path project_dir = std::getenv("PROJECT_DIR");
    fs::path database_dir = project_dir / "database/e15190/microball/acceptance";
    fs::path path_config = database_dir / "config.dat";
    fs::path path_geometry = database_dir / "geometry.dat";
    fs::path path_threshold = database_dir / "fitted_threshold.dat";

    microball->ConfigurateSetup(reaction, path_config.string());
    microball->ReadGeometryMap(path_geometry.string());
    microball->ReadThresholdKinergyMap(path_threshold.string());
}

bool ReadMicroballParticle(Microball *&mb, const particle &part)
{
    bool pass_charge = mb->IsChargedParticle(part.Z);
    bool pass_coverage = mb->IsCovered(part.theta_lab, part.phi);
    bool pass_threshold = mb->IsAccepted(part.kinergy_lab, part.theta_lab, part.N + part.Z, part.Z);
    return pass_charge && pass_coverage && pass_threshold;
}

bool ReadHiRAParticle(HiRA *&hira, const particle &particle)
{
    return hira->PassAngularCut(particle.theta_lab, particle.phi) && hira->PassCharged(particle.Z) && hira->PassKinergyCut(particle.Z + particle.N, particle.Z, particle.kinergy_lab);
}

void correct_phi_value(particle &part, Microball *&microball)
{
    int ring = microball->GetRingID(part.theta_lab);
    if (ring == -1)
    {
        return;
    }
    double phi_min_in_ring = microball->GetPhiMinInRing(ring);
    double phi_max_in_ring = microball->GetPhiMaxInRing(ring);

    if (part.phi < phi_min_in_ring)
    {
        part.phi += 360;
    }
    if (part.phi > phi_max_in_ring)
    {
        part.phi -= 360;
    }
}

class ArgumentParser
{
public:
    // required arguments
    std::string reaction;
    std::vector<std::string> input_files;
    std::string output_file;

    ArgumentParser(int argc, char *argv[])
    {
        options = {
            {"help", no_argument, 0, 'h'},
            {"reaction", required_argument, 0, 'r'},
            {"input", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int opt;
        while ((opt = getopt_long(argc, argv, "hr:i:o:", options.data(), &option_index)) != -1)
        {
            switch (opt)
            {

            case 'r':
            {
                this->reaction = optarg;
                break;
            }
            case 'i':
            {
                std::string files = optarg;
                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;
                while ((pos = files.find(delimiter)) != std::string::npos)
                {
                    token = files.substr(0, pos);
                    this->input_files.push_back(token);
                    files.erase(0, pos + delimiter.length());
                }
                if (files.length() > 0)
                {
                    this->input_files.push_back(files);
                }
                break;
            }
            case 'o':
            {
                this->output_file = optarg;
                break;
            }
            case '?':
            {
                std::cout << "Got unknown option." << std::endl;
                break;
            }

            case 'h':
            {
                this->help();
            }
            default:
            {
                std::cout << "Got unknown parse returns: " << optarg << std::endl;
            }
            }
        }
    }
    void help()
    {
        std::cout << "" << std::endl;
    }

protected:
    std::vector<option> options;
};
