#include "Physics.cpp"
#include "Microball.cpp"
#include "HiRA.cpp"
#include "ProgressBar.cpp"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <getopt.h>
#include <filesystem>
namespace fs = std::filesystem;

std::map<std::array<int, 2>, double> AME_MASS_TABLE;
void ReadAMETable(std::map<std::array<int, 2>, double> &ame_mass_table, const std::string &filename = "")
{
    std::string path = filename;
    if (path.empty())
    {
        fs::path project_dir = std::getenv("PROJECT_DIR");
        fs::path dat_path = project_dir / "database/ame/ame_mass.txt";
        path = fs::absolute(dat_path);
    }
    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');
    int Z, A;
    double mass;
    while (infile >> Z)
    {
        infile >> A >> mass;
        ame_mass_table[{Z, A}] = mass;
    }
}

struct particle
{
    int N, Z;
    double px, py, pz_cms;
    double kinergy_lab, pz_lab;
    double theta_lab, phi;
    void initialize(const double &betacms);
};

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

    this->pz_lab = gamma * (this->pz_cms + betacms * (kinergy_cms + mass));
    double pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz_lab, 2.));

    this->theta_lab = TMath::ATan2(pmag_trans, this->pz_lab) * TMath::RadToDeg();

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
        reaction = "";
        input_files = {};
        output_file = "";

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
                this->input_files.push_back(optarg);
                while (optind < argc && *argv[optind] != '-')
                {
                    this->input_files.push_back(argv[optind++]);
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

        if (this->reaction.empty())
        {
            std::cout << "Reaction tag is required." << std::endl;
            this->help();
        }
        if (this->input_files.empty())
        {
            std::cout << "Input files are required." << std::endl;
            this->help();
        }
        if (this->output_file.empty())
        {
            std::cout << "Output file is required." << std::endl;
            this->help();
        }
        for (auto pth : this->input_files)
        {
            if (!fs::exists(pth))
            {
                std::cout << "Input file " << pth << " does not exist." << std::endl;
                this->help();
            }
            std::cout << "Input file: " << pth << std::endl;
        }
    }
    void help()
    {
        const char *msg = R"(
            -r      reaction tag, e.g. Ca48Ni64E140
            -i      a list of input ROOT files, separated by space.
            -o      ROOT file output path.
            -h      Print help message.
        )";
        std::cout << msg << std::endl;
    }

protected:
    std::vector<option> options;
};