#include "AME.hh"
#include "Particle.hh"
#include "Physics.hh"
#include "ProgressBar.cpp"
#include "BaseHistograms.hh"

#include <array>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

const int NDECAYS = 10; // 10 decay events generated from each primary event

struct AMD
{
    const static int NMAX = 128;
    int multi, Nc;
    double b;
    std::array<double, NMAX> px;
    std::array<double, NMAX> py;
    std::array<double, NMAX> pz;
    std::array<int, NMAX> N;
    std::array<int, NMAX> Z;
};

AMD amd;
void Initialize_TChain(TChain *&chain, const std::vector<std::string> &input_pths, const std::string &analysis = "filtered", const std::string &mode = "3");

class ArgumentParser
{
    // check this out : https://www.gnu.org/software/libc/manual/html_node/Using-Getopt.html#Using-Getopt
public:
    // required arguments
    std::string reaction = "";
    std::vector<std::string> input_files = {};
    std::string output_file = "";
    std::string mode = "filtered";
    std::string table = "3";
    std::array<int, 2> cut_on_multiplicity = {0, 128};
    std::array<double, 2> cut_on_impact_parameter = {0., 3.};

    std::string beam;
    std::string target;
    int beamA, beamZ, targetA, targetZ, beam_energy;

    ArgumentParser(int argc, char *argv[])
    {

        options = {
            {"help", no_argument, 0, 'h'},
            {"reaction", required_argument, 0, 'r'},
            {"input", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {"mode", required_argument, 0, 'm'},
            {"table", required_argument, 0, 't'},
            {"cut_on_multiplicity", required_argument, 0, 'c'},
            {"cut_on_impact_parameter", required_argument, 0, 'b'},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int opt;
        while ((opt = getopt_long(argc, argv, "hr:i:o:c:b:t:m:", options.data(), &option_index)) != -1)
        {
            switch (opt)
            {
            case 'r':
            {
                this->reaction = optarg;
                this->_Initialize_ReactionSystem(this->reaction);
                break;
            }
            case 'o':
            {
                this->output_file = optarg;
                break;
            }
            case 'c':
            {
                std::string cut_on_multiplicity = optarg;
                std::istringstream iss(cut_on_multiplicity);
                iss >> this->cut_on_multiplicity[0] >> this->cut_on_multiplicity[1];
                break;
            }
            case 'i':
            {
                // split optarg by space
                std::istringstream iss(optarg);
                std::string token;
                while (std::getline(iss, token, ' '))
                {
                    this->input_files.push_back(token);
                }
                break;
            }

            case 'b':
            {
                std::string cut_on_bfm = optarg;
                std::istringstream iss(cut_on_bfm);
                iss >> this->cut_on_impact_parameter[0] >> this->cut_on_impact_parameter[1];
                break;
            }
            case 'm':
            {
                this->mode = optarg;
                break;
            }
            case 't':
            {
                this->table = optarg;
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
                exit(1);
            }

            default:
            {
                std::cout << "Got unknown parse returns: " << optarg << std::endl;
            }
            }
        }

        if (this->table == "21" && this->mode == "filtered")
        {
            std::cout << "Table 21 is not available for filtered mode." << std::endl;
            exit(1);
        }

        for (auto f : this->input_files)
        {
            if (!fs::exists(f))
            {
                std::cerr << "File " << f.c_str() << " does not exist." << std::endl;
                std::exit(1);
            }
        }
    }
    void help()
    {
        const char *msg = R"(
            -r      reaction tag, e.g. `Ca48Ni64E140`
            -i      a list of input ROOT files, separated by space.
            -o      ROOT file output path.
            -c      cut on uball charged particles, e.g. `0 128`
            -b      cut on impact parameter, e.g. `0. 3.`
            -m      mode, either `filtered` or `raw`
            -t      table number, either `21` or `3` (implement 21t later)
            -h      Print help message.
        )";
        std::cout << msg << std::endl;
    }
    void _Initialize_ReactionSystem(const std::string &reaction)
    {
        {
            std::regex pattern("[A-Z][a-z]");
            std::vector<std::string> tokens;
            std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
            std::sregex_iterator end;

            while (iter != end)
            {
                std::smatch match = *iter;
                tokens.push_back(match.str());
                ++iter;
            }

            this->beam = tokens[0];
            this->target = tokens[1];

            auto GetZ = [](const std::string &nuclei) -> double
            {
                if (nuclei == "Ca")
                    return 20;
                else if (nuclei == "Ni")
                    return 28;
                else if (nuclei == "Sn")
                    return 50;
                else
                    return 0;
            };
            this->beamZ = GetZ(tokens[0]);
            this->targetZ = GetZ(tokens[1]);
        }

        {
            std::regex pattern("[0-9]+");
            std::vector<int> tokens;

            std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
            std::sregex_iterator end;

            while (iter != end)
            {
                std::smatch match = *iter;
                tokens.push_back(std::stoi(match.str()));
                ++iter;
            }
            this->beamA = tokens[0];
            this->targetA = tokens[1];
            this->beam_energy = tokens[2];
        }
    }

protected:
    std::vector<option> options;
};

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &input_pths, const std::string &analysis, const std::string &mode)
{
    for (auto &pth : input_pths)
    {
        if (!fs::exists(pth))
        {
            std::string msg = Form("%s does not exist.", pth.c_str());
            throw std::invalid_argument(msg.c_str());
        }
        chain->Add(pth.c_str());
    }

    if (analysis == "filtered" && mode == "3")
    {
        chain->SetBranchAddress("uball_multi", &amd.Nc);
        chain->SetBranchAddress("hira_multi", &amd.multi);
        chain->SetBranchAddress("b", &amd.b);
        chain->SetBranchAddress("hira_px", &amd.px[0]);
        chain->SetBranchAddress("hira_py", &amd.py[0]);
        chain->SetBranchAddress("hira_pz", &amd.pz[0]);
        chain->SetBranchAddress("hira_N", &amd.N[0]);
        chain->SetBranchAddress("hira_Z", &amd.Z[0]);
    }

    else if (analysis == "raw")
    {
        chain->SetBranchAddress("multi", &amd.multi);
        chain->SetBranchAddress("b", &amd.b);
        chain->SetBranchAddress("px", &amd.px[0]);
        chain->SetBranchAddress("py", &amd.py[0]);
        chain->SetBranchAddress("pz", &amd.pz[0]);
        chain->SetBranchAddress("N", &amd.N[0]);
        chain->SetBranchAddress("Z", &amd.Z[0]);
    }
}
