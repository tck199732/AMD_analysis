#include "AME.hh"
#include "HiRA.hh"
#include "Particle.hh"
#include "Physics.hh"
#include "Microball.hh"
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

class ArgumentParser
{
public:
    // required arguments
    std::string reaction;
    std::vector<std::string> input_files;
    std::string output_file;

    std::string beam;
    std::string target;
    int beamA, beamZ, targetA, targetZ, beam_energy;

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
                this->_Initialize_ReactionSystem(this->reaction);
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
                std::exit(1);
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