#include "../src/Physics.cpp"

#include <array>
#include <vector>
#include <string>
#include <map>
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

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &input_pths, const std::string &analysis = "filtered", const std::string &mode = "3")
{
    for (auto &pth : input_pths)
    {
        if (!fs::exists(pth))
        {
            throw std::invalid_argument("path_data does not exist.");
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

class ArgumentParser
{
    // check this out : https://www.gnu.org/software/libc/manual/html_node/Using-Getopt.html#Using-Getopt
public:
    // required arguments
    std::string reaction = "";
    std::string output_file = "";
    std::vector<std::string> input_files_table21;
    std::vector<std::string> input_files_table3;
    int cut_on_uball_charged_particles = 0;
    int verbose = 0;
    ArgumentParser(int argc, char *argv[])
    {
        // If the value of this variable is nonzero, then getopt prints an error message to the standard error stream if it encounters an unknown option character or an option with a missing required argument. This is the default behavior. If you set this variable to zero, getopt does not print any messages, but it still returns the character ? to indicate an error
        // opterr = 0;

        // When getopt encounters an unknown option character or an option with a missing required argument, it stores that option character in this variable. You can use this for providing your own diagnostic messages.
        // optopt;
        int opt;
        while ((opt = getopt(argc, argv, "hvr:f:p:o:c:")) != -1)
        {
            switch (opt)
            {

            case 'r':
            {
                this->reaction = optarg;
                break;
            }
            case 'o':
            {
                this->output_file = optarg;
                break;
            }
            case 'c':
            {
                this->cut_on_uball_charged_particles = std::stoi(optarg);
                break;
            }
            case 'f':
            {
                std::string files = optarg;
                std::string delimiter = " ";

                size_t pos = 0;
                std::string token;
                while ((pos = files.find(delimiter)) != std::string::npos)
                {
                    token = files.substr(0, pos);
                    std::cout << token << std::endl;
                    this->input_files_table3.push_back(token);
                    files.erase(0, pos + delimiter.length());
                }
                if (files.length() > 0)
                {
                    this->input_files_table3.push_back(files);
                }

                break;
            }
            case 'p':
            {
                std::string files = optarg;
                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;
                while ((pos = files.find(delimiter)) != std::string::npos)
                {
                    token = files.substr(0, pos);
                    this->input_files_table21.push_back(token);
                    files.erase(0, pos + delimiter.length());
                }
                if (files.length() > 0)
                {
                    this->input_files_table21.push_back(files);
                }
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
            case 'v':
            {
                this->verbose = 1;
                break;
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
};
