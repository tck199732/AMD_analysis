#include "../Microball.cpp"
#include "../Physics.cpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <getopt.h>

#include "TChain.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"

struct AMD
{
    const static int MAX_MULTI = 128;
    // base
    int multi;
    double b;
    std::array<int, MAX_MULTI> N;
    std::array<int, MAX_MULTI> Z;
    std::array<double, MAX_MULTI> px;
    std::array<double, MAX_MULTI> py;
    std::array<double, MAX_MULTI> pz;
};

AMD amd;
std::map<std::array<int, 2>, double> ame_mass_table;

void Initialize_TChain(TChain *&chain)
{
    chain->SetBranchAddress("multi", &amd.multi);
    chain->SetBranchAddress("b", &amd.b);
    chain->SetBranchAddress("N", &amd.N[0]);
    chain->SetBranchAddress("Z", &amd.Z[0]);
    chain->SetBranchAddress("px", &amd.px[0]);
    chain->SetBranchAddress("py", &amd.py[0]);
    chain->SetBranchAddress("pz", &amd.pz[0]);

    chain->SetMakeClass(1);
    chain->SetBranchStatus("*", false);
    chain->SetBranchStatus("multi", true);
    chain->SetBranchStatus("b", true);
    chain->SetBranchStatus("N", true);
    chain->SetBranchStatus("Z", true);
    chain->SetBranchStatus("px", true);
    chain->SetBranchStatus("py", true);
    chain->SetBranchStatus("pz", true);
    return;
}
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

class ArgumentParser
{
public:
    // required arguments
    std::string reaction;
    std::string input_file;
    std::string output_file;

    bool is_apply_cut_charged_particle = false;
    bool is_apply_cut_coverage = false;
    bool is_apply_cut_kinergy = false;
    bool is_apply_cut_multiple_hit = false;

    ArgumentParser(int argc, char *argv[])
    {
        options = {
            {"help", no_argument, 0, 'h'},
            {"reaction", required_argument, 0, 'r'},
            {"input", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {"is_apply_cut_charged_particle", no_argument, 0, 'n'},
            {"is_apply_cut_coverage", no_argument, 0, 'f'},
            {"is_apply_cut_kinergy", no_argument, 0, 't'},
            {"is_apply_cut_multiple_hit", no_argument, 0, 'm'},
            {0, 0, 0, 0},
        };

        int option_index = 0;
        int opt;
        while ((opt = getopt_long(argc, argv, "hr:i:o:n:f:t:m", options.data(), &option_index)) != -1)
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
                this->input_file = optarg;
                break;
            }
            case 'o':
            {
                this->output_file = optarg;
                break;
            }
            case 'n':
            {
                this->is_apply_cut_charged_particle = true;
                break;
            }
            case 'f':
            {
                this->is_apply_cut_coverage = true;
                break;
            }
            case 't':
            {
                this->is_apply_cut_kinergy = true;
                break;
            }
            case 'm':
            {
                this->is_apply_cut_multiple_hit = true;
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
