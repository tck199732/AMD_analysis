#include "../src/Physics.cpp"

#include <array>
#include <vector>
#include <string>
#include <map>
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
        chain->SetBranchAddress("Nc", &amd.Nc);
    }
    chain->SetBranchAddress("multi", &amd.multi);
    chain->SetBranchAddress("b", &amd.b);
    chain->SetBranchAddress("px", &amd.px[0]);
    chain->SetBranchAddress("py", &amd.py[0]);
    chain->SetBranchAddress("pz", &amd.pz[0]);

    chain->SetBranchAddress("N", &amd.N[0]);
    chain->SetBranchAddress("Z", &amd.Z[0]);
}
