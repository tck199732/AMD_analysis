#include "../src/system_info.cpp"
#include "../src/particle.cpp"

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

struct data_structure
{
    const static int NMAX = 128;
    int multi, Nc;
    double b;
    std::array<double, NMAX> px;
    std::array<double, NMAX> py;
    std::array<double, NMAX> pz;
    std::array<int, NMAX> N;
    std::array<int, NMAX> Z;

    std::array<double, NMAX> t;
    std::array<double, NMAX> x;
    std::array<double, NMAX> y;
    std::array<double, NMAX> z;
};

data_structure DATASTRUCT;

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &input_pths, const std::string &analysis = "default", const std::string &mode = "3")
{

    for (auto &pth : input_pths)
    {
        if (!fs::exists(pth))
        {
            throw std::invalid_argument("path_data does not exist.");
        }
        chain->Add(pth.c_str());
    }
    if (analysis == "filtered")
    {
        chain->SetBranchAddress("Nc", &DATASTRUCT.Nc);
    }
    chain->SetBranchAddress("multi", &DATASTRUCT.multi);
    chain->SetBranchAddress("b", &DATASTRUCT.b);
    chain->SetBranchAddress("px", &DATASTRUCT.px[0]);
    chain->SetBranchAddress("py", &DATASTRUCT.py[0]);
    chain->SetBranchAddress("pz", &DATASTRUCT.pz[0]);

    chain->SetBranchAddress("N", &DATASTRUCT.N[0]);
    chain->SetBranchAddress("Z", &DATASTRUCT.Z[0]);
    if (mode == "21t")
    {
        chain->SetBranchAddress("t", &DATASTRUCT.t[0]);
        chain->SetBranchAddress("x", &DATASTRUCT.x[0]);
        chain->SetBranchAddress("y", &DATASTRUCT.y[0]);
        chain->SetBranchAddress("z", &DATASTRUCT.z[0]);
    }
}
