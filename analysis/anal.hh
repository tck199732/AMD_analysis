#include "../src/Physics.cpp"
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

struct AMD_Structure
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

AMD_Structure AMD;

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
        chain->SetBranchAddress("Nc", &AMD.Nc);
    }
    chain->SetBranchAddress("multi", &AMD.multi);
    chain->SetBranchAddress("b", &AMD.b);
    chain->SetBranchAddress("px", &AMD.px[0]);
    chain->SetBranchAddress("py", &AMD.py[0]);
    chain->SetBranchAddress("pz", &AMD.pz[0]);

    chain->SetBranchAddress("N", &AMD.N[0]);
    chain->SetBranchAddress("Z", &AMD.Z[0]);

    if (mode == "21t")
    {
        chain->SetBranchAddress("t", &AMD.t[0]);
        chain->SetBranchAddress("x", &AMD.x[0]);
        chain->SetBranchAddress("y", &AMD.y[0]);
        chain->SetBranchAddress("z", &AMD.z[0]);
    }
}
