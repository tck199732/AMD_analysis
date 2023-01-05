#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <filesystem>
namespace fs = std::filesystem;

#include "TFile.h"
#include "TTree.h"
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

    // table21
    std::array<double, MAX_MULTI> ENG;
    std::array<double, MAX_MULTI> LANG;
    std::array<double, MAX_MULTI> JX;
    std::array<double, MAX_MULTI> JY;
    std::array<double, MAX_MULTI> JZ;

    // table3
    std::array<double, MAX_MULTI> J;
    std::array<double, MAX_MULTI> M;
    std::array<double, MAX_MULTI> WEIGHT;
    std::array<int, MAX_MULTI> iFRG;

    // table21t
    std::array<double, MAX_MULTI> t;
    std::array<double, MAX_MULTI> x;
    std::array<double, MAX_MULTI> y;
    std::array<double, MAX_MULTI> z;
};

AMD amd;

void Initialize_Tree(TTree *&tree, const std::string &mode)
{
    // Set base branches
    tree->Branch("multi", &amd.multi, "multi/I");
    tree->Branch("b", &amd.b, "b/D");
    tree->Branch("N", &amd.N[0], "N[multi]/I");
    tree->Branch("Z", &amd.Z[0], "Z[multi]/I");
    tree->Branch("px", &amd.px[0], "px[multi]/D");
    tree->Branch("py", &amd.py[0], "py[multi]/D");
    tree->Branch("pz", &amd.pz[0], "pz[multi]/D");

    if (mode == "21" || mode == "21t")
    {
        tree->Branch("ENG", &amd.ENG[0], "ENG[multi]/D");
        tree->Branch("LANG", &amd.LANG[0], "LANG[multi]/D");
        tree->Branch("JX", &amd.JX[0], "JX[multi]/D");
        tree->Branch("JY", &amd.JY[0], "JY[multi]/D");
        tree->Branch("JZ", &amd.JZ[0], "JZ[multi]/D");
    }

    if (mode == "3")
    {
        tree->Branch("J", &amd.J[0], "J[multi]/D");
        tree->Branch("M", &amd.M[0], "M[multi]/D");
        tree->Branch("WEIGHT", &amd.WEIGHT[0], "WEIGHT[multi]/D");
        tree->Branch("iFRG", &amd.iFRG[0], "iFRG[multi]/I");
    }

    if (mode == "21t")
    {
        tree->Branch("t", &amd.t[0], "t[multi]/D");
        tree->Branch("x", &amd.x[0], "x[multi]/D");
        tree->Branch("y", &amd.y[0], "y[multi]/D");
        tree->Branch("z", &amd.z[0], "z[multi]/D");
    }
    return;
}
