#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <filesystem>
namespace fs = std::filesystem;
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "../src/root_io.cpp"

int main()
{
    std::string mode = "3";
    RootReader *reader = new RootReader("AMD");
    // e.g.
    reader->add_file("/data/amd/dec2021/b3fm/Ca48Ni64_En140MeV_SkM_110k_table3.root");
    std::vector<branch> branches = {
        {"multi", "int"},
        {"b", "double"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
    };
    if (mode == "3")
    {
        branches.push_back({"J", "double[multi]"});
        branches.push_back({"M", "double[multi]"});
        branches.push_back({"WEIGHT", "double[multi]"});
        branches.push_back({"iFRG", "int[multi]"});
    }
    for (auto &br : branches)
    {
        br.autofill();
    }
    reader->set_branches(branches);
    int nevents3 = reader->tree->GetEntries();

    std::ofstream outfile("test3.dat");
    outfile << "          10  MONTE CARLO PARTICLES PER FRAGMENT\n";

    auto write = [&outfile](auto &&...args) -> void
    {
        ((outfile << args << "\t"), ...) << "\n";
    };

    int fmulti;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;
    double *J;
    double *M;
    double *WEIGHT;
    int *iFRAG;

    for (int ievt = 0; ievt < nevents3; ievt++)
    {
        std::map<std::string, std::any> map = reader->get_entry(ievt);
        try
        {
            fmulti = std::any_cast<int>(map["multi"]);
            bimp = std::any_cast<double>(map["b"]);
            fn = std::any_cast<int *>(map["N"]);
            fz = std::any_cast<int *>(map["Z"]);
            px = std::any_cast<double *>(map["px"]);
            py = std::any_cast<double *>(map["py"]);
            pz = std::any_cast<double *>(map["pz"]);
            J = std::any_cast<double *>(map["J"]);
            M = std::any_cast<double *>(map["M"]);
            WEIGHT = std::any_cast<double *>(map["WEIGHT"]);
            iFRAG = std::any_cast<int *>(map["iFRG"]);
        }
        catch (const std::bad_any_cast &e)
        {
            std::cout << e.what() << '\n';
        }

        int sim = ievt % (nevents3 / 10) + 1;
        for (unsigned int i = 0; i < fmulti; i++)
        {
            write(fz[i], fn[i], px[i], py[i], pz[i], J[i], M[i], WEIGHT[i], bimp, sim, iFRAG[i]);
            outfile << "\n";
            // outfile << fz[i] << "\t" << fn[i] << "\t" << px[i] << "\t" << py[i] << "\t" << pz[i] << "\t" << J[i] << "\t" << M[i] << "\t" << WEIGHT[i] << "\t" << bimp << "\t" << sim << "\t" << iFRAG[i] << "\n";
        }
    }
}