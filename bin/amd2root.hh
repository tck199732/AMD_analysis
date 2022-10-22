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

struct particle
{
    int Z, N;
    double px, py, pz;
    double t, x, y, z;
    void report();
};

struct event
{
    int eventID;
    double b;
    std::vector<particle> particles;
    void report();
};

struct manager
{
    std::string reaction = "";
    std::string mode = ""; // 21 = table21, 3 = table3
    fs::path path_data, path_out;
    fs::path path_amdgid, path_coll_hist;

    std::string tr_name = "AMD";
    void init();
    void init_writer();
    void compile();
    void compile21();
    void compile21t();
    void compile3();
    void finish();
    void fill(const event &event);
    int amass;
    RootWriter *writer;
};