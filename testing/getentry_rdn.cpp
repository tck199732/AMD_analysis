#include <chrono>
#include <algorithm>
#include <random>

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "../src/particle.cpp"
#include "../src/root_io.cpp"
#include "../src/system_info.cpp"

struct manager
{
    int max_event;
    std::string reaction;
    fs::path path_data3;
    RootReader *reader3;
    int nevents3;
    void init();
    void run();
};

int main()
{
    int nevents = 100;
    manager manager = {nevents, "Ca48Ni64E140", "/data/amd/dec2021/b3fm/filtered/Ca48Ni64E140_SkM_table3.root"};
    manager.init();
    manager.run();
}

void manager::init()
{
    std::vector<branch> branches = {
        {"Nc", "int"},
        {"multi", "int"},
        {"b", "double"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
    };

    for (auto &br : branches)
    {
        br.autofill();
    }
    this->reader3 = new RootReader("AMD");
    reader3->add_file(fs::absolute(this->path_data3));
    reader3->set_branches(branches);
    this->nevents3 = this->reader3->tree->GetEntries();
}

void manager::run()
{
    int Nc, multi;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    std::vector<int> eid(this->nevents3);
    for (int ievt3 = 0; ievt3 < nevents3; ievt3++)
    {
        eid[ievt3] = ievt3;
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(eid.begin(), eid.end(), std::default_random_engine(seed));

    for (int ievt3 = 0; ievt3 < this->max_event; ievt3++)
    {
        // get entry in ascending order
        auto start0 = std::chrono::steady_clock::now();
        std::map<std::string, std::any> map0 = this->reader3->get_entry(ievt3);
        auto end0 = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds0 = end0 - start0;

        // get entry in random order
        int ievt = eid[ievt3];
        auto start = std::chrono::steady_clock::now();
        std::map<std::string, std::any> map = this->reader3->get_entry(ievt);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        if (ievt3 % 10 == 0)
        {
            std::cout << ievt3 << "\t" << ievt << "\telapsed time (ascending vs random): " << elapsed_seconds0.count() << "\t" << elapsed_seconds.count() << std::endl;
        }
    }
}
