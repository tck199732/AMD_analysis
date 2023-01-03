#include "TFile.h"
#include "TH2D.h"

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "../src/particle.cpp"
#include "../src/root_io.cpp"
#include "../src/system_info.cpp"

struct histograms
{
    std::string reaction;
    double betacms, beam_rapidity;

    std::vector<std::string> particlenames = {
        "n", "p", "d", "t", "3He", "4He"};

    double norm = 0.0;
    std::map<std::string, TH2D *> h2_time_momentum;
    std::map<std::string, TH2D *> h2_time_position;
    void init();
    void fill(const particle &particle, const double &weight);
    void normalize();
    void write();
};

void histograms::init()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_time_momentum[pn] = new TH2D(("h2_time_momentum_" + pn).c_str(), "", 600, 0., 600., 500, 0, 500);
        this->h2_time_position[pn] = new TH2D(("h2_time_position_" + pn).c_str(), "", 400, 0., 40., 500, 0, 500);

        this->h2_time_position[pn]->Sumw2();
        this->h2_time_position[pn]->SetDirectory(0);
        this->h2_time_momentum[pn]->Sumw2();
        this->h2_time_momentum[pn]->SetDirectory(0);
    }
}

void histograms::write()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_time_momentum[pn]->Write();
        this->h2_time_position[pn]->Write();
    }
}

void histograms::normalize()
{
    std::cout << "normalizing histograms : " << this->norm << std::endl;
    for (const auto &pn : this->particlenames)
    {
        this->h2_time_momentum[pn]->Scale(1. / this->norm);
        this->h2_time_position[pn]->Scale(1. / this->norm);
    }
}

void histograms::fill(const particle &particle, const double &weight)
{
    if (this->h2_time_momentum.count(particle.name) == 1)
    {
        this->h2_time_momentum[particle.name]->Fill(particle.pcms, particle.t, 1.);
        this->h2_time_position[particle.name]->Fill(particle.r, particle.t, 1.);
    }
    return;
}