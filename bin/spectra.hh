#include "TFile.h"
#include "TH2D.h"

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "../src/particle.cpp"
#include "../src/root_io.cpp"
#include "../src/system_info.cpp"

struct event
{
    int Nc;
    double bimp;
    std::vector<particle> particles;
};

struct eventcut
{
    std::array<int, 2> Nccut;
    std::array<double, 2> bcut;
    bool pass(const event &event)
    {
        return (event.Nc >= this->Nccut[0] && event.Nc <= this->Nccut[1] && event.bimp >= this->bcut[0] && event.bimp < this->bcut[1]);
    }
};

struct histograms
{
    std::string reaction, tag;
    double betacms, beam_rapidity;

    std::vector<std::string> particlenames = {
        "n", "p", "d", "t", "3He", "4He", "coal_n", "coal_p"};

    double norm = 0.0;
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;
    void init();
    void fill(const event &event, const double &weight);
    void normalize();
    void write();
};

void histograms::init()
{
    for (const auto &pn : this->particlenames)
    {
        std::string hname = "h2_pta_rapidity_lab_" + this->tag + "_" + pn;
        this->h2_pta_rapidity_lab[pn] = new TH2D(hname.c_str(), "", 100, 0., 1., 600, 0, 600);
        this->h2_pta_rapidity_lab[pn]->Sumw2();
        this->h2_pta_rapidity_lab[pn]->SetDirectory(0);
    }
}

void histograms::write()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_pta_rapidity_lab[pn]->Write();
    }
}

void histograms::normalize()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_pta_rapidity_lab[pn]->Scale(1. / this->norm);
    }
}

void histograms::fill(const event &event, const double &weight)
{
    for (auto &par : event.particles)
    {
        if (this->h2_pta_rapidity_lab.count(par.name) == 1)
        {
            this->h2_pta_rapidity_lab[par.name]->Fill(par.rapidity_lab / this->beam_rapidity, par.pt, weight);

            if (h2_pta_rapidity_lab.count("coal_p") == 1)
            {
                this->h2_pta_rapidity_lab["coal_p"]->Fill(par.rapidity_lab / this->beam_rapidity, par.pt, weight * par.zid);
            }
            if (h2_pta_rapidity_lab.count("coal_n") == 1)
            {
                this->h2_pta_rapidity_lab["coal_n"]->Fill(par.rapidity_lab / this->beam_rapidity, par.pt, weight * par.nid);
            }
        }
    }
}