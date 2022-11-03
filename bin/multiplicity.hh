#include "TFile.h"
#include "TH2D.h"

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "../src/particle.cpp"
#include "../src/root_io.cpp"
#include "../src/system_info.cpp"

const int ndecays = 10;
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
    std::string reaction, mode;
    double betacms, rapidity_beam;
    TH1D *h1_multi_H;
    TH1D *h1_multi_He;

    double norm = 0.;
    double weight = 1.;
    void init();
    void fill(const event &event);
    void normalize();
    void write();
};

void histograms::init()
{
    if (this->mode == "3")
    {
        this->weight = 1. / ndecays;
    }
    this->h1_multi_H = new TH1D(("h1_multi_H_mode" + this->mode).c_str(), "", 10, 0, 10);
    this->h1_multi_He = new TH1D(("h1_multi_He_mode" + this->mode).c_str(), "", 10, 0, 10);
    this->h1_multi_H->Sumw2();
    this->h1_multi_He->Sumw2();
    this->h1_multi_H->SetDirectory(0);
    this->h1_multi_He->SetDirectory(0);
}

void histograms::write()
{
    this->h1_multi_H->Write();
    this->h1_multi_He->Write();
}

void histograms::normalize()
{
    this->h1_multi_H->Scale(1. / this->norm);
    this->h1_multi_He->Scale(1. / this->norm);
}

void histograms::fill(const event &event)
{
    for (const auto &par : event.particles)
    {
        if (par.zid == 1)
        {
            this->h1_multi_H->Fill(par.aid, this->weight);
        }
        else if (par.zid == 2)
        {
            this->h1_multi_He->Fill(par.aid, this->weight);
        }
    }
    this->norm += this->weight;
    return;
}