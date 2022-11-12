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
    TH1D *h1_multi_H_rapidity_cut;
    TH1D *h1_multi_He_rapidity_cut;

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

    this->h1_multi_H_rapidity_cut = new TH1D(("h1_multi_H_rapidity_cut_mode" + this->mode).c_str(), "", 10, 0, 10);
    this->h1_multi_He_rapidity_cut = new TH1D(("h1_multi_He_rapidity_cut_mode" + this->mode).c_str(), "", 10, 0, 10);

    this->h1_multi_H->Sumw2();
    this->h1_multi_He->Sumw2();
    this->h1_multi_H->SetDirectory(0);
    this->h1_multi_He->SetDirectory(0);

    this->h1_multi_H_rapidity_cut->Sumw2();
    this->h1_multi_He_rapidity_cut->Sumw2();
    this->h1_multi_H_rapidity_cut->SetDirectory(0);
    this->h1_multi_He_rapidity_cut->SetDirectory(0);
}

void histograms::write()
{
    this->h1_multi_H->Write();
    this->h1_multi_He->Write();
    this->h1_multi_H_rapidity_cut->Write();
    this->h1_multi_He_rapidity_cut->Write();
}

void histograms::normalize()
{
    std::cout << "normalizing histograms : " << this->norm << std::endl;

    this->h1_multi_H_rapidity_cut->Scale(1. / this->norm);
    this->h1_multi_He_rapidity_cut->Scale(1. / this->norm);

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
            if (par.rapidity_lab / this->rapidity_beam >= 0.4 && par.rapidity_lab / this->rapidity_beam <= 0.6)
            {
                this->h1_multi_H_rapidity_cut->Fill(par.aid, this->weight);
            }
        }
        else if (par.zid == 2)
        {
            this->h1_multi_He->Fill(par.aid, this->weight);
            if (par.rapidity_lab / this->rapidity_beam >= 0.4 && par.rapidity_lab / this->rapidity_beam <= 0.6)
            {
                this->h1_multi_He_rapidity_cut->Fill(par.aid, this->weight);
            }
        }
    }
    this->norm += this->weight;
    return;
}