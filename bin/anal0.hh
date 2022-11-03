
#include "../src/system_info.cpp"
#include "../src/detector_uball.cpp"
#include "../src/particle.cpp"
#include "../src/root_io.cpp"

#include <vector>
#include <string>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

fs::path PROJECT_DIR = "/home/kin/Desktop/amd_analysis";
const int ndecays = 10;

struct event
{
    int multi;
    double bimp;
    std::vector<particle> particles;
};

// struct histograms21
// {
//     std::string reaction;
//     TH1D *h1_multi;
//     TH1D *h1_multi_H;
//     TH1D *h1_multi_He;
//     TH1D *h1_Z;
//     double weight = 1.;
//     double norm = 0.0;
//     void init();
//     void fill(const event &event);
//     void normalize();
//     void write();
// };

struct histograms
{
    std::string reaction, mode;
    TH1D *h1_multi;
    TH1D *h1_multi_H;
    TH1D *h1_multi_He;
    TH1D *h1_Z_prim;
    TH1D *h1_Z_seq;

    double weight = 1.;
    double norm = 0.0;
    void init();
    void fill(const event &event);
    void fill(const event &event21, const event &event3);
    void normalize();
    void write();
};

void histograms::init()
{
    this->h1_multi = new TH1D(("h1_multi_mode" + this->mode).c_str(), "", 80, 0, 80);
    this->h1_multi_H = new TH1D(("h1_multi_H_mode" + this->mode).c_str(), "", 10, 0, 10);
    this->h1_multi_He = new TH1D(("h1_multi_He_mode" + this->mode).c_str(), "", 10, 0, 10);
    this->h1_Z_prim = new TH1D(("h1_Z_prim_mode" + this->mode).c_str(), "", 50, 0, 50);

    this->h1_multi->Sumw2();
    this->h1_multi_H->Sumw2();
    this->h1_multi_He->Sumw2();
    this->h1_Z_prim->Sumw2();
    this->h1_multi->SetDirectory(0);
    this->h1_multi_H->SetDirectory(0);
    this->h1_multi_He->SetDirectory(0);
    this->h1_Z_prim->SetDirectory(0);

    if (this->mode == "3")
    {
        this->weight = 1. / ndecays;
        this->h1_Z_seq = new TH1D(("h1_Z_seq_mode" + this->mode).c_str(), "", 50, 0, 50);
        this->h1_Z_seq->Sumw2();
        this->h1_Z_seq->SetDirectory(0);
    }
}

void histograms::fill(const event &event)
{
    this->h1_multi->Fill(event.multi, this->weight);

    for (auto &par : event.particles)
    {
        if (par.zid == 1)
        {
            this->h1_multi_H->Fill(par.aid, this->weight);
        }
        if (par.zid == 2)
        {
            this->h1_multi_He->Fill(par.aid, this->weight);
        }
        if (this->mode == "21")
        {
            this->h1_Z_prim->Fill(par.zid, this->weight);
        }
    }

    this->norm += this->weight;
}
void histograms::fill(const event &event21, const event &event3)
{
    if (this->weight != 1. / ndecays)
    {
        throw std::invalid_argument("should not use this fill function");
    }
    this->fill(event3);

    for (auto &par3 : event3.particles)
    {
        bool isPrim = 0;
        for (auto &par21 : event21.particles)
        {
            if (par3.nid == par21.nid && par3.aid == par21.aid)
            {
                if (pow(par3.px - par21.px, 2.) + pow(par3.py - par21.py, 2.) + pow(par3.pz - par21.pz, 2.) <= 0.05)
                {
                    isPrim = 1;
                }
            }
        }
        if (isPrim)
        {
            this->h1_Z_prim->Fill(par3.zid, this->weight);
        }
        else
        {
            this->h1_Z_seq->Fill(par3.zid, this->weight);
        }
    }
}

void histograms::normalize()
{
    std::cout << "mode : " << this->mode << "\tnormalization : " << this->norm << std::endl;
    this->h1_multi->Scale(1. / this->norm);
    this->h1_multi_H->Scale(1. / this->norm);
    this->h1_multi_He->Scale(1. / this->norm);
    this->h1_Z_prim->Scale(1. / this->norm);
    if (this->mode == "3")
    {
        this->h1_Z_seq->Scale(1. / this->norm);
    }
}
void histograms::write()
{
    this->h1_multi->Write();
    this->h1_multi_H->Write();
    this->h1_multi_He->Write();
    this->h1_Z_prim->Write();
    if (this->mode == "3")
    {
        this->h1_Z_seq->Write();
    }
}
