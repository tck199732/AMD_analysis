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

struct histograms
{
    std::string reaction, mode;
    double norm = 0.;
    TH2D *h2_multi_b;
    void init();
    void fill(const event &event);
    void normalize();
    void write();
};

void histograms::init()
{
    this->h2_multi_b = new TH2D(
        Form("h2_multi_b_mode%s", this->mode.c_str()), "", 80, -0.5, 80.5, 1000, 0., 10.);
    this->h2_multi_b->SetDirectory(0);
}
void histograms::write()
{
    this->h2_multi_b->Write();
}

void histograms::fill(const event &event)
{
    double weight = 1.;
    if (this->mode == "3")
    {
        weight = 1. / ndecays;
    }
    this->h2_multi_b->Fill(event.Nc, event.bimp, weight);
    this->norm += 1.;
}

void histograms::normalize()
{
    if (this->mode == "3")
    {
        this->norm /= ndecays;
    }
    this->h2_multi_b->Scale(1. / this->norm);
}
struct manager
{
    std::string reaction, mode;
    fs::path path_data;
    fs::path path_out;

    system_info *sys_info;
    double betacms, rapidity_beam;

    RootReader *reader;
    histograms hist;
    TFile *outputfile;

    void init();
    void read();
    void fill(const event &event);
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string mode = argv[2];
    std::string path_data = argv[3];
    std::string path_out = argv[4];
    manager manager = {reaction, mode, path_data, path_out};
    manager.init();
    manager.read();
    manager.finish();
}

void manager::init()
{
    if (!fs::exists(this->path_data))
    {
        throw std::invalid_argument("path_data does not exist.");
    }

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

    this->reader = new RootReader("AMD");
    reader->add_file(fs::absolute(this->path_data));
    reader->set_branches(branches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->hist = {this->reaction, this->mode};
    this->hist.init();
}

void manager::read()
{
    int nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;

    int Nc;
    int multi;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    std::vector<particle> particles;
    for (int ievt = 0; ievt < nevents; ievt++)
    {
        std::map<std::string, std::any> map = this->reader->get_entry(ievt);
        try
        {
            Nc = std::any_cast<int>(map["Nc"]);
            multi = std::any_cast<int>(map["multi"]);
            bimp = std::any_cast<double>(map["b"]);
            fn = std::any_cast<int *>(map["N"]);
            fz = std::any_cast<int *>(map["Z"]);
            px = std::any_cast<double *>(map["px"]);
            py = std::any_cast<double *>(map["py"]);
            pz = std::any_cast<double *>(map["pz"]);
        }

        catch (const std::bad_any_cast &e)
        {
            std::cout << e.what() << '\n';
        }

        for (unsigned int i = 0; i < multi; i++)
        {
            particle particle{fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            particles.push_back(particle);
        }
        event event = {Nc, bimp, particles};

        this->fill(event);
        particles.clear();
    }
}

void manager::fill(const event &event)
{
    if (event.Nc == 0)
    {
        return;
    }
    this->hist.fill(event);
}

void manager::finish()
{
    this->outputfile = new TFile(this->path_out.c_str(), "RECREATE");
    this->hist.normalize();
    this->hist.write();
    this->outputfile->Write();
    std::cout << "DONE" << std::endl;
}