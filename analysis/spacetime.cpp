
#include "anal.hh"
#include <chrono>

struct Histograms
{
    std::vector<std::string> particlenames = {
        "n", "p", "d", "t", "3He", "4He"};

    std::map<std::string, TH2D *> h2_time_momentum;
    std::map<std::string, TH2D *> h2_time_position;
    void init();
    void fill(const particle &particle, const double &weight);
    void normalize(const double &norm);
    void write();
};

void Histograms::init()
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

void Histograms::write()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_time_momentum[pn]->Write();
        this->h2_time_position[pn]->Write();
    }
}

void Histograms::normalize(const double &norm)
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_time_momentum[pn]->Scale(1. / norm);
        this->h2_time_position[pn]->Scale(1. / norm);
    }
}

void Histograms::fill(const particle &particle, const double &weight)
{
    if (this->h2_time_momentum.count(particle.name) == 1)
    {
        this->h2_time_momentum[particle.name]->Fill(particle.pcms, particle.t, 1.);
        this->h2_time_position[particle.name]->Fill(particle.r, particle.t, 1.);
    }
    return;
}

int main(int argc, char *argv[])
{
    // parse args
    std::string reaction = argv[1];
    std::string output_pth = argv[2];

    std::vector<std::string> input_pths;
    for (int i = 3; i < argc; i++)
    {
        input_pths.push_back(argv[i]);
    }

    // initialize TChain and structures
    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, input_pths, "default", "21t");

    Histograms histograms;
    histograms.init();

    double betacms = Physics::GetReactionBeta(reaction);
    double beam_rapidity = Physics::GetBeamRapidity(reaction);

    // count the number particles with emission time <= 0. (tested all are 0 not negative);
    long count_neg_time = 0;

    // main anal
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        for (unsigned int i = 0; i < AMD.multi; i++)
        {
            particle particle = {AMD.N[i], AMD.Z[i], AMD.px[i], AMD.py[i], AMD.pz[i]};
            particle.set_xyzt(AMD.x[i], AMD.y[i], AMD.z[i], AMD.t[i]);
            particle.autofill(betacms);

            if (particle.t <= 0.0)
            {
                count_neg_time += 1;
                continue;
            }
            histograms.fill(particle, 1.);
        }
    }
    histograms.normalize(chain->GetEntries());

    std::cout << "number fo unphysical emission time : " << count_neg_time << std::endl;

    // saving result
    TFile *outputfile = new TFile(output_pth.c_str(), "RECREATE");
    outputfile->cd();
    histograms.write();
    outputfile->Write();
}