#include "Histograms.hh"
Histograms::Histograms(const std::string &suffix)
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        std::string hname = Form("h2_pt_rapidity_%s_%s", suffix.c_str(), pn.c_str());
        this->h2_pta_rapidity_lab[pn] = new TH2D(hname.c_str(), "", 300, -1.5, 1.5, 800, 0, 800);
        this->h2_pta_rapidity_lab[pn]->Sumw2();
        this->h2_pta_rapidity_lab[pn]->SetDirectory(0);
    }
}

void Histograms::Fill(const Particle &particle, const double &weight)
{
    std::string name = this->NUCLEINAMES[{particle.Z, particle.Z + particle.N}];

    double A = particle.Z + particle.N;
    if (this->h2_pta_rapidity_lab.count(name) == 1)
    {
        this->h2_pta_rapidity_lab[name]->Fill(particle.rapidity_lab_normed, particle.pmag_trans / A, weight);
    }

    // fill coalescence
    if (this->h2_pta_rapidity_lab.count("coal_p") == 1)
    {
        this->h2_pta_rapidity_lab["coal_p"]->Fill(particle.rapidity_lab_normed, particle.pmag_trans / A, weight * particle.Z);
    }
    if (h2_pta_rapidity_lab.count("coal_n") == 1)
    {
        this->h2_pta_rapidity_lab["coal_n"]->Fill(particle.rapidity_lab_normed, particle.pmag_trans / A, weight * particle.N);
    }
    return;
}

void Histograms::Write()
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        this->h2_pta_rapidity_lab[pn]->Write();
    }
}

void Histograms::Normalize(const double &scale)
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        this->h2_pta_rapidity_lab[pn]->Scale(1. / scale);
    }
}
