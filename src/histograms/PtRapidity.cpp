#include "BaseHistograms.hh"

PtRapidity::PtRapidity(const std::string &suffix) : BaseHistograms(suffix)
{
    for (auto &pn : this->PARTICLENAMES)
    {
        this->Histogram2D_Collection[pn] = new TH2D(("h2_PtRapidity_" + suffix + "_" + pn).c_str(), "", 100, 0, 1., 800, 0, 800);
        this->Histogram2D_Collection[pn]->Sumw2();
    }
}

void PtRapidity::Fill(const Particle &particle, const double &weight)
{
    std::string name = this->NUCLEINAMES[{particle.Z, particle.A}];

    double pta = particle.pmag_trans / particle.A;
    double normed_rapidity = particle.rapidity_lab_normed;

    if (this->Histogram2D_Collection.count(name) == 1)
    {
        this->Histogram2D_Collection[name]->Fill(normed_rapidity, pta, weight);
    }

    // fill coalescence
    if (this->Histogram2D_Collection.count("coal_p") == 1)
    {
        this->Histogram2D_Collection["coal_p"]->Fill(normed_rapidity, pta, weight * particle.Z);
    }
    if (Histogram2D_Collection.count("coal_n") == 1)
    {
        this->Histogram2D_Collection["coal_n"]->Fill(normed_rapidity, pta, weight * particle.N);
    }
    return;
}
