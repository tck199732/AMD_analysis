#include "BaseHistograms.hh"

PmagEmissionTime::PmagEmissionTime(const std::string &suffix) : BaseHistograms(suffix)
{
    this->name = "h2_PmagEmissionTime_" + suffix;
    for (auto &pn : this->PARTICLENAMES)
    {
        std::string hname = name + "_" + pn;

        this->Histogram2D_Collection[pn] = new TH2D(hname.c_str(), "", 800, 0, 800., 500, 0, 500);
        this->Histogram2D_Collection[pn]->Sumw2();
    }
}

void PmagEmissionTime::Fill(const Particle &particle, const double &weight)
{
    std::string pn = this->NUCLEINAMES[{particle.Z, particle.A}];
    if (this->Histogram2D_Collection.count(pn) == 1)
    {
        this->Histogram2D_Collection[pn]->Fill(particle.pmag_cms / particle.A, particle.t_cms, weight);
    }
    return;
}
