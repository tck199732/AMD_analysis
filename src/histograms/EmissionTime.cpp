#include "BaseHistograms.hh"

EmissionTime::EmissionTime(const std::string &suffix) : BaseHistograms(suffix)
{
    for (auto &pn : this->PARTICLENAMES)
    {
        std::string hname = Form("h2_PmagEmissionTime_%s_%s", suffix.c_str(), pn.c_str());

        this->Histogram2D_Collection[pn] = new TH2D(hname.c_str(), "", 800, 0, 800., 500, 0, 500);
        this->Histogram2D_Collection[pn]->Sumw2();
    }
}

void EmissionTime::Fill(const Particle &particle, const double &weight)
{
    std::string name = this->NUCLEINAMES[{particle.Z, particle.A}];
    if (this->Histogram2D_Collection.count(name) == 1)
    {
        this->Histogram2D_Collection[name]->Fill(particle.pmag_cms / particle.A, particle.t_cms, weight);
    }
    return;
}
