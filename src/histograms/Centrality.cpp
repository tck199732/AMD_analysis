#include "BaseHistograms.hh"

ImpactParameterMultiplicity::ImpactParameterMultiplicity(const std::string &suffix) : BaseHistograms(suffix)
{
    this->name = Form("h2_ImpactParam_Multi_%s", suffix.c_str());
    this->Histogram2D_Collection[name] = new TH2D(name.c_str(), "", 80, -0.5, 79.5, 100., 0., 10.);
    this->Histogram2D_Collection[name]->Sumw2();
}

void ImpactParameterMultiplicity::Fill(const int &multi, const double &b, const double &weight)
{
    this->Histogram2D_Collection[this->name]->Fill(multi, b, weight);
    return;
}