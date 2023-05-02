#include "BaseHistograms.hh"
BaseHistograms::BaseHistograms(const std::string &suffix)
{
    this->name = suffix;
}

BaseHistograms::~BaseHistograms()
{
    for (auto &[name, h1] : this->Histogram1D_Collection)
    {
        delete h1;
    }

    for (auto &[name, h2] : this->Histogram2D_Collection)
    {
        delete h2;
    }

    for (auto &[name, h3] : this->Histogram3D_Collection)
    {
        delete h3;
    }
}

void BaseHistograms::Fill(const Particle &particle, const double &weight)
{
    throw std::runtime_error("Histograms::Fill() is not implemented.");
    std::exit(1);
}

void BaseHistograms::Write()
{
    for (auto &[name, h1] : this->Histogram1D_Collection)
    {
        h1->Write();
    }

    for (auto &[name, h2] : this->Histogram2D_Collection)
    {
        h2->Write();
    }

    for (auto &[name, h3] : this->Histogram3D_Collection)
    {
        h3->Write();
    }
}

void BaseHistograms::Normalize(const double &scale)
{
    for (auto &[name, h1] : this->Histogram1D_Collection)
    {
        h1->Scale(1. / scale);
    }

    for (auto &[name, h2] : this->Histogram2D_Collection)
    {
        h2->Scale(1. / scale);
    }

    for (auto &[name, h3] : this->Histogram3D_Collection)
    {
        h3->Scale(1. / scale);
    }
}
