#include "BaseHistograms.hh"

KinergyTheta::KinergyTheta(const std::string &suffix, const std::string &frame) : BaseHistograms(suffix)
{
    this->frame = frame;
    this->name = Form("h2_KinergyTheta_%s_%s", frame.c_str(), suffix.c_str());
    for (auto &pn : this->PARTICLENAMES)
    {
        std::string hname = Form("h2_KinergyTheta_%s_%s_%s", frame.c_str(), suffix.c_str(), pn.c_str());
        this->Histogram2D_Collection[pn] = new TH2D(hname.c_str(), "", 180, 0, 180., 400, 0, 400);
        this->Histogram2D_Collection[pn]->Sumw2();
    }
}

void KinergyTheta::Fill(const Particle &particle, const double &weight)
{
    std::string name = this->NUCLEINAMES[{particle.Z, particle.A}];

    double theta_deg, kinergy;
    if (this->frame == "cms")
    {
        theta_deg = particle.theta_cms * TMath::RadToDeg();
        kinergy = particle.kinergy_cms / particle.A;
    }
    else if (this->frame == "lab")
    {
        theta_deg = particle.theta_lab * TMath::RadToDeg();
        kinergy = particle.kinergy_lab / particle.A;
    }
    else
    {
        std::cerr << "KinergyTheta::Fill: frame " << this->frame << " not supported." << std::endl;
        exit(1);
    }

    if (this->Histogram2D_Collection.count(name) == 1)
    {
        this->Histogram2D_Collection[name]->Fill(theta_deg, kinergy, weight);
    }

    // fill coalescence
    if (this->Histogram2D_Collection.count("coal_p") == 1)
    {
        this->Histogram2D_Collection["coal_p"]->Fill(theta_deg, kinergy, weight * particle.Z);
    }
    if (Histogram2D_Collection.count("coal_n") == 1)
    {
        this->Histogram2D_Collection["coal_n"]->Fill(theta_deg, kinergy, weight * particle.N);
    }
    return;
}
