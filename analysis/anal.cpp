#include "anal.hh"

class Histograms
{
public:
    Histograms(const std::string &suffix);
    ~Histograms() { ; }
    void fill(const particle &particle, const double &weight);
    void normalize(const double &scale);
    void write();
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;

protected:
    std::vector<std::string> PARTICLENAMES = {
        "n", "p", "d", "t", "3He", "4He", "coal_n", "coal_p"};
};

void Replace_Errorbars(Histograms *&hist, Histograms *&hist_one_decay)
{
    for (auto &[pn, h2] : hist->h2_pta_rapidity_lab)
    {
        for (int i = 0; i < h2->GetNbinsX(); i++)
        {
            for (int j = 0; j < h2->GetNbinsY(); j++)
            {
                double error = hist_one_decay->h2_pta_rapidity_lab[pn]->GetBinError(i + 1, j + 1);
                h2->SetBinError(i + 1, j + 1, error);
            }
        }
    }
    return;
}

void analyze_table3(TChain *&chain, Histograms *&hist, const ArgumentParser &argparser);
void analyze_table21(TChain *&chain, Histograms *&hist, const ArgumentParser &argparser);

int main(int argc, char *argv[])
{
    ArgumentParser argparser(argc, argv);
    ReadAMETable(AME_MASS_TABLE);

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, argparser.input_files, argparser.mode, argparser.table);

    Histograms *hist = 0;
    if (argparser.table == "3")
    {
        hist = new Histograms("secondary");
        analyze_table3(chain, hist, argparser);
    }
    else if (argparser.table == "21")
    {
        hist = new Histograms("primary");
        analyze_table21(chain, hist, argparser);
    }

    // saving results
    TFile *outputfile = new TFile(argparser.output_file.c_str(), "RECREATE");
    outputfile->cd();
    hist->write();
    outputfile->Write();
}

void analyze_table3(TChain *&chain, Histograms *&hist, const ArgumentParser &argparser)
{
    Histograms *hist_one_decay = new Histograms("secondary_one_decay");
    double betacms = Physics::GetReactionBeta(argparser.reaction);
    double beam_rapidity = Physics::GetBeamRapidity(argparser.reaction);

    int nevents = chain->GetEntries();
    double norm = 0.;
    double norm_one_decay = 0.;

    for (int ievt = 0; ievt < nevents; ievt++)
    {
        chain->GetEntry(ievt);
        int multi;
        if (argparser.mode == "filtered")
            multi = amd.Nc;
        else if (argparser.mode == "raw")
            multi = amd.multi;

        if (multi >= argparser.cut_on_multiplicity[0] && multi <= argparser.cut_on_multiplicity[1] && amd.b >= argparser.cut_on_impact_parameter[0] && amd.b <= argparser.cut_on_impact_parameter[1])
        {
            norm += 1. / NDECAYS;
            if (ievt < nevents / NDECAYS)
            {
                norm_one_decay += 1.;
            }
        }
        else
        {
            continue;
        }

        for (int i = 0; i < amd.multi; i++)
        {
            particle particle = {amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i]};
            particle.initialize(betacms, beam_rapidity);
            if (ievt < nevents / NDECAYS)
            {
                hist_one_decay->fill(particle, 1.);
            }
            hist->fill(particle, 1. / NDECAYS);
        }
    }
    hist->normalize(norm);
    hist_one_decay->normalize(norm_one_decay);
    Replace_Errorbars(hist, hist_one_decay);
}

void analyze_table21(TChain *&chain, Histograms *&hist, const ArgumentParser &argparser)
{
    double betacms = Physics::GetReactionBeta(argparser.reaction);
    double beam_rapidity = Physics::GetBeamRapidity(argparser.reaction);

    int nevents = chain->GetEntries();
    double norm = 0.;

    for (int ievt = 0; ievt < nevents; ievt++)
    {
        chain->GetEntry(ievt);

        if (amd.multi >= argparser.cut_on_multiplicity[0] && amd.multi <= argparser.cut_on_multiplicity[1] && amd.b >= argparser.cut_on_impact_parameter[0] && amd.b <= argparser.cut_on_impact_parameter[1])
        {
            norm += 1.;
        }
        else
        {
            continue;
        }

        for (int i = 0; i < amd.multi; i++)
        {
            particle particle = {amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i]};
            particle.initialize(betacms, beam_rapidity);
            hist->fill(particle, 1.);
        }
    }
    hist->normalize(norm);
}

Histograms::Histograms(const std::string &suffix)
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        std::string hname = Form("h2_pt_rapidity_%s_%s", suffix.c_str(), pn.c_str());
        this->h2_pta_rapidity_lab[pn] = new TH2D(hname.c_str(), "", 100, 0., 1., 800, 0, 800);
        this->h2_pta_rapidity_lab[pn]->Sumw2();
        this->h2_pta_rapidity_lab[pn]->SetDirectory(0);
    }
}

void Histograms::fill(const particle &particle, const double &weight)
{
    std::string name = Physics::GetNucleiName(particle.Z, particle.Z + particle.N);

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

void Histograms::write()
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        this->h2_pta_rapidity_lab[pn]->Write();
    }
}

void Histograms::normalize(const double &scale)
{
    for (const auto &pn : this->PARTICLENAMES)
    {
        this->h2_pta_rapidity_lab[pn]->Scale(1. / scale);
    }
}