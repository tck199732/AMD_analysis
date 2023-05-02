#include "anal.hh"

void Replace_Errorbars(PtRapidity *&hist, PtRapidity *&hist_one_decay);
void analyze_table3(TChain *&chain, PtRapidity *&hist, const ArgumentParser &argparser);
void analyze_table21(TChain *&chain, PtRapidity *&hist, const ArgumentParser &argparser);

AME *ame;
int main(int argc, char *argv[])
{
    ame = new AME();
    ArgumentParser argparser(argc, argv);

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, argparser.input_files, argparser.mode, argparser.table);

    PtRapidity *hist = 0;
    if (argparser.table == "3")
    {
        hist = new PtRapidity("secondary");
        analyze_table3(chain, hist, argparser);
    }
    else if (argparser.table == "21")
    {
        hist = new PtRapidity("primary");
        analyze_table21(chain, hist, argparser);
    }

    // saving results
    TFile *outputfile = new TFile(argparser.output_file.c_str(), "RECREATE");
    outputfile->cd();
    hist->Write();
    outputfile->Write();
}

void analyze_table3(TChain *&chain, PtRapidity *&hist, const ArgumentParser &argparser)
{
    PtRapidity *hist_one_decay = new PtRapidity("secondary_one_decay");
    double beam_mass = ame->GetMass(argparser.beamZ, argparser.beamA);
    double target_mass = ame->GetMass(argparser.targetZ, argparser.targetA);
    double betacms = Physics::GetReactionBeta(beam_mass, target_mass, argparser.beam_energy, argparser.beamA);
    double rapidity_beam = Physics::GetBeamRapidity(beam_mass, target_mass, argparser.beam_energy, argparser.beamA);

    int nevents = chain->GetEntries();
    double norm = 0.;
    double norm_one_decay = 0.;

    ProgressBar bar(nevents, argparser.reaction);
    for (int ievt = 0; ievt < nevents; ievt++)
    {
        chain->GetEntry(ievt);
        int multi;
        std::string frame = "";
        if (argparser.mode == "filtered")
        {
            multi = amd.Nc;
            frame = "lab";
        }

        else if (argparser.mode == "raw")
        {
            multi = amd.multi;
            frame = "cms";
        }

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
            bar.Update();
            continue;
        }

        for (int i = 0; i < amd.multi; i++)
        {
            double A = amd.N[i] + amd.Z[i];
            double mass = ame->GetMass(amd.Z[i], A);

            Particle particle(amd.N[i], amd.Z[i], amd.px[i] / A, amd.py[i] / A, amd.pz[i] / A, mass, frame);
            particle.Initialize(betacms, rapidity_beam);

            if (ievt < nevents / NDECAYS)
            {
                hist_one_decay->Fill(particle, 1.);
            }
            hist->Fill(particle, 1. / NDECAYS);
        }
        bar.Update();
    }
    hist->Normalize(norm);
    hist_one_decay->Normalize(norm_one_decay);
    Replace_Errorbars(hist, hist_one_decay);
}

void analyze_table21(TChain *&chain, PtRapidity *&hist, const ArgumentParser &argparser)
{
    double beam_mass = ame->GetMass(argparser.beamZ, argparser.beamA);
    double target_mass = ame->GetMass(argparser.targetZ, argparser.targetA);
    double betacms = Physics::GetReactionBeta(beam_mass, target_mass, argparser.beam_energy, argparser.beamA);
    double rapidity_beam = Physics::GetBeamRapidity(beam_mass, target_mass, argparser.beam_energy, argparser.beamA);
    int nevents = chain->GetEntries();
    double norm = 0.;

    ProgressBar bar(nevents, argparser.reaction);
    for (int ievt = 0; ievt < nevents; ievt++)
    {
        chain->GetEntry(ievt);
        if (amd.multi >= argparser.cut_on_multiplicity[0] && amd.multi <= argparser.cut_on_multiplicity[1] && amd.b >= argparser.cut_on_impact_parameter[0] && amd.b <= argparser.cut_on_impact_parameter[1])
        {
            norm += 1.;
        }
        else
        {
            bar.Update();
            continue;
        }

        for (int i = 0; i < amd.multi; i++)
        {
            double mass = ame->GetMass(amd.Z[i], amd.Z[i] + amd.N[i]);
            Particle particle(amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i], mass);
            particle.Initialize(betacms, rapidity_beam);
            hist->Fill(particle, 1.);
        }
        bar.Update();
    }
    hist->Normalize(norm);
}

void Replace_Errorbars(PtRapidity *&hist, PtRapidity *&hist_one_decay)
{
    for (auto &[pn, h2] : hist->Histogram2D_Collection)
    {
        for (int i = 0; i < h2->GetNbinsX(); i++)
        {
            for (int j = 0; j < h2->GetNbinsY(); j++)
            {
                double error = hist_one_decay->Histogram2D_Collection[pn]->GetBinError(i + 1, j + 1);
                h2->SetBinError(i + 1, j + 1, error);
            }
        }
    }
    return;
}
