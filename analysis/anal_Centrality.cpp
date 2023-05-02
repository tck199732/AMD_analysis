#include "anal.hh"

void Replace_Errorbars(EmissionTime *&hist, EmissionTime *&hist_one_decay);
void analyze_table3(TChain *&chain, EmissionTime *&hist, const ArgumentParser &argparser);
void analyze_table21(TChain *&chain, EmissionTime *&hist, const ArgumentParser &argparser);

int main(int argc, char *argv[])
{
    ArgumentParser argparser(argc, argv);

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, argparser.input_files, argparser.mode, argparser.table);

    ImpactParameterMultiplicity *hist = new ImpactParameterMultiplicity(
        "table" + argparser.table);

    if (argparser.table == "3")
    {
        analyze_table3(chain, hist, argparser);
    }
    else if (argparser.table == "21")
    {
        analyze_table21(chain, hist, argparser);
    }

    // saving results
    TFile *outputfile = new TFile(argparser.output_file.c_str(), "RECREATE");
    outputfile->cd();
    hist->Write();
    outputfile->Write();
    outputfile->Close();
}

void analyze_table3(TChain *&chain, EmissionTime *&hist, const ArgumentParser &argparser)
{
    EmissionTime *hist_one_decay = new EmissionTime("table3_one_decay");
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        int multi = amd.Nc;
        if (argparser.mode == "raw")
        {
            multi = amd.multi;
        }
        if (multi == 0)
        {
            continue;
        }

        // event cut
        if (!(amd.b > argparser.cut_on_impact_parameter[0] && amd.b < argparser.cut_on_impact_parameter[1] && multi > argparser.cut_on_multiplicity[0] && multi < argparser.cut_on_multiplicity[1]))
        {
            continue;
        }

        hist->Fill(multi, amd.b, 1. / NDECAYS);

        if (ievt < chain->GetEntries() / NDECAYS)
        {
            hist_one_decay->Fill(multi, amd.b, 1.);
        }
    }
    hist_one_decay->Normalize(chain->GetEntries() / NDECAYS);
    hist->Normalize(chain->GetEntries() / NDECAYS);
    Replace_Errorbars(hist, hist_one_decay);
    return;
}

void analyze_table21(TChain *&chain, EmissionTime *&hist, const ArgumentParser &argparser)
{
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);
        int multi = amd.multi;

        // event cut
        if (!(amd.b > argparser.cut_on_impact_parameter[0] && amd.b < argparser.cut_on_impact_parameter[1] && multi > argparser.cut_on_multiplicity[0] && multi < argparser.cut_on_multiplicity[1]))
        {
            continue;
        }
        if (multi == 0)
        {
            continue;
        }
        hist->Fill(multi, amd.b);
    }
    hist->Normalize(chain->GetEntries());
    return;
}

void Replace_Errorbars(EmissionTime *&hist, EmissionTime *&hist_one_decay)
{
    TH2D *h2 = hist->Histogram2D_Collection[hist->name];
    TH2D *htemp = hist_one_decay->Histogram2D_Collection[hist_one_decay->name];

    for (int i = 1; i <= h2->GetNbinsX(); i++)
    {
        for (int j = 1; j <= h2->GetNbinsY(); j++)
        {
            double error = htemp->GetBinError(i, j);
            h2->SetBinError(i, j, error);
        }
    }
}
