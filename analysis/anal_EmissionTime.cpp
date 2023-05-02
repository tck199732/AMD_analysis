#include "anal.hh"
AME *ame;
int main(int argc, char *argv[])
{
    ArgumentParser argparser(argc, argv);
    ame = new AME();
    if (argparser.table != "21t" || argparser.mode != "raw")
    {
        std::cout << "This analysis is only for table21t and raw mode." << std::endl;
        std::exit(1);
    }

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, argparser.input_files, argparser.mode, argparser.table);

    PmagEmissionTime *hist = new PmagEmissionTime("table21t");

    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);
        // event cut
        if (!(amd.b >= argparser.cut_on_impact_parameter[0] && amd.b <= argparser.cut_on_impact_parameter[1] && amd.multi >= argparser.cut_on_multiplicity[0] && amd.multi <= argparser.cut_on_multiplicity[1]))
        {
            continue;
        }
        for (int ip = 0; ip < amd.multi; ip++)
        {
            double A = amd.N[ip] + amd.Z[ip];
            double mass = ame->GetMass(amd.Z[ip], A);
            Particle particle(amd.N[ip], amd.Z[ip], amd.px[ip], amd.py[ip], amd.pz[ip], mass);
            particle.SetXYZT(amd.x[ip], amd.y[ip], amd.z[ip], amd.t[ip], "cms");
            particle.Initialize(amd.betacms, amd.rapidity_beam);
            hist->Fill(particle);
        }
    }

    // saving results
    TFile *outputfile = new TFile(argparser.output_file.c_str(), "RECREATE");
    outputfile->cd();
    hist->Normalize(chain->GetEntries());
    hist->Write();
    outputfile->Write();
    outputfile->Close();
}