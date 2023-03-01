#include "filter_e15190.hh"

struct AMD
{
    const static int MAX_MULTI = 128;
    // base
    int multi;
    double b;
    std::array<int, MAX_MULTI> N;
    std::array<int, MAX_MULTI> Z;
    std::array<double, MAX_MULTI> px;
    std::array<double, MAX_MULTI> py;
    std::array<double, MAX_MULTI> pz;

    // table3
    std::array<double, MAX_MULTI> J;
    std::array<double, MAX_MULTI> M;
    std::array<double, MAX_MULTI> WEIGHT;
    std::array<int, MAX_MULTI> iFRG;
};
struct E15190
{
    const static int MAX_MULTI = 128;
    // impact parameter
    double b;

    // microball
    int uball_multi;
    std::array<int, MAX_MULTI> uball_N;
    std::array<int, MAX_MULTI> uball_Z;
    std::array<double, MAX_MULTI> uball_px;
    std::array<double, MAX_MULTI> uball_py;
    std::array<double, MAX_MULTI> uball_pz;

    // hira
    int hira_multi;
    std::array<int, MAX_MULTI> hira_N;
    std::array<int, MAX_MULTI> hira_Z;
    std::array<double, MAX_MULTI> hira_px;
    std::array<double, MAX_MULTI> hira_py;
    std::array<double, MAX_MULTI> hira_pz;

    // veto wall
    // neutron wall
};

AMD amd;
E15190 filtered_amd;

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &pths)
{
    chain->SetBranchAddress("multi", &amd.multi);
    chain->SetBranchAddress("b", &amd.b);
    chain->SetBranchAddress("N", &amd.N[0]);
    chain->SetBranchAddress("Z", &amd.Z[0]);
    chain->SetBranchAddress("px", &amd.px[0]);
    chain->SetBranchAddress("py", &amd.py[0]);
    chain->SetBranchAddress("pz", &amd.pz[0]);

    chain->SetMakeClass(1);
    chain->SetBranchStatus("*", false);
    chain->SetBranchStatus("multi", true);
    chain->SetBranchStatus("b", true);
    chain->SetBranchStatus("N", true);
    chain->SetBranchStatus("Z", true);
    chain->SetBranchStatus("px", true);
    chain->SetBranchStatus("py", true);
    chain->SetBranchStatus("pz", true);

    for (auto &pth : pths)
    {
        if (!fs::exists(pth))
        {
            std::string msg = Form("%s does not exists.", pth.c_str());
            throw std::invalid_argument(msg.c_str());
        }
        chain->Add(pth.c_str());
    }
}

void Initialize_TTree(TTree *&tree)
{
    // impact parameter
    tree->Branch("b", &filtered_amd.b, "b/D");

    // microball
    tree->Branch("uball_multi", &filtered_amd.uball_multi, "uball_multi/I");
    tree->Branch("uball_N", &filtered_amd.uball_N[0], "uball_N[uball_multi]/I");
    tree->Branch("uball_Z", &filtered_amd.uball_Z[0], "uball_Z[uball_multi]/I");
    tree->Branch("uball_px", &filtered_amd.uball_px[0], "uball_px[uball_multi]/D");
    tree->Branch("uball_py", &filtered_amd.uball_py[0], "uball_py[uball_multi]/D");
    tree->Branch("uball_pz", &filtered_amd.uball_pz[0], "uball_pz[uball_multi]/D");

    // hira
    tree->Branch("hira_multi", &filtered_amd.hira_multi, "hira_multi/I");
    tree->Branch("hira_N", &filtered_amd.hira_N[0], "hira_N[hira_multi]/I");
    tree->Branch("hira_Z", &filtered_amd.hira_Z[0], "hira_Z[hira_multi]/I");
    tree->Branch("hira_px", &filtered_amd.hira_px[0], "hira_px[hira_multi]/D");
    tree->Branch("hira_py", &filtered_amd.hira_py[0], "hira_py[hira_multi]/D");
    tree->Branch("hira_pz", &filtered_amd.hira_pz[0], "hira_pz[hira_multi]/D");
}

int main(int argc, char *argv[])
{
    ArgumentParser argparser(argc, argv);

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, argparser.input_files);

    double betacms = Physics::GetReactionBeta(argparser.reaction);
    double rapidity_beam = Physics::GetBeamRapidity(argparser.reaction); // not needed

    Microball *microball = new Microball();
    Initialize_MicroBall(microball, argparser.reaction);
    HiRA *hira = new HiRA();

    TTree *tree = new TTree("AMD", "");
    Initialize_TTree(tree);

    ProgressBar bar(chain->GetEntries(), argparser.reaction);
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        microball->ResetCsIHitMap();
        hira->ResetCounter();
        for (unsigned int i = 0; i < amd.multi; i++)
        {
            particle particle = {amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i]};
            particle.initialize(betacms);
            // phi is calculated according to microball detector, if the particle is not covered by microball, phi is not correct and should be in the range of [-pi, pi].
            correct_phi_value(particle, microball);

            int uball_multi = microball->GetCsIHits();
            int hira_multi = hira->GetCountPass();

            if (ReadMicroballParticle(microball, particle))
            {
                filtered_amd.uball_N[uball_multi] = particle.N;
                filtered_amd.uball_Z[uball_multi] = particle.Z;
                filtered_amd.uball_px[uball_multi] = particle.px;
                filtered_amd.uball_py[uball_multi] = particle.py;
                filtered_amd.uball_pz[uball_multi] = particle.pz_lab;

                microball->AddCsIHit(particle.theta_lab, particle.phi);
            }

            if (ReadHiRAParticle(hira, particle))
            {
                filtered_amd.hira_N[hira_multi] = particle.N;
                filtered_amd.hira_Z[hira_multi] = particle.Z;
                filtered_amd.hira_px[hira_multi] = particle.px;
                filtered_amd.hira_py[hira_multi] = particle.py;
                filtered_amd.hira_pz[hira_multi] = particle.pz_lab;
                hira->CountPass();
            }
        }

        filtered_amd.hira_multi = hira->GetCountPass();
        filtered_amd.uball_multi = microball->GetCsIHits();
        filtered_amd.b = amd.b;
        // if microball multi is 0, in experiment we don't see the event. We still keep the event here as this data can be easily removed in the analysis.
        tree->Fill();
        bar.Update();
    }

    TFile *outputfile = new TFile(argparser.output_file.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}
