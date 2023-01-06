#include "filter_e15190.hh"

struct AMD
{
    const static int MAX_MULTI = 128;
    // base
    int multi, Nc;
    double b;
    std::array<int, MAX_MULTI> N;
    std::array<int, MAX_MULTI> Z;
    std::array<double, MAX_MULTI> px;
    std::array<double, MAX_MULTI> py;
    std::array<double, MAX_MULTI> pz;

    // // table21
    // std::array<double, MAX_MULTI> ENG;
    // std::array<double, MAX_MULTI> LANG;
    // std::array<double, MAX_MULTI> JX;
    // std::array<double, MAX_MULTI> JY;
    // std::array<double, MAX_MULTI> JZ;

    // table3
    std::array<double, MAX_MULTI> J;
    std::array<double, MAX_MULTI> M;
    std::array<double, MAX_MULTI> WEIGHT;
    std::array<int, MAX_MULTI> iFRG;

    // // table21t
    // std::array<double, MAX_MULTI> t;
    // std::array<double, MAX_MULTI> x;
    // std::array<double, MAX_MULTI> y;
    // std::array<double, MAX_MULTI> z;
};

AMD amd;
AMD filtered_amd;

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &pths)
{
    chain->SetBranchAddress("multi", &amd.multi);
    chain->SetBranchAddress("b", &amd.b);
    chain->SetBranchAddress("N", &amd.N[0]);
    chain->SetBranchAddress("Z", &amd.Z[0]);
    chain->SetBranchAddress("px", &amd.px[0]);
    chain->SetBranchAddress("py", &amd.py[0]);
    chain->SetBranchAddress("pz", &amd.pz[0]);

    // table3
    chain->SetBranchAddress("J", &amd.J[0]);
    chain->SetBranchAddress("M", &amd.M[0]);
    chain->SetBranchAddress("WEIGHT", &amd.WEIGHT[0]);
    chain->SetBranchAddress("iFRG", &amd.iFRG[0]);

    chain->SetMakeClass(1);
    chain->SetBranchStatus("*", false);
    chain->SetBranchStatus("multi", true);
    chain->SetBranchStatus("b", true);
    chain->SetBranchStatus("N", true);
    chain->SetBranchStatus("Z", true);
    chain->SetBranchStatus("px", true);
    chain->SetBranchStatus("py", true);
    chain->SetBranchStatus("pz", true);

    // table3
    chain->SetBranchStatus("J", true);
    chain->SetBranchStatus("M", true);
    chain->SetBranchStatus("WEIGHT", true);
    chain->SetBranchStatus("iFRG", true);

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
    tree->Branch("multi", &filtered_amd.multi, "multi/I");
    tree->Branch("Nc", &filtered_amd.Nc, "Nc/I");
    tree->Branch("b", &filtered_amd.b, "b_fm/D");
    tree->Branch("N", &filtered_amd.N[0], "N[Nc]/I");
    tree->Branch("Z", &filtered_amd.Z[0], "Z[Nc]/I");
    tree->Branch("px", &filtered_amd.px[0], "px[Nc]/D");
    tree->Branch("py", &filtered_amd.py[0], "py[Nc]/D");
    tree->Branch("pz", &filtered_amd.pz[0], "pz[Nc]/D");
    // table3
    tree->Branch("J", &filtered_amd.J[0], "J[Nc]/D");
    tree->Branch("M", &filtered_amd.M[0], "M[Nc]/D");
    tree->Branch("WEIGHT", &filtered_amd.WEIGHT[0], "WEIGHT[Nc]/D");
    tree->Branch("iFRG", &filtered_amd.iFRG[0], "iFRG[Nc]/D");
}

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];

    std::string path_out = argv[2];
    int narg = 3;
    int number_files = std::stoi(argv[3]);
    std::vector<std::string> input_pths;
    for (int _ = 0; _ < number_files; _++)
    {
        input_pths.push_back(argv[++narg]);
    }

    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, input_pths);

    double betacms = Physics::GetReactionBeta(reaction);
    double rapidity_beam = Physics::GetBeamRapidity(reaction);

    Microball *uball_detector = GetMicroBall(reaction);
    HiRA *hira_detector = new HiRA();

    TTree *tree = new TTree("AMD", "");
    Initialize_TTree(tree);

    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        uball_detector->ResetCsI();

        int hira_counts = 0;
        for (unsigned int i = 0; i < amd.multi; i++)
        {
            particle particle = {amd.N[i], amd.Z[i], amd.px[i], amd.py[i], amd.pz[i]};
            particle.autofill(betacms, rapidity_beam);

            ReadMicroballParticle(uball_detector, particle);
            if (ReadHiRAParticle(hira_detector, particle))
            {
                filtered_amd.N[hira_counts] = amd.N[i];
                filtered_amd.Z[hira_counts] = amd.Z[i];
                filtered_amd.px[hira_counts] = amd.px[i];
                filtered_amd.py[hira_counts] = amd.py[i];
                filtered_amd.pz[hira_counts] = amd.pz[i];

                // table3
                filtered_amd.J[hira_counts] = amd.J[i];
                filtered_amd.M[hira_counts] = amd.M[i];
                filtered_amd.WEIGHT[hira_counts] = amd.WEIGHT[i];
                filtered_amd.iFRG[hira_counts] = amd.iFRG[i];
                hira_counts++;
            }
        }

        filtered_amd.multi = hira_counts;
        filtered_amd.Nc = uball_detector->GetCsIHits();
        filtered_amd.b = amd.b;
        tree->Fill();
    }

    TFile *outputfile = new TFile(path_out.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}
