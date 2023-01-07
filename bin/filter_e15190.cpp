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
    std::array<double, MAX_MULTI> uball_J;
    std::array<double, MAX_MULTI> uball_M;
    std::array<double, MAX_MULTI> uball_WEIGHT;
    std::array<int, MAX_MULTI> uball_iFRG;

    // hira
    int hira_multi;
    std::array<int, MAX_MULTI> hira_N;
    std::array<int, MAX_MULTI> hira_Z;
    std::array<double, MAX_MULTI> hira_px;
    std::array<double, MAX_MULTI> hira_py;
    std::array<double, MAX_MULTI> hira_pz;
    std::array<double, MAX_MULTI> hira_J;
    std::array<double, MAX_MULTI> hira_M;
    std::array<double, MAX_MULTI> hira_WEIGHT;
    std::array<int, MAX_MULTI> hira_iFRG;

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
    // impact parameter
    tree->Branch("b", &filtered_amd.b, "b/D");

    // microball
    tree->Branch("uball_multi", &filtered_amd.uball_multi, "uball_multi/I");
    tree->Branch("uball_N", &filtered_amd.uball_N[0], "uball_N[uball_multi]/I");
    tree->Branch("uball_Z", &filtered_amd.uball_Z[0], "uball_Z[uball_multi]/I");
    tree->Branch("uball_px", &filtered_amd.uball_px[0], "uball_px[uball_multi]/D");
    tree->Branch("uball_py", &filtered_amd.uball_py[0], "uball_py[uball_multi]/D");
    tree->Branch("uball_pz", &filtered_amd.uball_pz[0], "uball_pz[uball_multi]/D");
    tree->Branch("uball_J", &filtered_amd.uball_J[0], "uball_J[uball_multi]/D");
    tree->Branch("uball_M", &filtered_amd.uball_M[0], "uball_M[uball_multi]/D");
    tree->Branch("uball_WEIGHT", &filtered_amd.uball_WEIGHT[0], "uball_WEIGHT[uball_multi]/D");
    tree->Branch("uball_iFRG", &filtered_amd.uball_iFRG[0], "uball_iFRG[uball_multi]/D");

    // hira
    tree->Branch("hira_multi", &filtered_amd.hira_multi, "hira_multi/I");
    tree->Branch("hira_N", &filtered_amd.hira_N[0], "hira_N[hira_multi]/I");
    tree->Branch("hira_Z", &filtered_amd.hira_Z[0], "hira_Z[hira_multi]/I");
    tree->Branch("hira_px", &filtered_amd.hira_px[0], "hira_px[hira_multi]/D");
    tree->Branch("hira_py", &filtered_amd.hira_py[0], "hira_py[hira_multi]/D");
    tree->Branch("hira_pz", &filtered_amd.hira_pz[0], "hira_pz[hira_multi]/D");
    tree->Branch("hira_J", &filtered_amd.hira_J[0], "hira_J[hira_multi]/D");
    tree->Branch("hira_M", &filtered_amd.hira_M[0], "hira_M[hira_multi]/D");
    tree->Branch("hira_WEIGHT", &filtered_amd.hira_WEIGHT[0], "hira_WEIGHT[hira_multi]/D");
    tree->Branch("hira_iFRG", &filtered_amd.hira_iFRG[0], "hira_iFRG[hira_multi]/D");
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
    double rapidity_beam = Physics::GetBeamRapidity(reaction); // not needed

    Microball *uball_detector = GetMicroBall(reaction);
    HiRA *hira_detector = new HiRA();

    TTree *tree = new TTree("AMD", "");
    Initialize_TTree(tree);

    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        uball_detector->ResetCsI();
        hira_detector->ResetCounter();
        for (unsigned int i = 0; i < amd.multi; i++)
        {
            int A = amd.N[i] + amd.Z[i];
            double px = amd.px[i] * A;
            double py = amd.py[i] * A;
            double pz_cms = amd.pz[i] * A;

            particle particle = {amd.N[i], amd.Z[i], px, py, pz_cms};
            particle.initialize(betacms);

            int uball_multi = uball_detector->GetCsIHits();
            int hira_multi = hira_detector->GetCountPass();

            if (ReadMicroballParticle(uball_detector, particle))
            {
                filtered_amd.uball_N[uball_multi] = amd.N[i];
                filtered_amd.uball_Z[uball_multi] = amd.Z[i];
                filtered_amd.uball_px[uball_multi] = px;
                filtered_amd.uball_py[uball_multi] = py;
                filtered_amd.uball_pz[uball_multi] = pz_cms;

                // table3
                filtered_amd.uball_J[uball_multi] = amd.J[i];
                filtered_amd.uball_M[uball_multi] = amd.M[i];
                filtered_amd.uball_WEIGHT[uball_multi] = amd.WEIGHT[i];
                filtered_amd.uball_iFRG[uball_multi] = amd.iFRG[i];

                uball_detector->AddCsIHit(particle.theta_lab, particle.phi);
            }

            if (ReadHiRAParticle(hira_detector, particle))
            {
                filtered_amd.hira_N[hira_multi] = amd.N[i];
                filtered_amd.hira_Z[hira_multi] = amd.Z[i];
                filtered_amd.hira_px[hira_multi] = px;
                filtered_amd.hira_py[hira_multi] = py;
                filtered_amd.hira_pz[hira_multi] = pz_cms;

                // table3
                filtered_amd.hira_J[hira_multi] = amd.J[i];
                filtered_amd.hira_M[hira_multi] = amd.M[i];
                filtered_amd.hira_WEIGHT[hira_multi] = amd.WEIGHT[i];
                filtered_amd.hira_iFRG[hira_multi] = amd.iFRG[i];
                hira_detector->CountPass();
            }
        }

        filtered_amd.hira_multi = hira_detector->GetCountPass();
        filtered_amd.uball_multi = uball_detector->GetCsIHits();
        filtered_amd.b = amd.b;
        tree->Fill();
    }

    TFile *outputfile = new TFile(path_out.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}
