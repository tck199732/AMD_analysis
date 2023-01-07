
#include "filter_e15190.hh"

struct IMQMD
{
    const static int MAX_MULTI = 128;
    int multi, Nc;
    double b_fm;
    std::array<int, MAX_MULTI> id;
    std::array<int, MAX_MULTI> A;
    std::array<int, MAX_MULTI> Z;
    std::array<double, MAX_MULTI> x_fm;
    std::array<double, MAX_MULTI> y_fm;
    std::array<double, MAX_MULTI> z_fm;
    std::array<double, MAX_MULTI> px_MeV;
    std::array<double, MAX_MULTI> py_MeV;
    std::array<double, MAX_MULTI> pz_MeV;
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
    std::array<double, MAX_MULTI> uball_id;
    std::array<double, MAX_MULTI> uball_x;
    std::array<double, MAX_MULTI> uball_y;
    std::array<double, MAX_MULTI> uball_z;

    // hira
    int hira_multi;
    std::array<int, MAX_MULTI> hira_N;
    std::array<int, MAX_MULTI> hira_Z;
    std::array<double, MAX_MULTI> hira_px;
    std::array<double, MAX_MULTI> hira_py;
    std::array<double, MAX_MULTI> hira_pz;
    std::array<double, MAX_MULTI> hira_id;
    std::array<double, MAX_MULTI> hira_x;
    std::array<double, MAX_MULTI> hira_y;
    std::array<double, MAX_MULTI> hira_z;

    // veto wall
    // neutron wall
};

IMQMD imqmd;
E15190 filtered_imqmd;

void Initialize_TChain(TChain *&chain, const std::vector<std::string> &pths)
{
    chain->SetBranchAddress("multi", &imqmd.multi);
    chain->SetBranchAddress("b_fm", &imqmd.b_fm);
    chain->SetBranchAddress("id", &imqmd.id[0]);
    chain->SetBranchAddress("A", &imqmd.A[0]);
    chain->SetBranchAddress("Z", &imqmd.Z[0]);
    chain->SetBranchAddress("x_fm", &imqmd.x_fm[0]);
    chain->SetBranchAddress("y_fm", &imqmd.y_fm[0]);
    chain->SetBranchAddress("z_fm", &imqmd.z_fm[0]);
    chain->SetBranchAddress("px_MeV", &imqmd.px_MeV[0]);
    chain->SetBranchAddress("py_MeV", &imqmd.py_MeV[0]);
    chain->SetBranchAddress("pz_MeV", &imqmd.pz_MeV[0]);

    chain->SetMakeClass(1);
    chain->SetBranchStatus("*", false);
    chain->SetBranchStatus("multi", true);
    chain->SetBranchStatus("b_fm", true);
    chain->SetBranchStatus("id", true);
    chain->SetBranchStatus("A", true);
    chain->SetBranchStatus("Z", true);
    chain->SetBranchStatus("x_fm", true);
    chain->SetBranchStatus("y_fm", true);
    chain->SetBranchStatus("z_fm", true);
    chain->SetBranchStatus("px_MeV", true);
    chain->SetBranchStatus("py_MeV", true);
    chain->SetBranchStatus("pz_MeV", true);

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
    tree->Branch("b", &filtered_imqmd.b, "b/D");

    // microball
    tree->Branch("uball_multi", &filtered_imqmd.uball_multi, "uball_multi/I");
    tree->Branch("uball_N", &filtered_imqmd.uball_N[0], "uball_N[uball_multi]/I");
    tree->Branch("uball_Z", &filtered_imqmd.uball_Z[0], "uball_Z[uball_multi]/I");
    tree->Branch("uball_px", &filtered_imqmd.uball_px[0], "uball_px[uball_multi]/D");
    tree->Branch("uball_py", &filtered_imqmd.uball_py[0], "uball_py[uball_multi]/D");
    tree->Branch("uball_pz", &filtered_imqmd.uball_pz[0], "uball_pz[uball_multi]/D");

    tree->Branch("uball_id", &filtered_imqmd.uball_id[0], "uball_id[uball_multi]/I");
    tree->Branch("uball_x", &filtered_imqmd.uball_x[0], "uball_x[uball_multi]/D");
    tree->Branch("uball_y", &filtered_imqmd.uball_y[0], "uball_y[uball_multi]/D");
    tree->Branch("uball_z", &filtered_imqmd.uball_z[0], "uball_z[uball_multi]/D");

    // hira
    tree->Branch("hira_multi", &filtered_imqmd.hira_multi, "hira_multi/I");
    tree->Branch("hira_N", &filtered_imqmd.hira_N[0], "hira_N[hira_multi]/I");
    tree->Branch("hira_Z", &filtered_imqmd.hira_Z[0], "hira_Z[hira_multi]/I");
    tree->Branch("hira_px", &filtered_imqmd.hira_px[0], "hira_px[hira_multi]/D");
    tree->Branch("hira_py", &filtered_imqmd.hira_py[0], "hira_py[hira_multi]/D");
    tree->Branch("hira_pz", &filtered_imqmd.hira_pz[0], "hira_pz[hira_multi]/D");

    tree->Branch("hira_id", &filtered_imqmd.hira_id[0], "hira_id[hira_multi]/I");
    tree->Branch("hira_x", &filtered_imqmd.hira_x[0], "hira_x[hira_multi]/D");
    tree->Branch("hira_y", &filtered_imqmd.hira_y[0], "hira_y[hira_multi]/D");
    tree->Branch("hira_z", &filtered_imqmd.hira_z[0], "hira_z[hira_multi]/D");
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

    TChain *chain = new TChain("cluster");
    Initialize_TChain(chain, input_pths);

    double betacms = Physics::GetReactionBeta(reaction);
    double rapidity_beam = Physics::GetBeamRapidity(reaction);

    Microball *uball_detector = GetMicroBall(reaction);
    HiRA *hira_detector = new HiRA();

    TTree *tree = new TTree("cluster", "");
    Initialize_TTree(tree);

    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        uball_detector->ResetCsI();
        hira_detector->ResetCounter();

        for (unsigned int i = 0; i < imqmd.multi; i++)
        {
            particle particle = {imqmd.A[i] - imqmd.Z[i], imqmd.Z[i], imqmd.px_MeV[i], imqmd.py_MeV[i], imqmd.pz_MeV[i]};
            particle.initialize(betacms);

            int uball_multi = uball_detector->GetCsIHits();
            int hira_multi = hira_detector->GetCountPass();

            if (ReadMicroballParticle(uball_detector, particle))
            {
                filtered_imqmd.uball_id[uball_multi] = imqmd.id[i];
                filtered_imqmd.uball_N[uball_multi] = imqmd.A[i] - imqmd.Z[i];
                filtered_imqmd.uball_Z[uball_multi] = imqmd.Z[i];
                filtered_imqmd.uball_x[uball_multi] = imqmd.x_fm[i];
                filtered_imqmd.uball_y[uball_multi] = imqmd.y_fm[i];
                filtered_imqmd.uball_z[uball_multi] = imqmd.z_fm[i];
                filtered_imqmd.uball_px[uball_multi] = imqmd.px_MeV[i];
                filtered_imqmd.uball_py[uball_multi] = imqmd.py_MeV[i];
                filtered_imqmd.uball_pz[uball_multi] = imqmd.pz_MeV[i];

                uball_detector->AddCsIHit(particle.theta_lab, particle.phi);
            }

            if (ReadHiRAParticle(hira_detector, particle))
            {
                filtered_imqmd.hira_id[hira_multi] = imqmd.id[i];
                filtered_imqmd.hira_N[hira_multi] = imqmd.A[i] - imqmd.Z[i];
                filtered_imqmd.hira_Z[hira_multi] = imqmd.Z[i];
                filtered_imqmd.hira_x[hira_multi] = imqmd.x_fm[i];
                filtered_imqmd.hira_y[hira_multi] = imqmd.y_fm[i];
                filtered_imqmd.hira_z[hira_multi] = imqmd.z_fm[i];
                filtered_imqmd.hira_px[hira_multi] = imqmd.px_MeV[i];
                filtered_imqmd.hira_py[hira_multi] = imqmd.py_MeV[i];
                filtered_imqmd.hira_pz[hira_multi] = imqmd.pz_MeV[i];
                hira_detector->CountPass();
            }
        }

        filtered_imqmd.b = imqmd.b_fm;
        filtered_imqmd.hira_multi = hira_detector->GetCountPass();
        filtered_imqmd.uball_multi = uball_detector->GetCsIHits();
        tree->Fill();
    }

    TFile *outputfile = new TFile(path_out.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}
