
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

IMQMD imqmd;
IMQMD filtered_imqmd;

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
    tree->Branch("multi", &filtered_imqmd.multi, "multi/I");
    tree->Branch("Nc", &filtered_imqmd.Nc, "Nc/I");
    tree->Branch("b_fm", &filtered_imqmd.b_fm, "b_fm/D");
    tree->Branch("id", &filtered_imqmd.id[0], "z_fm[Nc]/I");
    tree->Branch("A", &filtered_imqmd.A[0], "A[Nc]/I");
    tree->Branch("Z", &filtered_imqmd.Z[0], "Z[Nc]/I");
    tree->Branch("x_fm", &filtered_imqmd.x_fm[0], "x_fm[Nc]/D");
    tree->Branch("y_fm", &filtered_imqmd.y_fm[0], "y_fm[Nc]/D");
    tree->Branch("z_fm", &filtered_imqmd.z_fm[0], "z_fm[Nc]/D");
    tree->Branch("px_MeV", &filtered_imqmd.px_MeV[0], "px_MeV[Nc]/D");
    tree->Branch("py_MeV", &filtered_imqmd.py_MeV[0], "py_MeV[Nc]/D");
    tree->Branch("pz_MeV", &filtered_imqmd.pz_MeV[0], "pz_MeV[Nc]/D");
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

    system_info *ReactionInfo = new system_info();
    double betacms = ReactionInfo->get_betacms(reaction);
    double rapidity_beam = ReactionInfo->get_rapidity_beam(reaction);

    MicroBall *uball_detector = new MicroBall();
    uball_detector->Configurate(reaction);

    TTree *tree = new TTree("cluster", "");
    Initialize_TTree(tree);

    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);
        filtered_imqmd.multi = imqmd.multi;
        filtered_imqmd.b_fm = imqmd.b_fm;
        uball_detector->Reset();

        int hira_counts = 0;
        for (unsigned int i = 0; i < imqmd.multi; i++)
        {
            particle particle = {imqmd.A[i] - imqmd.Z[i], imqmd.Z[i], imqmd.px_MeV[i] / imqmd.A[i], imqmd.py_MeV[i] / imqmd.A[i], imqmd.pz_MeV[i] / imqmd.A[i]};

            particle.autofill(betacms, rapidity_beam);
            uball_detector->ReadParticle(particle);

            if (hira_detector.pass_angle(particle) && hira_detector.pass_ekinlab(particle))
            {
                filtered_imqmd.id[hira_counts] = imqmd.id[i];
                filtered_imqmd.A[hira_counts] = imqmd.A[i];
                filtered_imqmd.Z[hira_counts] = imqmd.Z[i];
                filtered_imqmd.x_fm[hira_counts] = imqmd.x_fm[i];
                filtered_imqmd.y_fm[hira_counts] = imqmd.y_fm[i];
                filtered_imqmd.z_fm[hira_counts] = imqmd.z_fm[i];
                filtered_imqmd.px_MeV[hira_counts] = imqmd.px_MeV[i];
                filtered_imqmd.py_MeV[hira_counts] = imqmd.py_MeV[i];
                filtered_imqmd.pz_MeV[hira_counts] = imqmd.pz_MeV[i];
                hira_counts++;
            }
        }

        filtered_imqmd.Nc = uball_detector->GetNumberOfChargedParticles();
        tree->Fill();
    }

    TFile *outputfile = new TFile(path_out.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}
