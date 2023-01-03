void centrality_temp()
{

    TChain *c = new TChain("AMD");
    // c->Add("/data/kin/amd/feb2022/b10fm/Ca40Ni58_En140MeV_SkM_30k_table3.root");
    c->Add("/data/kin/amd/feb2022/b10fm/Ca48Ni64_En140MeV_SkM_37k_table3.root");
    int multi;
    double b;

    c->SetBranchAddress("multi", &multi);
    c->SetBranchAddress("b", &b);

    TH2D *h2_multi_b = new TH2D("h2_multi_b", "", 80, -0.5, 79.5, 100, 0, 10);

    int nevents = c->GetEntries();
    for (int n = 0; n < nevents; n++)
    {
        c->GetEntry(n);
        h2_multi_b->Fill(multi, b);
    }

    TFile *f = new TFile("h2_multi_b_Ca48Ni64_En140MeV_SkM.root", "RECREATE");
    h2_multi_b->Write();
    f->Write();
    f->Close();
}