#include "microball_test.hh"

struct particle
{
    int N, Z;
    double px, py, pz_cms;

    double kinergy_lab;
    double theta_lab, phi;
    double rapidity_lab;
    void initialize(const double &betacms);
};

void correct_phi_value(particle &part, Microball *&microball)
{
    int ring = microball->GetRingID(part.theta_lab);
    if (ring == -1)
    {
        return;
    }
    double phi_min_in_ring = microball->GetPhiMinInRing(ring);
    double phi_max_in_ring = microball->GetPhiMaxInRing(ring);

    if (part.phi < phi_min_in_ring)
    {
        part.phi += 360;
    }
    if (part.phi > phi_max_in_ring)
    {
        part.phi -= 360;
    }
}

bool PassMicroballParticle(Microball *mb, particle &part)
{
    bool pass_charge = mb->IsChargedParticle(part.Z);
    bool pass_coverage = mb->IsCovered(part.theta_lab, part.phi);
    bool pass_threshold = mb->IsAccepted(part.kinergy_lab, part.theta_lab, part.N + part.Z, part.Z);

    if (mb->Get_Is_apply_cut_charged_particle() == false && !pass_charge)
    {
        std::cout << "incorrect charge cut" << std::endl;
    }
    if (mb->Get_Is_apply_cut_coverage() == false && !pass_coverage)
    {
        std::cout << "incorrect coverage cut" << std::endl;
    }
    if (mb->Get_Is_apply_cut_kinergy() == false && !pass_threshold)
    {
        std::cout << "incorrect threshold cut" << std::endl;
    }
    return pass_charge && pass_coverage && pass_threshold;
}

int main(int argc, char *argv[])
{
    ArgumentParser argparser = ArgumentParser(argc, argv);

    std::string reaction = argparser.reaction;
    double betacms = Physics::GetReactionBeta(reaction);
    double beam_rapidity = Physics::GetBeamRapidity(reaction);
    std::cout << "reaction: " << reaction << std::endl;
    std::cout << "beta: " << betacms << std::endl;
    std::cout << "beam rapidity: " << beam_rapidity << std::endl;

    ReadAMETable("../../database/ame/ame_mass.txt", ame_mass_table);

    TChain *chain = new TChain("AMD", "");
    chain->Add(argparser.input_file.c_str());
    Initialize_TChain(chain);

    Microball *microball = new Microball();
    microball->Set_Is_apply_cut_charged_particle(argparser.is_apply_cut_charged_particle);
    microball->Set_Is_apply_cut_coverage(argparser.is_apply_cut_coverage);
    microball->Set_Is_apply_cut_kinergy(argparser.is_apply_cut_kinergy);
    microball->Set_Is_apply_cut_multiple_hit(argparser.is_apply_cut_multiple_hit);

    microball->ConfigurateSetup(reaction, "../../database/e15190/microball/acceptance/config.dat");
    microball->ReadGeometryMap("../../database/e15190/microball/acceptance/geometry.dat");
    microball->ReadThresholdKinergyMap("../../database/e15190/microball/acceptance/fitted_threshold.dat");

    // microball->ViewDetectorSetupMap();
    // microball->ViewGeometryMap();
    // microball->ViewThresholdKinergyMap();

    TH2D *h2_multi_b = new TH2D("h2_multi_b", "", 80, -0.5, 79.5, 100, 0., 10.);
    TH2D *h2_theta_phi = new TH2D("h2_theta_phi", "", 390, -30., 360., 180, 0., 180.);
    h2_multi_b->Sumw2();
    h2_theta_phi->Sumw2();

    for (int i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);
        microball->ResetCsIHitMap();
        for (int j = 0; j < amd.multi; j++)
        {
            particle part = {amd.N[j], amd.Z[j], amd.px[j], amd.py[j], amd.pz[j]};
            part.initialize(betacms);
            if (microball->Get_Is_apply_cut_coverage())
            {
                correct_phi_value(part, microball);
            }
            else
            {
                if (part.phi < 0)
                {
                    part.phi += 360;
                }
            }

            if (PassMicroballParticle(microball, part))
            {
                h2_theta_phi->Fill(part.phi, part.theta_lab);
                microball->AddCsIHit(part.theta_lab, part.phi);
            }
        }
        int Mch = microball->GetCsIHits();
        h2_multi_b->Fill(Mch, amd.b);
    }

    TFile *f = new TFile(argparser.output_file.c_str(), "RECREATE");
    h2_multi_b->Write();
    h2_theta_phi->Write();
    f->Close();

    return 0;
}

void particle::initialize(const double &betacms)
{
    double mass = ame_mass_table[{this->Z, this->N + this->Z}];
    double A = this->N + this->Z;
    this->px = this->px * A;
    this->py = this->py * A;
    this->pz_cms = this->pz_cms * A;

    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));
    double pmag_trans = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    double pmag_cms = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz_cms, 2.));

    double kinergy_cms = TMath::Sqrt(pow(pmag_cms, 2.) + pow(mass, 2.)) - mass;

    double pz_lab = gamma * (this->pz_cms + betacms * (kinergy_cms + mass));
    double pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(pz_lab, 2.));

    this->theta_lab = TMath::ATan2(pmag_trans, pz_lab) * TMath::RadToDeg();

    this->kinergy_lab = TMath::Sqrt(pow(pmag_lab, 2.) + pow(mass, 2.)) - mass;
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg(); // from -180 to 180 here, need to convert to uball coordinate for each ring

    this->rapidity_lab = 0.5 * TMath::Log((this->kinergy_lab + pz_lab + mass) / (this->kinergy_lab - pz_lab + mass));
}
