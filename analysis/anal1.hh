
#include "../src/system_info.cpp"
#include "../src/detector_uball.cpp"
#include "../src/particle.cpp"
#include "../src/root_io.cpp"

#include <vector>
#include <string>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

fs::path PROJECT_DIR = "/home/kin/Desktop/amd_analysis";
const int ndecays = 10;

struct event
{
    int multi;
    double bimp;
    std::vector<particle> particles;
};

struct histograms
{
    std::string reaction, mode;
    TH1D *h1_invariant_mass_5Li;
    TH1D *h1_invariant_mass_6Li;
    TH1D *h1_invariant_mass_7Li;
    TH1D *h1_invariant_mass_7Be;
    TH1D *h1_invariant_mass_8Be;

    std::map<std::pair<int, int>, double> mass_table; // nz to mass

    double mass_5Li = 4669.149399085722;
    double mass_6Li = 5603.051494942865;
    double mass_7Li = 6535.36582155001;
    double mass_7Be = 6536.22771492001;
    double mass_8Be = 7456.894491337155;

    double weight = 1. / ndecays;
    double norm = 0.0;
    void init();
    void fill(const event &event);
    void fill_pairs(TH1D *h1, const std::vector<particle> collection1, const std::vector<particle> collection2, const double &m0);
    void fill_pairs(TH1D *h1, const std::vector<particle> collection, const double &m0);
    void normalize();
    void write();
};

void histograms::init()
{

    mass_table[{0, 1}] = 938.7830734811444;
    mass_table[{1, 1}] = 1876.1239277292887;
    mass_table[{2, 1}] = 2809.432118151433;
    mass_table[{1, 2}] = 2809.4135261314327;
    mass_table[{2, 2}] = 3728.4013255385776;

    this->h1_invariant_mass_8Be = new TH1D("h1_invariant_mass_8Be", "", 600, -10, 20);
    this->h1_invariant_mass_7Be = new TH1D("h1_invariant_mass_7Be", "", 600, -10, 20);
    this->h1_invariant_mass_7Li = new TH1D("h1_invariant_mass_7Li", "", 600, -10, 20);
    this->h1_invariant_mass_6Li = new TH1D("h1_invariant_mass_6Li", "", 600, -10, 20);
    this->h1_invariant_mass_5Li = new TH1D("h1_invariant_mass_5Li", "", 600, -10, 20);

    this->h1_invariant_mass_8Be->Sumw2();
    this->h1_invariant_mass_7Be->Sumw2();
    this->h1_invariant_mass_7Li->Sumw2();
    this->h1_invariant_mass_6Li->Sumw2();
    this->h1_invariant_mass_5Li->Sumw2();

    this->h1_invariant_mass_8Be->SetDirectory(0);
    this->h1_invariant_mass_7Be->SetDirectory(0);
    this->h1_invariant_mass_7Li->SetDirectory(0);
    this->h1_invariant_mass_6Li->SetDirectory(0);
    this->h1_invariant_mass_5Li->SetDirectory(0);
}

void histograms::fill(const event &event)
{
    std::vector<particle> proton_particles;
    std::vector<particle> deuteron_particles;
    std::vector<particle> triton_particles;
    std::vector<particle> He3_particles;
    std::vector<particle> alpha_particles;

    for (auto &par : event.particles)
    {
        if (par.zid == 1 && par.nid == 0)
        {
            proton_particles.push_back(par);
        }
        else if (par.zid == 1 && par.nid == 1)
        {
            deuteron_particles.push_back(par);
        }
        else if (par.zid == 1 && par.nid == 2)
        {
            triton_particles.push_back(par);
        }
        else if (par.zid == 2 && par.nid == 1)
        {
            He3_particles.push_back(par);
        }
        else if (par.zid == 2 && par.nid == 2)
        {
            alpha_particles.push_back(par);
        }
    }

    this->fill_pairs(this->h1_invariant_mass_5Li, proton_particles, alpha_particles, this->mass_5Li);

    this->fill_pairs(this->h1_invariant_mass_6Li, deuteron_particles, alpha_particles, this->mass_6Li);

    this->fill_pairs(this->h1_invariant_mass_7Li, alpha_particles, triton_particles, this->mass_7Li);

    this->fill_pairs(this->h1_invariant_mass_7Be, He3_particles, alpha_particles, this->mass_7Be);

    this->fill_pairs(this->h1_invariant_mass_8Be, alpha_particles, this->mass_8Be);

    this->norm += this->weight;
}

void histograms::fill_pairs(TH1D *h1, const std::vector<particle> collection, const double &m0)
{
    int size = collection.size();
    if (size <= 1)
    {
        return;
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = i + 1; j < size; j++)
        {

            int zid1 = collection[i].zid;
            int zid2 = collection[j].zid;
            int nid1 = collection[i].nid;
            int nid2 = collection[j].nid;
            int aid1 = collection[i].aid;
            int aid2 = collection[j].aid;

            double px1 = collection[i].px * aid1;
            double px2 = collection[j].px * aid2;
            double py1 = collection[i].py * aid1;
            double py2 = collection[j].py * aid2;
            double pz1 = collection[i].pz * aid1;
            double pz2 = collection[j].pz * aid2;

            double m1 = this->mass_table[{nid1, zid1}];
            double m2 = this->mass_table[{nid2, zid2}];

            double e1 = collection[i].ekincms * aid1 + m1;
            double e2 = collection[j].ekincms * aid2 + m2;

            double p1p2 = px1 * px2 + py1 * py2 + pz1 * pz2;
            double M2 = pow(m1, 2.) + pow(m2, 2.) + 2 * (e1 * e2 - p1p2);
            h1->Fill(TMath::Sqrt(M2) - m0);
        }
    }
}

void histograms::fill_pairs(TH1D *h1, const std::vector<particle> collection1, const std::vector<particle> collection2, const double &m0)
{

    int size1 = collection1.size();
    int size2 = collection2.size();
    if (size1 == 0 || size2 == 0)
    {
        return;
    }

    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size2; j++)
        {

            int zid1 = collection1[i].zid;
            int zid2 = collection2[j].zid;
            int nid1 = collection1[i].nid;
            int nid2 = collection2[j].nid;
            int aid1 = collection1[i].aid;
            int aid2 = collection2[j].aid;

            double px1 = collection1[i].px * aid1;
            double px2 = collection2[j].px * aid2;
            double py1 = collection1[i].py * aid1;
            double py2 = collection2[j].py * aid2;
            double pz1 = collection1[i].pz * aid1;
            double pz2 = collection2[j].pz * aid2;

            double m1 = this->mass_table[{nid1, zid1}];
            double m2 = this->mass_table[{nid2, zid2}];

            double e1 = collection1[i].ekincms * aid1 + m1;
            double e2 = collection2[j].ekincms * aid2 + m2;

            double p1p2 = px1 * px2 + py1 * py2 + pz1 * pz2;
            double M2 = pow(m1, 2.) + pow(m2, 2.) + 2 * (e1 * e2 - p1p2);
            h1->Fill(TMath::Sqrt(M2) - m0);
        }
    }
}

void histograms::normalize()
{
    std::cout << "mode : " << this->mode << "\tnormalization : " << this->norm << std::endl;
    this->h1_invariant_mass_8Be->Scale(1. / this->norm);
    this->h1_invariant_mass_7Be->Scale(1. / this->norm);
    this->h1_invariant_mass_7Li->Scale(1. / this->norm);
    this->h1_invariant_mass_6Li->Scale(1. / this->norm);
    this->h1_invariant_mass_5Li->Scale(1. / this->norm);
}
void histograms::write()
{
    this->h1_invariant_mass_8Be->Write();
    this->h1_invariant_mass_7Be->Write();
    this->h1_invariant_mass_7Li->Write();
    this->h1_invariant_mass_6Li->Write();
    this->h1_invariant_mass_5Li->Write();
}
