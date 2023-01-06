#include "anal.hh"
struct Event
{
    int multi;
    double bimp;
    std::vector<particle> particles;
};

struct EventCut
{
    std::array<int, 2> Nccut;
    std::array<double, 2> bcut;
    bool pass(const Event &event)
    {
        return (event.multi >= this->Nccut[0] && event.multi <= this->Nccut[1] && event.bimp >= this->bcut[0] && event.bimp < this->bcut[1]);
    }
};

struct Histograms
{
    std::string mode;
    std::vector<std::string> particlenames = {
        "n", "p", "d", "t", "3He", "4He", "coal_n", "coal_p"};

    double norm = 0.0;
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;
    void init();
    void fill(const particle &particle, const double &weight);
    void fill(const Event &event, const double &weight);
    void normalize();
    void write();
};

int main(int argc, char *argv[])
{
    // passe args
    std::string reaction = argv[1];
    std::string output_pth = argv[2];

    int nfiles_input = std::stoi(argv[3]);
    std::vector<std::string> input_pths21;
    std::vector<std::string> input_pths3;

    int narg = 4;
    for (int _ = 0; _ < nfiles_input; _++)
    {
        input_pths21.push_back(argv[narg]);
        narg++;
    }
    for (int _ = 0; _ < nfiles_input; _++)
    {
        input_pths3.push_back(argv[narg]);
        narg++;
    }

    // read args for event cut
    int Ncmin = std::stoi(argv[narg]);
    int Ncmax = std::stoi(argv[narg++]);
    double bmin = std::stod(argv[narg++]);
    double bmax = std::stod(argv[narg++]);

    // initialize Tchains and other structures
    TChain *chain21 = new TChain("AMD");
    TChain *chain3 = new TChain("AMD");
    Initialize_TChain(chain21, input_pths21, "filtered", "21");
    Initialize_TChain(chain3, input_pths3, "filtered", "3");

    Histograms histograms21 = {"21"};
    Histograms histograms3 = {"3"};
    Histograms histograms3_single_decay = {"3-single-decay"};
    histograms21.init();
    histograms3.init();
    histograms3_single_decay.init();

    double betacms = Physics::GetReactionBeta(reaction);
    double beam_rapidity = Physics::GetBeamRapidity(reaction);

    EventCut event_cut = {{Ncmin, Ncmax}, {bmin, bmax}};

    // counter of weights of each primary event
    int nevents21 = chain21->GetEntries();
    int nevents3 = chain3->GetEntries();
    std::vector<double> weights(nevents21, 0.);

    auto start3 = std::chrono::steady_clock::now();
    for (int ievt3 = 0; ievt3 < nevents3; ievt3++)
    {
        chain3->GetEntry(ievt3);

        Event event3 = {AMD.Nc, AMD.b};
        if (event_cut.pass(event3))
        {
            histograms3.norm += 1. / NDECAYS;
            if (ievt3 < nevents3 / NDECAYS)
            {
                histograms3_single_decay.norm += 1.;
            }
        }
        else
        {
            continue;
        }
        weights[ievt3 % nevents21] += 1.;

        for (unsigned int i = 0; i < event3.multi; i++)
        {
            particle particle = {AMD.N[i], AMD.Z[i], AMD.px[i], AMD.py[i], AMD.pz[i]};
            particle.autofill(betacms, beam_rapidity);

            if (ievt3 < nevents3 / NDECAYS)
            {
                histograms3_single_decay.fill(particle, 1.);
            }
            histograms3.fill(particle, 1. / NDECAYS);
        }
    }
    auto end3 = std::chrono::steady_clock::now();
    std::chrono::duration<double> timer3 = end3 - start3;

    auto start21 = std::chrono::steady_clock::now();
    for (int ievt21 = 0; ievt21 < nevents21; ievt21++)
    {
        chain21->GetEntry(ievt21);
        weights[ievt21] /= 1.0 * NDECAYS;

        Event event21 = {AMD.Nc, AMD.b};

        for (unsigned int i = 0; i < event21.multi; i++)
        {
            particle particle = {AMD.N[i], AMD.Z[i], AMD.px[i], AMD.py[i], AMD.pz[i]};
            particle.autofill(betacms, beam_rapidity);
            histograms21.fill(particle, weights[ievt21]);
        }

        histograms21.norm += weights[ievt21];
    }
    auto end21 = std::chrono::steady_clock::now();
    std::chrono::duration<double> timer21 = end21 - start21;

    // error correction for sampling 10 decay events
    for (auto &[pn, h2] : histograms3.h2_pta_rapidity_lab)
    {
        for (int i = 0; i < h2->GetNbinsX(); i++)
        {
            for (int j = 0; j < h2->GetNbinsY(); j++)
            {
                double error = histograms3_single_decay.h2_pta_rapidity_lab[pn]->GetBinError(i + 1, j + 1);
                h2->SetBinError(i + 1, j + 1, error);
            }
        }
    }

    // saving results
    TFile *outputfile = new TFile(output_pth.c_str(), "RECREATE");
    outputfile->cd();
    histograms21.normalize();
    histograms3.normalize();
    histograms3_single_decay.normalize();
    histograms21.write();
    histograms3.write();
    histograms3_single_decay.write();
    outputfile->Write();
    std::cout << "Elapsed Time reading table3 : " << timer3.count() << std::endl;
    std::cout << "Elapsed Time reading table21 : " << timer21.count() << std::endl;
}

void Histograms::init()
{
    for (const auto &pn : this->particlenames)
    {
        std::string hname = "h2_pta_rapidity_lab_" + this->mode + "_" + pn;
        this->h2_pta_rapidity_lab[pn] = new TH2D(hname.c_str(), "", 100, 0., 1., 600, 0, 600);
        this->h2_pta_rapidity_lab[pn]->Sumw2();
        this->h2_pta_rapidity_lab[pn]->SetDirectory(0);
    }
}

void Histograms::write()
{
    for (const auto &pn : this->particlenames)
    {
        this->h2_pta_rapidity_lab[pn]->Write();
    }
}

void Histograms::normalize()
{
    std::cout << "normalizing Histograms : " << this->norm << std::endl;
    for (const auto &pn : this->particlenames)
    {
        this->h2_pta_rapidity_lab[pn]->Scale(1. / this->norm);
    }
}

void Histograms::fill(const particle &particle, const double &weight)
{
    if (this->h2_pta_rapidity_lab.count(particle.name) == 1)
    {
        this->h2_pta_rapidity_lab[particle.name]->Fill(particle.rapidity_lab_normed, particle.pt, weight);
    }

    if (h2_pta_rapidity_lab.count("coal_p") == 1)
    {
        this->h2_pta_rapidity_lab["coal_p"]->Fill(particle.rapidity_lab_normed, particle.pt, weight * particle.zid);
    }
    if (h2_pta_rapidity_lab.count("coal_n") == 1)
    {
        this->h2_pta_rapidity_lab["coal_n"]->Fill(particle.rapidity_lab_normed, particle.pt, weight * particle.nid);
    }
    return;
}

void Histograms::fill(const Event &event, const double &weight)
{
    for (const auto &par : event.particles)
    {
        this->fill(par, weight);
    }
}
