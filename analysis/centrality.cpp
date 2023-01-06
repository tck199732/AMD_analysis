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
    bool pass(const Event &event);
};

bool EventCut::pass(const Event &event)
{
    return (event.multi >= this->Nccut[0] && event.multi <= this->Nccut[1] && event.bimp >= this->bcut[0] && event.bimp < this->bcut[1]);
}

struct Histograms
{
    TH2D *h2_multi_b;
    void init();
    void fill(const Event &event, const double &weight = 1.);
    void normalize(const double &norm);
    void write();
};

void Histograms::init()
{
    this->h2_multi_b = new TH2D("h2_multi_b", "", 80, -0.5, 80.5, 1000, 0., 10.);
    this->h2_multi_b->Sumw2();
    this->h2_multi_b->SetDirectory(0);
}
void Histograms::write()
{
    this->h2_multi_b->Write();
}

void Histograms::fill(const Event &event, const double &weight)
{
    this->h2_multi_b->Fill(event.multi, event.bimp, weight);
}

void Histograms::normalize(const double &norm)
{
    this->h2_multi_b->Scale(1. / norm);
}

int main(int argc, char *argv[])
{
    // parse args
    std::string analysis = argv[1]; // default or filtered
    std::string reaction = argv[2];
    std::string mode = argv[3];
    std::string output_pth = argv[4];

    std::vector<std::string> input_pths;

    for (int i = 5; i < argc; i++)
    {
        input_pths.push_back(argv[i]);
    }

    // set up TChain and the data structure
    TChain *chain = new TChain("AMD");
    Initialize_TChain(chain, input_pths, analysis, mode);

    // set up histogram and event cut
    Histograms histograms;
    histograms.init();
    EventCut event_cut = {{1, 100}, {0., 10.}};

    // counter and weight of histogram
    int event_counter = 0;
    double weight = mode == "3" ? 1. / NDECAYS : 1.;

    // main analysis
    for (int ievt = 0; ievt < chain->GetEntries(); ievt++)
    {
        chain->GetEntry(ievt);

        int multi = analysis == "default" ? AMD.multi : AMD.Nc;
        Event event = {multi, AMD.b};
        if (event_cut.pass(event))
        {
            histograms.fill(event, weight);
            event_counter++;
        }
    }
    histograms.normalize(event_counter / NDECAYS);

    // saving result
    TFile *output_file = new TFile(output_pth.c_str(), "RECREATE");
    output_file->cd();
    histograms.write();
    output_file->Write();
    output_file->Close();
}