
#include "multiplicity.hh"

struct manager
{
    std::string reaction, mode;
    fs::path path_data;
    fs::path path_out;
    int Ncmin, Ncmax;
    double bmin, bmax;
    system_info *sys_info;
    double betacms, rapidity_beam;

    RootReader *reader;
    histograms hist;
    eventcut event_cut;
    TFile *outputfile;

    int nevents;
    void init();
    void read();
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string mode = argv[2];
    std::string path_data = argv[3];
    std::string path_out = argv[4];
    int Ncmin = 1;
    int Ncmax = 25;
    double bmin = 0.0;
    double bmax = 3.0;
    if (argc >= 6)
    {
        Ncmin = std::stoi(argv[5]);
    }
    if (argc >= 7)
    {
        Ncmax = std::stoi(argv[6]);
    }
    if (argc >= 8)
    {
        bmin = std::stod(argv[7]);
    }
    if (argc >= 9)
    {
        bmax = std::stod(argv[8]);
    }

    manager manager = {reaction, mode, path_data, path_out, Ncmin, Ncmax, bmin, bmax};
    manager.init();
    manager.read();
    manager.finish();
    return 0;
}

void manager::init()
{
    if (!fs::exists(this->path_data))
    {
        throw std::invalid_argument("path_data3 does not exist.");
    }

    std::vector<branch> branches = {
        {"Nc", "int"},
        {"multi", "int"},
        {"b", "double"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
    };

    for (auto &br : branches)
    {
        br.autofill();
    }

    this->reader = new RootReader("AMD");
    reader->add_file(fs::absolute(this->path_data));
    reader->set_branches(branches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->event_cut = {{this->Ncmin, this->Ncmax}, {this->bmin, this->bmax}};

    this->hist = {this->reaction, this->mode, this->betacms, this->rapidity_beam};
    this->hist.init();
    this->nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;
}

void manager::read()
{
    int Nc, multi;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    for (int ievt = 0; ievt < this->nevents; ievt++)
    {
        std::map<std::string, std::any> map = this->reader->get_entry(ievt);
        try
        {
            Nc = std::any_cast<int>(map["Nc"]);
            bimp = std::any_cast<double>(map["b"]);
            multi = std::any_cast<int>(map["multi"]);
            fn = std::any_cast<int *>(map["N"]);
            fz = std::any_cast<int *>(map["Z"]);
            px = std::any_cast<double *>(map["px"]);
            py = std::any_cast<double *>(map["py"]);
            pz = std::any_cast<double *>(map["pz"]);
        }

        catch (const std::bad_any_cast &e)
        {
            std::cout << e.what() << '\n';
        }
        event event = {Nc, bimp};
        if (!this->event_cut.pass(event))
        {
            continue;
        }

        for (unsigned int i = 0; i < multi; i++)
        {
            particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            event.particles.push_back(particle);
        }
        this->hist.fill(event);
    }
}

void manager::finish()
{
    this->outputfile = new TFile(this->path_out.c_str(), "RECREATE");
    this->outputfile->cd();
    this->hist.normalize();
    this->hist.write();
    this->outputfile->Write();
}