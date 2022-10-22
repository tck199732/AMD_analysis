

#include "ExpFilter.hh"
#include "TFile.h"
#include "TH2D.h"
struct histograms
{
    std::string reaction, mode;
    std::vector<std::string> particlenames;

    std::map<std::string, TH2D *> h2_pta_rapidity_lab_all;

    void init();
    void fill(const event &event);
    void write();
};

void histograms::init()
{

    for (const auto &pn : this->particlenames)
    {
        this->
    }
}
void histograms::write()
{
    this->h2_multi_b->Write();
}

void histograms::fill(const event &event)
{
    this->h2_multi_b->Fill(event.Nc, event.bimp);
}

struct manager
{
    std::string reaction, mode;
    fs::path path_list;
    fs::path path_out;

    std::vector<fs::path> data_paths;
    system_info *sys_info;
    double betacms, rapidity_beam;

    RootReader *reader;
    histograms histograms;
    TFile *outputfile;

    void init();
    void read();
    void fill(const event &event);
    void finish();
};

void main(int argc, char *argv[])
{

    std::string reaction = argv[1];
    std::string mode = argv[2];
    std::string path_list = argv[3];
    std::string path_out = argv[4];
    manager manager = {reaction, mode, path_list, path_out};
    manager.init();
    manager.read();
    manager.finish();
}

void manager::init()
{
    std::string pth_name;
    std::ifstream pth_stream(fs::absolute(this->path_list));
    while (pth_stream >> pth_name)
    {
        fs::path pth = DIR_DATA / (pth_name + "_table" + this->mode + ".root");
        std::cout << pth.string();
        if (fs::exists(pth))
        {
            std::cout << " exist" << std::endl;
            this->data_paths.push_back(pth);
        }
        else
        {
            std::cout << " does not exist" << std::endl;
        }
    };

    std::vector<branch> branches = {
        {"multi", "int"},
        {"BSim", "double"},
        {"pxCMS", "double[]"},
        {"pyCMS", "double[]"},
        {"pzCMS", "double[]"},
        {"N", "int[]"},
        {"Z", "int[]"},
    };

    for (auto &br : branches)
    {
        br.autofill();
    }

    this->reader = new RootReader("AMD");
    for (auto &pth : this->data_paths)
    {
        reader->add_file(fs::absolute(pth));
    }
    reader->set_branches(branches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    histograms.init();
}

void manager::read()
{
    int nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;

    int Nc;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    std::vector<particle> particles;
    for (int ievt = 0; ievt < nevents; ievt++)
    {
        std::map<std::string, std::any> map = this->reader->get_entry(ievt);
        try
        {
            Nc = std::any_cast<int>(map["multi"]);
            bimp = std::any_cast<double>(map["BSim"]);
            fn = std::any_cast<int *>(map["N"]);
            fz = std::any_cast<int *>(map["Z"]);
            px = std::any_cast<double *>(map["pxCMS"]);
            py = std::any_cast<double *>(map["pyCMS"]);
            pz = std::any_cast<double *>(map["pzCMS"]);
        }

        catch (const std::bad_any_cast &e)
        {
            std::cout << e.what() << '\n';
        }

        for (unsigned int i = 0; i < ; i++)
        {
            particle particle{fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            particles.push_back(particle);
        }

        event event = {Nc, bimp, particles};

        this->fill(event);
        particles.clear();
    }
}

void manager::fill(const event &event)
{
    this->histograms.fill(event);
}

void manager::finish()
{
    this->outputfile = new TFile(fs::absolute(this->path_out), "RECREATE");
    this->histograms->write();
    this->outputfile->write();
    std::cout << "DONE" << std::endl;
}