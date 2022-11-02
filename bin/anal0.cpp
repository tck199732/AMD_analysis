
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
    TH1D *h1_multi;
    TH1D *h1_multi_H;
    TH1D *h1_multi_He;
    TH1D *h1_Z_prim;
    TH1D *h1_Z_seq;

    double weight = 1.;
    double norm = 0.0;
    void init();
    void fill(const event &event);
    void fill(const event &event21, const event &event3);
    void normalize();
    void write();
};

void histograms::init()
{
    this->h1_multi = new TH1D(("h1_multi_mode" + this->mode).c_str(), "", 80, 0, 80);
    this->h1_multi_H = new TH1D(("h1_multi_H_mode" + this->mode).c_str(), "", 5, 0, 5);
    this->h1_multi_He = new TH1D(("h1_multi_He_mode" + this->mode).c_str(), "", 5, 0, 5);
    this->h1_Z_prim = new TH1D(("h1_Z_prim_mode" + this->mode).c_str(), "", 50, 0, 50);

    this->h1_multi->Sumw2();
    this->h1_multi_H->Sumw2();
    this->h1_multi_He->Sumw2();
    this->h1_Z_prim->Sumw2();
    this->h1_multi->SetDirectory(0);
    this->h1_multi_H->SetDirectory(0);
    this->h1_multi_He->SetDirectory(0);
    this->h1_Z_prim->SetDirectory(0);

    if (this->mode == "3")
    {
        this->weight = 1. / ndecays;
        this->h1_Z_seq = new TH1D(("h1_Z_seq_mode" + this->mode).c_str(), "", 50, 0, 50);
        this->h1_Z_seq->Sumw2();
        this->h1_Z_seq->SetDirectory(0);
    }
}

void histograms::fill(const event &event)
{
    this->h1_multi->Fill(event.multi, this->weight);

    for (auto &par : event.particles)
    {
        int zid = par.aid - par.nid;
        if (zid == 1)
        {
            this->h1_multi_H->Fill(par.aid, this->weight);
        }
        if (zid == 2)
        {
            this->h1_multi_He->Fill(par.aid, this->weight);
        }
        if (this->mode == "21")
        {
            this->h1_Z_prim->Fill(zid, this->weight);
        }
    }

    this->norm += this->weight;
}
void histograms::fill(const event &event21, const event &event3)
{
    if (this->weight != 1. / ndecays)
    {
        throw std::invalid_argument("should not use this fill function");
    }
    this->fill(event3);

    for (auto &par3 : event3.particles)
    {
        bool isPrim = 0;
        for (auto &par21 : event21.particles)
        {
            // if it is the same particle
            if (par3.nid == par21.nid && par3.aid == par21.aid)
            {
                if (pow(par3.px - par21.px, 2.) + pow(par3.py - par21.py, 2.) + pow(par3.pz - par21.pz, 2.) <= 0.001)
                {
                    this->h1_Z_prim->Fill(par3.aid - par3.nid, this->weight);
                    isPrim = 1;
                }
            }
        }
        if (!isPrim)
        {
            this->h1_Z_seq->Fill(par3.aid - par3.nid, this->weight);
        }
    }
}

void histograms::normalize()
{
    this->h1_multi->Scale(1. / this->norm);
    this->h1_multi_H->Scale(1. / this->norm);
    this->h1_multi_He->Scale(1. / this->norm);
    this->h1_Z_prim->Scale(1. / this->norm);
    if (this->mode == "3")
    {
        this->h1_Z_seq->Scale(1. / this->norm);
    }
}
void histograms::write()
{
    this->h1_multi->Write();
    this->h1_multi_H->Write();
    this->h1_multi_He->Write();
    this->h1_Z_prim->Write();
    if (this->mode == "3")
    {
        this->h1_Z_seq->Write();
    }
}

struct manager
{
    std::string reaction;
    fs::path data_dir;
    fs::path path_list;
    fs::path path_out;

    std::vector<fs::path> data_paths21;
    std::vector<fs::path> data_paths3;
    system_info *sys_info;
    double betacms, rapidity_beam;
    RootReader *reader21;
    RootReader *reader3;
    histograms hist21;
    histograms hist3;
    void init();
    void read();
    void fill(const event &event);
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string data_dir = argv[2];
    std::string path_list = argv[3];
    std::string path_out = argv[4];

    manager manager = {reaction, data_dir, path_list, path_out};
    manager.init();
    manager.read();
    manager.finish();
    return 0;
}

void manager::init()
{
    std::string pth_name;
    std::ifstream pth_stream(fs::absolute(this->path_list));
    while (pth_stream >> pth_name)
    {
        fs::path pth21 = this->data_dir / (pth_name + "_table21.root");
        fs::path pth3 = this->data_dir / (pth_name + "_table3.root");
        std::cout << pth21.string();
        if (fs::exists(pth21))
        {
            std::cout << " exist" << std::endl;
            this->data_paths21.push_back(pth21);
        }
        else
        {
            std::cout << " does not exist" << std::endl;
        }
        std::cout << pth3.string();
        if (fs::exists(pth3))
        {
            std::cout << " exist" << std::endl;
            this->data_paths3.push_back(pth3);
        }
        else
        {
            std::cout << " does not exist" << std::endl;
        }
    };

    std::vector<branch> rbranches = {
        {"multi", "int"},
        {"b", "double"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
    };

    for (auto &br : rbranches)
    {
        br.autofill();
    }

    this->reader21 = new RootReader("AMD");
    for (auto &pth : this->data_paths21)
    {
        reader21->add_file(fs::absolute(pth));
    }
    reader21->set_branches(rbranches);

    this->reader3 = new RootReader("AMD");
    for (auto &pth : this->data_paths3)
    {
        reader3->add_file(fs::absolute(pth));
    }
    reader3->set_branches(rbranches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->hist21 = {this->reaction, "21"};
    this->hist21.init();
    this->hist3 = {this->reaction, "3"};
    this->hist3.init();
}

void manager::read()
{
    int nevents21 = this->reader21->tree->GetEntries();
    int nevents3 = this->reader3->tree->GetEntries();
    std::cout << "number of events = " << nevents21 << std::endl;

    int fmulti;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    std::map<std::string, std::any> map;
    for (int ibatch = 0; ibatch < ndecays; ibatch++)
    {
        int batch_evt0 = ibatch * nevents21;
        for (int ievt = 0; ievt < nevents21; ievt++)
        {
            map = this->reader21->get_entry(ievt);
            try
            {
                fmulti = std::any_cast<int>(map["multi"]);
                bimp = std::any_cast<double>(map["b"]);
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
            event event21 = {fmulti, bimp};
            for (unsigned int i = 0; i < fmulti; i++)
            {
                particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
                particle.autofill(this->betacms);
                event21.particles.push_back(particle);
            }
            this->hist21.fill(event21);

            int ievt3 = batch_evt0 + ievt;
            map = this->reader3->get_entry(ievt3);
            try
            {
                fmulti = std::any_cast<int>(map["multi"]);
                bimp = std::any_cast<double>(map["b"]);
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
            event event3 = {fmulti, bimp};
            for (unsigned int i = 0; i < fmulti; i++)
            {
                particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
                particle.autofill(this->betacms);
                event3.particles.push_back(particle);
            }
            this->hist3.fill(event21, event3);

            // std::cout << "reading " << ievt3 << std::endl;
        }
    }
}

void manager::finish()
{

    TFile *outf = new TFile(this->path_out.c_str(), "RECREATE");
    outf->cd();
    this->hist21.normalize();
    this->hist3.normalize();
    this->hist21.write();
    this->hist3.write();
    outf->Write();
    outf->Close();
    std::cout << "DONE" << std::endl;
}
