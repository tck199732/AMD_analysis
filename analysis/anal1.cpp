#include "anal1.hh"

struct manager
{
    std::string reaction;
    fs::path data_dir;
    fs::path path_list;
    fs::path path_out;

    std::vector<fs::path> data_paths3;
    std::map<std::string, RootReader *> reader3;

    system_info *sys_info;
    double betacms, rapidity_beam;

    histograms hist3;
    void init();
    void read();
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
        fs::path pth3 = this->data_dir / (pth_name + "_table3.root");

        if (fs::exists(pth3))
        {
            this->reader3[pth_name] = new RootReader("AMD");
            this->reader3[pth_name]->add_file(fs::absolute(pth3));
        }
        else
        {
            std::cout << pth3.string() << " does not exist" << std::endl;
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

    for (auto &[pth, reader] : this->reader3)
    {
        this->reader3[pth]->set_branches(rbranches);
    }

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->hist3 = {this->reaction, "3"};
    this->hist3.init();

    for (auto &[pth, reader] : this->reader3)
    {
        int nevents3 = this->reader3[pth]->tree->GetEntries();
        std::cout << "number of events3 = " << nevents3 << std::endl;
    }
}

void manager::read()
{
    int fmulti;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    std::map<std::string, std::any> map;

    for (auto &[pth, reader] : this->reader3)
    {
        for (int ievt3 = 0; ievt3 < this->reader3[pth]->tree->GetEntries(); ievt3++)
        {
            map = this->reader3[pth]->get_entry(ievt3);
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
            this->hist3.fill(event3);
        }
    }
}

void manager::finish()
{
    TFile *outf = new TFile(this->path_out.c_str(), "RECREATE");
    outf->cd();
    this->hist3.normalize();
    this->hist3.write();
    outf->Write();
    outf->Close();
    std::cout << "DONE" << std::endl;
}
