
#include "spacetime.hh"
#include <chrono>
struct manager
{
    std::string reaction;
    fs::path path_data;
    fs::path path_out;
    system_info *sys_info;
    double betacms, rapidity_beam;

    RootReader *reader;
    histograms hist;
    TFile *outputfile;

    int nevents;
    std::chrono::duration<double> timer;
    void init();
    void read();
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string path_data = argv[2];
    std::string path_out = argv[3];

    manager manager = {reaction, path_data, path_out};
    manager.init();
    manager.read();
    manager.finish();
    return 0;
}

void manager::init()
{
    if (!fs::exists(this->path_data))
    {
        throw std::invalid_argument("path_data does not exist.");
    }

    std::vector<branch> branches = {
        {"multi", "int"},
        {"b", "double"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
        {"x", "double[multi]"},
        {"y", "double[multi]"},
        {"z", "double[multi]"},
        {"t", "double[multi]"},
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

    this->hist = {this->reaction, this->betacms, this->rapidity_beam};

    this->hist.init();
    this->nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;
}

void manager::read()
{
    double count_neg_time = 0.;

    int multi;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;
    double *x;
    double *y;
    double *z;
    double *t;

    for (int ievt = 0; ievt < nevents; ievt++)
    {
        std::map<std::string, std::any> map = this->reader->get_entry(ievt);
        try
        {
            bimp = std::any_cast<double>(map["b"]);
            multi = std::any_cast<int>(map["multi"]);
            fn = std::any_cast<int *>(map["N"]);
            fz = std::any_cast<int *>(map["Z"]);
            px = std::any_cast<double *>(map["px"]);
            py = std::any_cast<double *>(map["py"]);
            pz = std::any_cast<double *>(map["pz"]);
            x = std::any_cast<double *>(map["x"]);
            y = std::any_cast<double *>(map["y"]);
            z = std::any_cast<double *>(map["z"]);
            t = std::any_cast<double *>(map["t"]);
        }

        catch (const std::bad_any_cast &e)
        {
            std::cout << e.what() << '\n';
        }

        for (unsigned int i = 0; i < multi; i++)
        {
            particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
            particle.set_xyzt(x[i], y[i], z[i], t[i]);
            particle.autofill(this->betacms);

            if (particle.t <= 0.0)
            {
                count_neg_time += 1.;
                std::cout << particle.t << std::endl;
            }
            this->hist.fill(particle, 1.);
        }

        this->hist.norm += 1.;
    }
    std::cout << count_neg_time << std::endl;
}

void manager::finish()
{
    this->outputfile = new TFile(this->path_out.c_str(), "RECREATE");
    this->outputfile->cd();
    this->hist.normalize();
    this->hist.write();
    this->outputfile->Write();
}