
#include "spectra.hh"
#include <chrono>
struct manager
{
    std::string reaction;
    fs::path path_data21;
    fs::path path_data3;
    fs::path path_out;
    int Ncmin, Ncmax;
    double bmin, bmax;
    system_info *sys_info;
    double betacms, rapidity_beam;

    RootReader *reader21;
    RootReader *reader3;
    histograms hist21;
    histograms hist3;
    histograms hist3_ndecay1;
    eventcut event_cut;
    TFile *outputfile;

    int nevents21, nevents3;
    std::chrono::duration<double> timer21, timer3;
    void init();
    void read();
    void replace_error();
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string path_data21 = argv[2];
    std::string path_data3 = argv[3];
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

    manager manager = {reaction, path_data21, path_data3, path_out, Ncmin, Ncmax, bmin, bmax};
    manager.init();
    manager.read();
    manager.replace_error();
    manager.finish();
    return 0;
}

void manager::init()
{
    if (!fs::exists(this->path_data21))
    {
        throw std::invalid_argument("path_data21 does not exist.");
    }
    if (!fs::exists(this->path_data3))
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

    this->reader21 = new RootReader("AMD");
    reader21->add_file(fs::absolute(this->path_data21));
    reader21->set_branches(branches);
    this->reader3 = new RootReader("AMD");
    reader3->add_file(fs::absolute(this->path_data3));
    reader3->set_branches(branches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->event_cut = {{this->Ncmin, this->Ncmax}, {this->bmin, this->bmax}};

    this->hist21 = {this->reaction, "prim", this->betacms, this->rapidity_beam};
    this->hist3 = {this->reaction, "seq", this->betacms, this->rapidity_beam};
    this->hist3_ndecay1 = {this->reaction, "seq1", this->betacms, this->rapidity_beam};

    this->hist21.init();
    this->hist3.init();
    this->hist3_ndecay1.init();
    this->nevents21 = this->reader21->tree->GetEntries();
    this->nevents3 = this->reader3->tree->GetEntries();
    std::cout << "number of events = " << nevents21 << std::endl;
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

    int ndecays = this->nevents3 / this->nevents21;
    std::vector<double> weights(this->nevents21, 0.);
    auto start3 = std::chrono::steady_clock::now();
    for (int ievt3 = 0; ievt3 < this->nevents3; ievt3++)
    {

        std::map<std::string, std::any> map = this->reader3->get_entry(ievt3);
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
        event event3 = {Nc, bimp};
        if (!this->event_cut.pass(event3))
        {
            continue;
        }

        weights[ievt3 % nevents21] += 1.;

        for (unsigned int i = 0; i < multi; i++)
        {
            particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            if (ievt3 < nevents3 / ndecays)
            {
                this->hist3_ndecay1.fill(particle, 1.);
                this->hist3_ndecay1.norm += 1.;
            }
            this->hist3.fill(particle, 1. / ndecays);
            this->hist3.norm += 1. / ndecays;
        }
    }
    auto end3 = std::chrono::steady_clock::now();
    this->timer3 = end3 - start3;

    auto start21 = std::chrono::steady_clock::now();
    for (int ievt21 = 0; ievt21 < nevents21; ievt21++)
    {
        weights[ievt21] /= 1.0 * ndecays;
        std::map<std::string, std::any> map = this->reader21->get_entry(ievt21);
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

        event event21 = {Nc, bimp};
        // if (!this->event_cut.pass(event21))
        // {
        //     continue;
        // }

        for (unsigned int i = 0; i < multi; i++)
        {
            particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            this->hist21.fill(particle, weights[ievt21]);
        }

        this->hist21.norm += weights[ievt21];
    }
    auto end21 = std::chrono::steady_clock::now();
    this->timer21 = end21 - start21;
}

void manager::replace_error()
{
    for (auto &[pn, h2] : this->hist3.h2_pta_rapidity_lab)
    {
        for (int i = 0; i < h2->GetNbinsX(); i++)
        {
            for (int j = 0; j < h2->GetNbinsY(); j++)
            {
                double error = this->hist3_ndecay1.h2_pta_rapidity_lab[pn]->GetBinError(i + 1, j + 1);
                h2->SetBinError(i + 1, j + 1, error);
            }
        }
    }
}

void manager::finish()
{
    this->outputfile = new TFile(this->path_out.c_str(), "RECREATE");
    this->outputfile->cd();
    this->hist21.normalize();
    this->hist3.normalize();
    this->hist3_ndecay1.normalize();
    this->hist21.write();
    this->hist3.write();
    this->hist3_ndecay1.write();
    this->outputfile->Write();
    std::cout << "Elapsed Time reading table3 : " << this->timer3.count() << std::endl;
    std::cout << "Elapsed Time reading table21 : " << this->timer21.count() << std::endl;
    std::cout << "DONE" << std::endl;
}