
#include "ExpFilter.hh"

struct manager
{
    std::string reaction, mode;
    fs::path data_dir;
    fs::path path_list;
    fs::path path_out;

    std::vector<fs::path> data_paths;

    system_info *sys_info;
    double betacms, rapidity_beam;
    hira det_hira;
    uball det_uball;

    RootWriter *writer;
    RootReader *reader;

    void init();
    void read();
    void fill(const event &event);
    void finish();
};

int main(int argc, char *argv[])
{
    std::string reaction = argv[1];
    std::string mode = argv[2];
    std::string data_dir = argv[3];
    std::string path_list = argv[4];
    std::string path_out = argv[5];

    manager manager = {reaction, mode, data_dir, path_list, path_out};
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
        fs::path pth = this->data_dir / (pth_name + "_table" + this->mode + ".root");
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

    this->reader = new RootReader("AMD");
    for (auto &pth : this->data_paths)
    {
        reader->add_file(fs::absolute(pth));
    }
    reader->set_branches(rbranches);

    this->writer = new RootWriter(fs::absolute(this->path_out), "AMD");

    std::vector<branch> wbranches = {
        {"multi", "int"},
        {"b", "double"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
        {"Nc", "int"},
    };

    for (auto &br : wbranches)
    {
        br.autofill();
    }
    this->writer->set_branches("AMD", wbranches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->det_hira.init();

    for (auto &[pn, cut] : det_hira.Ekinlabcut)
    {
        std::cout << pn << "\t" << cut[0] << "\t" << cut[1] << std::endl;
    }

    this->det_uball.init();
    this->det_uball.config(this->reaction);
}

void manager::read()
{
    int nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;

    int fmulti;
    double bimp;
    int *fz;
    int *fn;
    double *px;
    double *py;
    double *pz;

    for (int ievt = 0; ievt < nevents; ievt++)
    {
        std::map<std::string, std::any> map = this->reader->get_entry(ievt);
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

        this->det_uball.reset();

        event event;
        for (unsigned int i = 0; i < fmulti; i++)
        {
            particle particle = {fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            this->det_uball.read_particle(particle);
            if (this->det_hira.pass_angle(particle) && this->det_hira.pass_ekinlab(particle))
            {
                event.particles.push_back(particle);
            }
        }

        event.Nc = this->det_uball.counter;
        event.bimp = bimp;
        this->fill(event);
    }
}

void manager::fill(const event &event)
{
    std::vector<double> px, py, pz;
    std::vector<int> N, Z;

    auto resize = [event](auto &&...args)
    {
        (args.resize(event.particles.size()), ...);
    };
    resize(px, py, pz, N, Z);

    for (unsigned int i = 0; i < event.particles.size(); i++)
    {
        px[i] = event.particles[i].px;
        py[i] = event.particles[i].py;
        pz[i] = event.particles[i].pz;
        N[i] = event.particles[i].nid;
        Z[i] = event.particles[i].zid;
    }

    this->writer->set("AMD", "Nc", event.Nc);
    this->writer->set("AMD", "multi", (int)event.particles.size());
    this->writer->set("AMD", "b", event.bimp);
    this->writer->set("AMD", "px", px);
    this->writer->set("AMD", "py", py);
    this->writer->set("AMD", "pz", pz);
    this->writer->set("AMD", "N", N);
    this->writer->set("AMD", "Z", Z);

    this->writer->fill("AMD");
}

void manager::finish()
{
    this->writer->write();
    std::cout << "DONE" << std::endl;
}
