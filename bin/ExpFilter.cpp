
#include "ExpFilter.hh"
fs::path DIR_DATA = "/data/AMD/root_data/dec2021";

struct manager
{
    std::string reaction, mode;
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
    std::string path_list = argv[3];
    std::string path_out = argv[4];

    manager manager = {reaction, mode, path_list, path_out};
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

    this->writer = new RootWriter(fs::absolute(this->path_out), "AMD");
    this->writer->set_branches("AMD", branches);

    this->sys_info = new system_info();
    this->betacms = sys_info->get_betacms(this->reaction);
    this->rapidity_beam = sys_info->get_rapidity_beam(this->reaction);

    this->det_uball.init();
    this->det_uball.config(this->reaction);
}

void manager::read()
{
    int nevents = this->reader->tree->GetEntries();
    std::cout << "number of events = " << nevents << std::endl;
    std::vector<particle> hira_particles;

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

        this->det_uball.reset();

        for (unsigned int i = 0; i < fmulti; i++)
        {
            particle particle{fn[i], fz[i], px[i], py[i], pz[i]};
            particle.autofill(this->betacms);
            this->det_uball.read_particle(particle);

            if (this->det_hira.pass_angle(particle) && this->det_hira.pass_ekinlab(particle))
            {
                hira_particles.push_back(particle);
            }
        }

        int Nc = this->det_uball.counter;
        event event = {Nc, bimp, hira_particles};

        this->fill(event);
        hira_particles.clear();
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

    int id = 0;
    for (const auto &particle : event.particles)
    {
        px[id] = particle.px;
        py[id] = particle.py;
        pz[id] = particle.pz;
        N[id] = particle.nid;
        Z[id] = particle.zid;
        id++;
    }

    this->writer->set("AMD", "multi", (int)event.particles.size());
    this->writer->set("AMD", "BSim", event.bimp);
    this->writer->set("AMD", "pxCMS", px);
    this->writer->set("AMD", "pyCMS", py);
    this->writer->set("AMD", "pzCMS", pz);
    this->writer->set("AMD", "N", N);
    this->writer->set("AMD", "Z", Z);
    this->writer->fill();
}

void manager::finish()
{
    this->writer->write();
    std::cout << "DONE" << std::endl;
}
