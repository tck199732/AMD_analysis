#include "amd2root.hh"

int main(int argc, char **argv)
{
    std::string reaction = argv[1];
    std::string mode = argv[2];
    std::string path_data = argv[3];
    std::string path_out = argv[4];
    std::string path_coll_hist = "";
    std::string path_amdgid = "";

    if (mode != "21" && mode != "3" && mode != "21t")
    {
        std::cout << "input mode = " << mode << std::endl;
        throw std::invalid_argument("invalid mode : 21, 21t or 3");
    }

    if (mode == "21t")
    {
        if (argc < 7)
        {
            throw std::invalid_argument("please provide path for amdgid and coll_hist.");
        }
        path_amdgid = argv[5];
        path_coll_hist = argv[6];
    }

    manager manager = {reaction, mode, path_data, path_out, path_amdgid, path_coll_hist};
    manager.init();
    manager.compile();
    manager.finish();
    return 0;
}

void manager::init()
{

    if (!fs::exists(this->path_data))
    {
        throw std::invalid_argument("data path does not exists.");
    }

    // find the total number of nucleons in the reaction system
    // e.g. Ca48Ni64E140 -> 48 + 64 = 112
    std::regex regexp("[0-9]+");
    std::smatch matches;
    this->amass = 0;

    std::string temp = this->reaction;
    for (auto _ = 0; _ < 2; _++)
    {
        std::regex_search(temp, matches, regexp);
        this->amass += stoi(matches[0]);
        temp = matches.suffix().str();
    }
    std::cout << "reaction system : " << this->reaction << std::endl;
    std::cout << "Nucleon Mass of Reaction : " << this->amass << std::endl;
    std::cout << "Operation Mode : " << this->mode << std::endl;
    std::cout << "path to data : " << fs::absolute(this->path_data) << std::endl;
    std::cout << "path to output : " << fs::absolute(this->path_out) << std::endl;

    if (this->mode == "21t")
    {
        if (!fs::exists(this->path_amdgid))
        {
            throw std::invalid_argument("amdgid.dat does not exist.");
        }
        if (!fs::exists(this->path_coll_hist))
        {
            throw std::invalid_argument("coll_hist.dat does not exist.");
        }
        std::cout << "path to amdgid : " << fs::absolute(this->path_amdgid) << std::endl;
        std::cout << "path to coll hist : " << fs::absolute(this->path_coll_hist) << std::endl;
    }

    this->init_writer();
}

void manager::init_writer()
{
    this->writer = new RootWriter(fs::absolute(this->path_out), this->tr_name, "RECREATE");

    std::vector<branch> branches = {
        {"multi", "int"},
        {"b", "double"},
        {"N", "int[multi]"},
        {"Z", "int[multi]"},
        {"px", "double[multi]"},
        {"py", "double[multi]"},
        {"pz", "double[multi]"},
    };

    if (this->mode == "21")
    {
        branches.push_back({"ENG", "double[multi]"});
        branches.push_back({"LANG", "double[multi]"});
        branches.push_back({"JX", "double[multi]"});
        branches.push_back({"JY", "double[multi]"});
        branches.push_back({"JZ", "double[multi]"});
    }

    if (this->mode == "3")
    {
        branches.push_back({"J", "double[multi]"});
        branches.push_back({"M", "double[multi]"});
        branches.push_back({"WEIGHT", "double[multi]"});
        branches.push_back({"iFRG", "int[multi]"});
    }

    if (this->mode == "21t")
    {
        branches.push_back({"t", "double[multi]"});
        branches.push_back({"x", "double[multi]"});
        branches.push_back({"y", "double[multi]"});
        branches.push_back({"z", "double[multi]"});
    }

    for (auto &br : branches)
    {
        br.autofill();
    }
    this->writer->set_branches(tr_name, branches);
}

void manager::compile()
{
    if (this->mode == "21")
    {
        std::cout << "extracting table21 to root file" << std::endl;
        this->compile21();
    }
    else if (this->mode == "21t")
    {
        std::cout << "extracting table21, with xyzt to root file" << std::endl;
        this->compile21t();
    }
    else if (this->mode == "3")
    {
        std::cout << "extracting table3 to root file" << std::endl;
        this->compile3();
    }
    return;
}

void manager::compile21()
{
    std::cout << "loading table21 : " << this->path_data << std::endl;
    std::ifstream file_table21(this->path_data.c_str());

    double px, py, pz;
    double ENG, LANG, JX, JY, JZ;
    double bimp;
    int N, Z;
    int eventID;
    int nucleons_count = 0;

    std::vector<particle> particles;
    event event;
    while (!file_table21.eof())
    {
        file_table21 >> Z >> N >> px >> py >> pz;
        if (Z == 0 && N == 0)
        {
            break;
        }
        file_table21 >> ENG >> LANG >> JX >> JY >> JZ;
        file_table21 >> bimp >> eventID;

        particle particle = {Z, N, px, py, pz};
        // particle.report();
        particles.push_back(particle);
        nucleons_count += Z + N;

        if (nucleons_count == this->amass)
        {
            event = {eventID, bimp, particles};
            this->fill(event);
            particles.clear();
            nucleons_count = 0;
        }
    }
    return;
}

void manager::compile21t()
{
    int event_processed = 0;
    std::cout << "loading table21 : " << fs::absolute(this->path_data) << std::endl;
    std::ifstream file_table21(fs::absolute(this->path_data));

    std::cout << "loading amdgid : " << fs::absolute(this->path_amdgid) << std::endl;
    std::ifstream file_amdgid(fs::absolute(this->path_amdgid));
    file_amdgid.ignore(99, '\n');

    std::cout << "loading collision history : " << fs::absolute(this->path_coll_hist) << std::endl;
    std::ifstream file_coll_hist(fs::absolute(this->path_coll_hist));
    file_coll_hist.ignore(99, '\n');

    // particle id -> [gids]
    std::map<int, std::vector<int>> gidmap;

    // gid -> [gid2, [t,x,y,z,px,py,pz]]
    std::map<int, std::vector<std::pair<int, std::vector<double>>>> coll_hist;

    int prim_pid, nuc, gid, N, Z, ievt;

    while (!file_amdgid.eof())
    {
        for (int _ = 1; _ <= this->amass; _++)
        {
            file_amdgid >> prim_pid >> nuc >> gid >> N >> Z >> ievt;
            // std::cout << prim_pid << "\t" << nuc << "\t" << gid << "\t" << N << "\t" << Z << "\t" << ievt << std::endl;
            if (gidmap.count(prim_pid) == 0)
            {
                gidmap[prim_pid] = std::vector<int>(1, gid);
            }
            else
            {
                gidmap[prim_pid].push_back(gid);
            }
        }

        // std::cout << gidmap.size() << std::endl;
        int current_evt = ievt;
        while (!file_coll_hist.eof())
        {
            int collid;
            std::vector<double> rp(7, 0.0);
            file_coll_hist >> ievt >> gid >> collid;
            for (int r = 0; r < 7; r++)
            {
                file_coll_hist >> rp[r];
            }
            std::pair<int, std::vector<double>> pair = std::make_pair(collid, rp);
            if (ievt == current_evt)
            {
                if (coll_hist.count(gid) == 0)
                {
                    std::vector<std::pair<int, std::vector<double>>> hist;
                    coll_hist[gid] = hist;
                }
                coll_hist[gid].push_back(pair);
            }
            else
            {
                std::vector<particle> particles;
                event event;
                double px, py, pz, bimp, buffer;
                for (unsigned int j = 1; j <= gidmap.size(); j++)
                {
                    std::vector<double> spacetime(4, 0.);
                    spacetime[0] = -1.;
                    for (const auto &id : gidmap[j])
                    {
                        std::vector<double> last_interaction(7, 0.); // for each nucleon in a prim. particle
                        for (const auto &hist : coll_hist[id])
                        {
                            int collid = hist.first;
                            std::vector<double> h = hist.second;
                            // if the interaction is not between two nucleons in the same prim. particle
                            if (std::find(gidmap[j].begin(), gidmap[j].end(), collid) - gidmap[j].begin() == gidmap[j].size())
                            {
                                if (h[0] >= last_interaction[0])
                                {
                                    last_interaction = h;
                                }
                            }
                        }
                        spacetime[1] += last_interaction[1];
                        spacetime[2] += last_interaction[2];
                        spacetime[3] += last_interaction[3];
                        spacetime[0] = std::max(spacetime[0], last_interaction[0]);
                    }
                    spacetime[1] /= gidmap[j].size();
                    spacetime[2] /= gidmap[j].size();
                    spacetime[3] /= gidmap[j].size();

                    file_table21 >> Z >> N >> px >> py >> pz;
                    for (auto _ = 0; _ < 5; _++)
                    {
                        file_table21 >> buffer;
                    }
                    file_table21 >> bimp >> ievt;

                    particle particle = {Z, N, px, py, pz, spacetime[0], spacetime[1], spacetime[2], spacetime[3]};
                    particles.push_back(particle);
                }

                event = {ievt, bimp, particles};
                // event.report();
                this->fill(event);
                particles.clear();

                gidmap.clear();
                coll_hist.clear();
                coll_hist[gid].push_back(pair);
                break;
            }
        }

        event_processed++;
        if (event_processed % 1000 == 0)
        {
            std::cout << event_processed << std::endl;
        }
    }
}

void manager::compile3()
{
    std::cout << "loading table3 : " << this->path_data << std::endl;
    std::ifstream file_table3(this->path_data.c_str());
    file_table3.ignore(99, '\n');

    double px, py, pz;
    double J, M, WEIGHT;
    double bimp;
    int N, Z, iFRG;
    int eventID;
    int nucleons_count = 0;

    std::vector<particle> particles;
    event event;
    while (!file_table3.eof())
    {
        file_table3 >> Z >> N >> px >> py >> pz;
        if (Z == 0 && N == 0)
        {
            break;
        }
        file_table3 >> J >> M >> WEIGHT;
        file_table3 >> bimp >> eventID;
        file_table3 >> iFRG;

        particle particle = {Z, N, px, py, pz};
        particles.push_back(particle);
        nucleons_count += Z + N;
        if (nucleons_count == this->amass)
        {
            event = {eventID, bimp, particles};
            this->fill(event);
            particles.clear();
            nucleons_count = 0;
        }
    }
    return;
}

void manager::fill(const event &event)
{
    std::vector<double> px, py, pz;
    std::vector<int> N, Z;
    std::vector<double> t, x, y, z;
    auto resize = [event](auto &&...args)
    {
        (args.resize(event.particles.size()), ...);
    };
    resize(px, py, pz, N, Z, t, x, y, z);
    int id = 0;
    for (const auto &particle : event.particles)
    {
        px[id] = particle.px;
        py[id] = particle.py;
        pz[id] = particle.pz;
        N[id] = particle.N;
        Z[id] = particle.Z;
        if (this->mode == "21t")
        {
            t[id] = particle.t;
            x[id] = particle.x;
            y[id] = particle.y;
            z[id] = particle.z;
        }
        id++;
    }

    // std::cout << sizeof(px[0]) * px.size() << std::endl;

    this->writer->set(this->tr_name, "multi", (int)event.particles.size());
    this->writer->set(this->tr_name, "b", event.b);
    this->writer->set(this->tr_name, "px", px);
    this->writer->set(this->tr_name, "py", py);
    this->writer->set(this->tr_name, "pz", pz);
    this->writer->set(this->tr_name, "N", N);
    this->writer->set(this->tr_name, "Z", Z);

    if (this->mode == "21t")
    {
        this->writer->set(this->tr_name, "t", t);
        this->writer->set(this->tr_name, "x", x);
        this->writer->set(this->tr_name, "y", y);
        this->writer->set(this->tr_name, "z", z);
    }

    this->writer->fill();
}

void manager::finish()
{
    this->writer->write();
    std::cout << "DONE" << std::endl;
}

void particle::report()
{
    auto print = [](auto &&...args)
    {
        ((std::cout << args << "\t"), ...) << std::endl;
    };
    print(this->Z, this->N, this->px, this->py, this->pz);
}

void event::report()
{

    std::cout << "event : " << this->eventID << '\t'
              << "b  = " << this->b << "[fm]" << std::endl;
    for (auto &particle : this->particles)
    {
        particle.report();
    }
}