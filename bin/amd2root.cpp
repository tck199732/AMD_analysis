#include "amd2root.hh"

void _check_path(const std::string &path)
{
    if (!fs::exists(path))
    {
        std::string msg = Form("%s does not exists.", path.c_str());
        throw std::invalid_argument(msg.c_str());
    }
    return;
}

int get_number_nucleons(const std::string reaction)
{
    std::regex regexp("[0-9]+");
    std::smatch matches;
    int amass = 0;

    std::string temp = reaction;
    for (auto _ = 0; _ < 2; _++)
    {
        std::regex_search(temp, matches, regexp);
        amass += stoi(matches[0]);
        temp = matches.suffix().str();
    }

    return amass;
}

void CompileTable21(TTree *&tree, const std::string &path, const int &amass);
void CompileTable3(TTree *&tree, const std::string &path, const int &amass);
void CompileTable21t(TTree *&tree, const std::string &path_table21, const std::string &path_amdgid, const std::string &path_coll_hist, const int &amass);

int main(int argc, char **argv)
{
    std::string reaction = argv[1];
    std::string mode = argv[2];

    // check arguments validity
    if (mode != "21" && mode != "3" && mode != "21t")
    {
        std::cout << "input mode = " << mode << std::endl;
        throw std::invalid_argument("acceptable modes : 21, 21t or 3");
    }

    std::string path_data = argv[3];
    std::string path_out = argv[4];
    std::string path_coll_hist = "";
    std::string path_amdgid = "";

    _check_path(path_data);
    if (mode == "21t")
    {
        path_amdgid = argv[5];
        path_coll_hist = argv[6];
        _check_path(path_amdgid);
        _check_path(path_coll_hist);
    }

    TTree *tree = new TTree("AMD", "AMD");
    Initialize_Tree(tree, mode);

    // find the total number of nucleons in the reaction system
    // e.g. Ca48Ni64E140 -> 48 + 64 = 112
    int amass = get_number_nucleons(reaction);

    if (mode == "21")
    {
        std::cout << "extracting table21 to root file" << std::endl;
        CompileTable21(tree, path_data, amass);
    }
    else if (mode == "21t")
    {
        std::cout << "extracting table21, with xyzt to root file" << std::endl;
        CompileTable21t(tree, path_data, path_amdgid, path_coll_hist, amass);
    }
    else if (mode == "3")
    {
        std::cout << "extracting table3 to root file" << std::endl;
        CompileTable3(tree, path_data, amass);
    }

    TFile *outputfile = new TFile(path_out.c_str(), "RECREATE");
    outputfile->cd();
    tree->Write();
    outputfile->Write();
    outputfile->Close();
}

void CompileTable21(TTree *&tree, const std::string &path, const int &amass)
{
    std::ifstream file_table21(path.c_str());

    int eventID;
    int nucleons_count = 0;
    int multi = 0;

    while (!file_table21.eof())
    {
        file_table21 >> amd.Z[multi] >> amd.N[multi] >> amd.px[multi] >> amd.py[multi] >> amd.pz[multi];
        if (amd.Z[multi] == 0 && amd.N[multi] == 0)
        {
            break;
        }
        file_table21 >> amd.ENG[multi] >> amd.LANG[multi] >> amd.JX[multi] >> amd.JY[multi] >> amd.JZ[multi];
        file_table21 >> amd.b >> eventID;

        nucleons_count += amd.Z[multi] + amd.N[multi];
        multi++;
        if (nucleons_count == amass)
        {
            amd.multi = multi;
            tree->Fill();
            multi = 0;
            nucleons_count = 0;
        }
    }
    return;
}

void CompileTable21t(TTree *&tree, const std::string &path_table21, const std::string &path_amdgid, const std::string &path_coll_hist, const int &amass)
{
    int event_processed = 0;
    std::ifstream file_table21(path_table21.c_str());

    std::ifstream file_amdgid(path_amdgid.c_str());
    file_amdgid.ignore(99, '\n');

    std::ifstream file_coll_hist(path_coll_hist.c_str());
    file_coll_hist.ignore(99, '\n');

    // particle id -> [gids]
    std::map<int, std::vector<int>> gidmap;

    // gid -> [gid2, [t,x,y,z,px,py,pz]]
    std::map<int, std::vector<std::pair<int, std::vector<double>>>> coll_hist;

    int prim_pid, nuc, gid, N, Z, ievt;

    while (!file_amdgid.eof())
    {
        for (int _ = 1; _ <= amass; _++)
        {
            file_amdgid >> prim_pid >> nuc >> gid >> N >> Z >> ievt;
            if (gidmap.count(prim_pid) == 0)
            {
                gidmap[prim_pid] = std::vector<int>(1, gid);
            }
            else
            {
                gidmap[prim_pid].push_back(gid);
            }
        }

        int current_evt = ievt;
        int multi = 0;
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
                double px, py, pz, bimp, buffer;
                for (unsigned int j = 1; j <= gidmap.size(); j++)
                {
                    std::vector<double> spacetime(4, 0.);
                    spacetime[0] = -1.;
                    for (const auto &id : gidmap[j])
                    {
                        // for each multi in a prim. particle
                        std::vector<double> last_interaction(7, 0.);
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

                    amd.t[multi] = spacetime[0];
                    amd.x[multi] = spacetime[1];
                    amd.y[multi] = spacetime[2];
                    amd.z[multi] = spacetime[3];

                    file_table21 >> amd.Z[multi] >> amd.N[multi] >> amd.px[multi] >> amd.py[multi] >> amd.pz[multi];
                    for (auto _ = 0; _ < 5; _++)
                    {
                        file_table21 >> buffer;
                    }
                    file_table21 >> amd.b >> ievt;
                }

                tree->Fill();
                gidmap.clear();
                coll_hist.clear();
                coll_hist[gid].push_back(pair);
                break;
            }
        }
        multi = 0;
        event_processed++;
    }
}

void CompileTable3(TTree *&tree, const std::string &path, const int &amass)
{
    std::ifstream file_table3(path.c_str());
    file_table3.ignore(99, '\n');

    int eventID;
    int nucleons_count = 0;
    int multi = 0;
    while (!file_table3.eof())
    {
        file_table3 >> amd.Z[multi] >> amd.N[multi] >> amd.px[multi] >> amd.py[multi] >> amd.pz[multi];
        if (amd.Z[multi] == 0 && amd.N[multi] == 0)
        {
            break;
        }
        file_table3 >> amd.J[multi] >> amd.M[multi] >> amd.WEIGHT[multi];
        file_table3 >> amd.b >> eventID;
        file_table3 >> amd.iFRG[multi];

        nucleons_count += amd.Z[multi] + amd.N[multi];
        multi++;

        if (nucleons_count == amass)
        {
            tree->Fill();
            nucleons_count = 0;
            multi = 0;
        }
    }
    return;
}

// void particle::report()
// {
//     auto print = [](auto &&...args)
//     {
//         ((std::cout << args << "\t"), ...) << std::endl;
//     };
//     print(this->Z, this->N, this->px, this->py, this->pz);
// }

// void event::report()
// {

//     std::cout << "event : " << this->eventID << '\t'
//               << "b  = " << this->b << "[fm]" << std::endl;
//     for (auto &particle : this->particles)
//     {
//         particle.report();
//     }
// }