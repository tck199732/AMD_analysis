#include "system_info.hh"
void coll_system::autofill()
{
    this->name = Form(
        "%s%sE%i", this->beam.c_str(), this->target.c_str(), this->energy);

    double mass_1u = 931.49410242;
    std::vector<std::string> beam_name = {"Ca40", "Ca48"};
    std::vector<double> beam_mass = {39.962590866, 47.95252290};
    std::vector<std::string> target_name = {"Ni58", "Ni64", "Sn112", "Sn124"};
    std::vector<double> target_mass = {57.935343, 63.9279660, 111.904818, 123.9052739};

    for (unsigned int i = 0; i < beam_mass.size(); i++)
    {
        beam_mass[i] *= mass_1u;
    }
    for (unsigned int i = 0; i < target_mass.size(); i++)
    {
        target_mass[i] *= mass_1u;
    }

    std::vector<int> beam_energy = {56, 140}; // MeV/A
    std::vector<int> beam_A = {40, 48};       // MeV/A

    auto is_contain = [](auto &arr, const auto &ele) -> bool
    {
        return std::find(arr.begin(), arr.end(), ele) != arr.end();
    };

    auto get_index = [](auto &arr, const auto &ele) -> int
    {
        return std::find(arr.begin(), arr.end(), ele) - arr.begin();
    };

    if (is_contain(beam_name, this->beam) && is_contain(target_name, this->target) && is_contain(beam_energy, this->energy))
    {

        int beam_id = get_index(beam_name, this->beam);
        int target_id = get_index(target_name, this->target);
        int en_id = get_index(beam_energy, this->energy);

        double beam_ke = beam_energy[en_id] * beam_A[beam_id];
        double beam_energy_tot = beam_ke + beam_mass[beam_id];
        double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * beam_mass[beam_id]);
        double gamma = beam_energy_tot / beam_mass[beam_id];

        this->betacms = mom_beam / (gamma * beam_mass[beam_id] + target_mass[target_id]);
        this->rapidity_beam = 0.5 * TMath::Log((beam_energy_tot + mom_beam) / (beam_energy_tot - mom_beam));
    }
}

system_info::system_info()
{

    std::vector<std::string> beam_name = {"Ca40", "Ca48"};
    std::vector<std::string> target_name = {"Ni58", "Ni64", "Sn112", "Sn124"};
    std::vector<int> beam_energy = {56, 140}; // MeV/A

    for (auto beam : beam_name)
    {
        for (auto target : target_name)
        {
            for (auto en : beam_energy)
            {
                coll_system sys = {beam, target, en};
                sys.autofill();

                this->systems[sys.name] = sys;
            }
        }
    }
}

system_info::~system_info()
{
    ;
}
