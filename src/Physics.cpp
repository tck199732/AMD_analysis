#include "Physics.hh"

namespace Physics
{
    std::map<std::string, double> particle_mass_old = {
        // light clusters
        {"n", 939.2},
        {"p", 938.272},
        {"d", 1875.610},
        {"t", 2808.9182},
        {"3He", 2808.3870},
        {"4He", 3727.374},
        // nuclei
        {"Ca40", 39.962590866},
        {"Ca48", 47.95252290},
        {"Ni58", 57.935343},
        {"Ni64", 63.9279660},
        {"Sn112", 111.904818},
        {"Sn124", 123.9052739},
    };

    // from ame2020
    std::map<std::string, double> particle_mass = {
        // light clusters
        {"n", 939.5654204771444},
        {"p", 938.7830734811444},
        {"d", 1876.1239277292887},
        {"t", 2809.432118151433},
        {"3He", 2809.4135261314327},
        {"4He", 3728.4013255385776},
        // nuclei
        {"Ca40", 37224.91769468577},
        {"Ca48", 44667.49204802292},
        {"Ni58", 53966.42906919436},
        {"Ni64", 59548.52352069724},
        {"Sn112", 104238.68442172016},
        {"Sn124", 115417.03722472588},
    };

    std::map<std::string, int> particle_A = {
        {"n", 1},
        {"p", 1},
        {"d", 2},
        {"t", 3},
        {"3He", 3},
        {"4He", 4},
        {"Ca40", 40},
        {"Ca48", 48},
        {"Ni58", 58},
        {"Ni64", 64},
        {"Sn112", 112},
        {"Sn124", 124},
    };

    std::map<std::string, int> particle_Z = {
        {"n", 0},
        {"p", 1},
        {"d", 1},
        {"t", 1},
        {"3He", 2},
        {"4He", 2},
        {"Ca40", 20},
        {"Ca48", 20},
        {"Ni58", 28},
        {"Ni64", 28},
        {"Sn112", 50},
        {"Sn124", 50},
    };
    double amu_MeV = 931.49410242;
};

std::string Physics::GetNucleiName(const int &Z, const int &A)
{
    for (auto &[pn, a] : Physics::particle_A)
    {
        if (A == a && Z == Physics::particle_Z[pn])
        {
            return pn;
        }
    }
    return "no-such-particle";
}

double Physics::GetNucleiMass(const std::string &symbol, const std::string &option)
{
    return option == "new" ? Physics::particle_mass[symbol] : Physics::particle_mass_old[symbol] * Physics::amu_MeV;
}

double Physics::GetReactionBeta(const std::string &reaction, const std::string &option)
{
    std::string input = reaction;
    std::regex pattern("([A-Za-z]+)|(\\d+)");
    std::smatch matches;
    std::vector<std::string> characters;
    std::vector<int> numbers;

    while (std::regex_search(input, matches, pattern))
    {
        if (matches[1].matched)
        {
            characters.push_back(matches[1].str());
        }
        else if (matches[2].matched)
        {
            numbers.push_back(std::stoi(matches[2].str()));
        }
        input = matches.suffix().str();
    }

    std::string beam = characters[0] + std::to_string(numbers[0]);
    std::string target = characters[1] + std::to_string(numbers[1]);
    int beamEperA = numbers[2];
    return Physics::GetReactionBeta(beam, target, beamEperA, option);
}

double Physics::GetReactionBeta(const std::string &beam, const std::string &target, const int &beamEperA, const std::string &option)
{
    double m1 = GetNucleiMass(beam, option);
    double m2 = GetNucleiMass(target, option);

    int beamA = std::stoi(beam.substr(2, beam.size() - 2));

    double beam_ke = beamEperA * 1.0 * beamA;
    double beam_energy_tot = beam_ke + m1;
    double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * m1);
    double gamma = beam_energy_tot / m1;

    return mom_beam / (gamma * m1 + m2);
}

double Physics::GetBeamRapidity(const std::string &reaction, const std::string &option)
{
    std::string input = reaction;
    std::regex pattern("([A-Za-z]+)|(\\d+)");
    std::smatch matches;
    std::vector<std::string> characters;
    std::vector<int> numbers;

    while (std::regex_search(input, matches, pattern))
    {
        if (matches[1].matched)
        {
            characters.push_back(matches[1].str());
        }
        else if (matches[2].matched)
        {
            numbers.push_back(std::stoi(matches[2].str()));
        }
        input = matches.suffix().str();
    }

    std::string beam = characters[0] + std::to_string(numbers[0]);
    std::string target = characters[1] + std::to_string(numbers[1]);
    int beamEperA = numbers[2];

    return Physics::GetBeamRapidity(beam, beamEperA, option);
}

double Physics::GetBeamRapidity(const std::string &beam, const int &beamEperA, const std::string &option)
{
    double m1 = GetNucleiMass(beam, option);
    int beamA = std::stoi(beam.substr(2, beam.size() - 2));

    double beam_ke = beamEperA * 1.0 * beamA;
    double beam_energy_tot = beam_ke + m1;
    double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * m1);

    return 0.5 * TMath::Log((beam_energy_tot + mom_beam) / (beam_energy_tot - mom_beam));
}