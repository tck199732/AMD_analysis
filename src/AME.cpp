#include "AME.hh"

AME::AME(const std::string &filename)
{
    this->ReadAMETable(filename);
}

void AME::ReadAMETable(const std::string &filename)
{

    std::string path = filename;
    fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
    fs::path default_path = PROJECT_DIR / DefaultPath;

    if (path.empty())
    {
        path = default_path;
    }

    MassTable.clear();
    ZATable.clear();

    std::ifstream infile(path.c_str());
    infile.ignore(99, '\n');

    std::string symbol;
    int Z, A;
    double mass;                       // MeV/c^2
    double binding_enregy_per_nucleon; // MeV
    while (infile >> Z >> A >> mass >> binding_enregy_per_nucleon)
    {
        std::pair pair = std::make_pair(Z, A);
        MassTable[pair] = mass;
        ZATable[symbol] = pair;
    }
}

double AME::GetMass(const int &Z, const int &A)
{
    if (this->MassTable.size() == 0)
    {
        this->ReadAMETable();
    }

    if (this->MassTable.count({Z, A}) == 0)
    {
        return this->_GetMassUnphysical(Z, A);
    }
    return this->MassTable[{Z, A}];
}

double AME::_GetMassUnphysical(const int &Z, const int &A)
{
    if (Z > 0 && A >= Z)
    {
        return A * this->NucleonMass;
    }
    else if (Z > A)
    {
        std::cout << "Z: " << Z << " A: " << A << std::endl;
        throw std::invalid_argument("Z and A are not physical!");
    }
    return -1e10;
}

double AME::GetMass(const std::string &symbol)
{
    if (this->ZATable.size() == 0)
    {
        this->ReadAMETable();
    }

    if (this->ZATable.count(symbol) == 0)
    {
        throw std::invalid_argument("Symbol is not valid!");
    }

    std::pair<int, int> ZA = this->ZATable[symbol];
    return this->GetMass(ZA.first, ZA.second);
}