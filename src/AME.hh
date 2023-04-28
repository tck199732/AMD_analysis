#ifndef AME_hh
#define AME_hh

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

class AME
{
public:
    AME(const std::string &filename = "");
    ~AME() { ; }

    void ReadAMETable(const std::string &filename = "");

    double GetMass(const std::string &symbol);
    double GetMass(const int &Z, const int &A);
    double _GetMassUnphysical(const int &Z, const int &A);

private:
    std::map<std::pair<int, int>, double> MassTable;
    std::map<std::string, std::pair<int, int>> ZATable;

protected:
    double NucleonMass = 938.272;
    fs::path DefaultPath = "database/ame/ame_mass.txt";
};

#endif