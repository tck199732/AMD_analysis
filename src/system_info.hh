#ifndef system_info_h
#define system_info_h

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "TMath.h"
#include "TString.h"

struct coll_system
{
    std::string beam, target;
    int energy;
    std::string name;
    double betacms, rapidity_beam;
    void autofill();
};

class system_info
{
private:
    std::map<std::string, coll_system> systems;

public:
    system_info();
    ~system_info();

    double get_rapidity_beam(const std::string &system_name)
    {
        return this->systems[system_name].rapidity_beam;
    }
    double get_betacms(const std::string &system_name)
    {
        return this->systems[system_name].betacms;
    }
};

#endif
