#ifndef particle_h
#define particle_h

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TMath.h"

struct particle
{
    int nid, zid;
    double px, py, pz, ekincms;
    double plab, pcms;
    double pt, pzlab, ekinlab, rapidity_lab, rapidity_cms, rapidity_lab_normed;
    double pseudo_rapidity_lab, pseudo_rapidity_cms;
    double thetalab, thetacms, phi;
    double etrans;
    int aid, pid;
    std::string name;
    std::string coordinate = "hira";
    void set_uball_coordinate()
    {
        this->coordinate = "uball";
    }
    void autofill(const double &betacms, const double &beam_rapidity = 1.);

    // for emission time analysis
    // are in cms frame presumbly...
    double x, y, z, t;
    double r;
    void set_xyzt(const double &x, const double &y, const double &z, const double &t);
};

#endif
