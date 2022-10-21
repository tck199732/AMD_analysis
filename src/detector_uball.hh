
#ifndef detector_uball_h
#define detector_uball_h 1

#include <iostream>
#include <fstream>
#include <array>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2D.h"

class detector_uball
{

protected:
    static const int NumRing = 10; // use ring 2, 3, 4, 5, 7, 8
    static const int NumDet = 15;
    static const int MaxA = 100;
    static const int MaxZ = 100;

public:
    detector_uball();
    ~detector_uball();

    void turn_off_mhit() { this->opt_csi_mhit = 0; }
    void turn_off_ekincut() { this->opt_ekincut = 0; }
    void turn_off_anglecut() { this->opt_anglecut = 0; }

    void read_geometry(const std::string &filename);
    void read_threshold(const int &ring, const std::string &filename);
    void remove_csi(const int &ring, const int &det);

    double get_thetamin(const int &ring, const int &det) { return this->theta_range[ring][det][0]; }
    double get_thetamax(const int &ring, const int &det) { return this->theta_range[ring][det][1]; }
    double get_phimin(const int &ring, const int &det) { return this->phi_range[ring][det][0]; }
    double get_phimax(const int &ring, const int &det) { return this->phi_range[ring][det][1]; }

    int get_ring_id(const double &thetalab);
    int get_det_id(const double &thetalab, const double &phi);
    double get_thres(const double &thetalab, const int &aid, const int &zid);
    double get_thres(const int &ring_id, const int &aid, const int &zid) { return this->ekinlabcut[ring_id][aid][zid][0]; }

    void reset_csi();
    int get_csi_hits(const int &ring, const int &det) { return this->csi_hits[ring][det]; }
    void add_csi_hit(const int &ring, const int &det) { this->csi_hits[ring][det]++; }
    void add_csi_hit(const double &thetalab, const double &phi);
    void hira_coordinate();
    void reset_phi_range();

    bool cover(const double &thetalab, const double &phi);
    bool punch_thr(const double &ekinlab, const double &thetalab, const int &aid, const int &zid);
    bool ready_csi(const double &thetalab, const double &phi);

    TH2D *get_hit_map();

private:
    int csi_hits[NumRing][NumDet];
    std::array<double, 2> theta_range[NumRing][NumDet];
    std::array<double, 2> phi_range[NumRing][NumDet];
    std::array<double, 2> ekinlabcut[NumRing][MaxA][MaxZ];

    bool opt_csi_mhit;
    bool opt_ekincut;
    bool opt_anglecut;
};

#endif
