#include "detector_uball.hh"

detector_uball::detector_uball()
{

    opt_csi_mhit = 1;
    opt_ekincut = 1;
    opt_anglecut = 1;

    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {

            this->theta_range[ring][det][0] = -1;
            this->phi_range[ring][det][0] = -1;
            this->theta_range[ring][det][1] = -1;
            this->phi_range[ring][det][1] = -1;
            this->csi_hits[ring][det] = 0;
        }

        for (int aid = 0; aid < this->MaxA; aid++)
        {
            for (int zid = 0; zid < this->MaxA; zid++)
            {
                this->ekinlabcut[ring][aid][zid][0] = 0.;
                this->ekinlabcut[ring][aid][zid][1] = 1000.;
            }
        }
    }
}

detector_uball::~detector_uball() { ; }

void detector_uball::read_geometry(const std::string &filename)
{

    std::ifstream infile(filename.c_str());

    if (!infile)
    {
        std::cout << "no such file." << std::endl;
        return;
    }

    infile.ignore(99, '\n');

    int ring, det;
    double theta1, theta2, phi1, phi2;
    while (infile >> ring)
    {
        infile >> det >> theta1 >> theta2 >> phi1 >> phi2;

        this->theta_range[ring][det][0] = theta1;
        this->theta_range[ring][det][1] = theta2;
        this->phi_range[ring][det][0] = phi1;
        this->phi_range[ring][det][1] = phi2;
    }

    infile.close();
    std::cout << "reading geometry : phi range = (-18, 342) ..." << std::endl;
}

void detector_uball::remove_csi(const int &ring, const int &det)
{
    this->theta_range[ring][det][0] = -1;
    this->phi_range[ring][det][0] = -1;
    this->theta_range[ring][det][1] = -1;
    this->phi_range[ring][det][1] = -1;
    this->csi_hits[ring][det] = 0;
}

void detector_uball::reset_csi()
{
    for (int i = 0; i < this->NumRing; i++)
    {
        for (int j = 0; j < this->NumDet; j++)
        {
            this->csi_hits[i][j] = 0;
        }
    }
}

void detector_uball::read_threshold(const int &ring, const std::string &filename)
{

    // threshold energy data are in the unit of MeV, not MeV/A

    std::ifstream infile(filename.c_str());
    if (!infile)
    {
        std::cout << "no such file." << std::endl;
        return;
    }

    std::cout << "reading threshold for ring " << ring << "..." << std::endl;
    infile.ignore(99, '\n');

    std::array<int, 2> buffer_int;
    double buffer_double;
    while (infile >> buffer_int[0])
    {
        infile >> buffer_int[1] >> buffer_double;
        this->ekinlabcut[ring][buffer_int[0]][buffer_int[1]][0] = buffer_double;
    }

    TF1 *fcn = new TF1("fcn", "[0]*x*x + [1]*x + [2]");

    std::vector<std::vector<double>> paras(this->MaxZ, std::vector<double>(3));

    for (int zid = 0; zid < this->MaxZ; zid++)
    {
        TGraph *gr = new TGraph();
        int npoints = 0;
        for (int aid = 0; aid < this->MaxA; aid++)
        {
            double thres = this->ekinlabcut[ring][aid][zid][0];
            if (thres > 0.)
            {
                gr->SetPoint(npoints, aid, thres);
                npoints++;
            }
        }
        if (npoints == 0)
        {
            delete gr;
            continue;
        }
        gr->Fit(fcn, "Q");
        paras[zid][0] = fcn->GetParameter(0);
        paras[zid][1] = fcn->GetParameter(1);
        paras[zid][2] = fcn->GetParameter(2);
        delete gr;
    }

    for (int zid = 0; zid < this->MaxZ; zid++)
    {
        for (int aid = 0; aid < this->MaxA; aid++)
        {
            if (zid > aid)
            {
                this->ekinlabcut[ring][aid][zid][0] = 9999.;
            }
            else if (this->ekinlabcut[ring][aid][zid][0] == 0.)
            {

                double num_a = 1.0 * aid;
                double fitted_thres = fcn->EvalPar(&num_a, &paras[zid][0]);

                if (fitted_thres > 0)
                {
                    this->ekinlabcut[ring][aid][zid][0] = fitted_thres;
                }
                else
                {
                    this->ekinlabcut[ring][aid][zid][0] = 0.001;
                }
            }
        }
    }

    infile.close();
    return;
}

void detector_uball::hira_coordinate()
{

    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            if (this->phi_range[ring][det][0] != -1. && this->phi_range[ring][det][1] != -1. && this->theta_range[ring][det][0] != -1. && this->theta_range[ring][det][1] != -1.)
            {
                this->phi_range[ring][det][0] += 90.;
                this->phi_range[ring][det][1] += 90.;

                if (this->phi_range[ring][det][0] >= 360.)
                {
                    this->phi_range[ring][det][0] -= 360.;
                }
                if (this->phi_range[ring][det][1] > 360.)
                {
                    this->phi_range[ring][det][1] -= 360.;
                }
                if (this->phi_range[ring][det][0] < 0.)
                {
                    this->phi_range[ring][det][0] += 360.;
                }
                if (this->phi_range[ring][det][1] <= 0.)
                {
                    this->phi_range[ring][det][1] += 360.;
                }
            }
        }
    }

    std::cout << "converting angles to hira coordinate : phi range = (0, 360)" << std::endl;
}

void detector_uball::reset_phi_range()
{
    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {

            if (this->phi_range[ring][det][0] != -1. && this->phi_range[ring][det][1] != -1. && this->theta_range[ring][det][0] != -1. && this->theta_range[ring][det][1] != -1.)
            {
                if (this->phi_range[ring][det][0] >= 360.)
                {
                    this->phi_range[ring][det][0] -= 360.;
                }
                else if (this->phi_range[ring][det][1] > 360.)
                {
                    this->phi_range[ring][det][1] -= 360.;
                }
                else if (this->phi_range[ring][det][0] < 0.)
                {
                    this->phi_range[ring][det][0] += 360.;
                }
                else if (this->phi_range[ring][det][1] <= 0.)
                {
                    this->phi_range[ring][det][1] += 360.;
                }
            }
        }
    }
}

int detector_uball::get_ring_id(const double &thetalab)
{

    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            if (thetalab >= this->theta_range[ring][det][0] && thetalab <= this->theta_range[ring][det][1])
            {
                return ring;
            }
        }
    }
    return -1;
}

int detector_uball::get_det_id(const double &thetalab, const double &phi)
{

    int ring_id = this->get_ring_id(thetalab);
    if (ring_id == -1)
    {
        return -1;
    }

    for (int det = 0; det < this->NumDet; det++)
    {

        if (this->phi_range[ring_id][det][0] < this->phi_range[ring_id][det][1])
        {
            if (phi >= this->phi_range[ring_id][det][0] && phi <= this->phi_range[ring_id][det][1])
            {
                return det;
            }
        }
        else if (this->phi_range[ring_id][det][0] > this->phi_range[ring_id][det][1])
        {
            if ((phi >= this->phi_range[ring_id][det][0] && phi <= 360) || (phi >= 0.0 && phi <= this->phi_range[ring_id][det][1]))
            {
                return det;
            }
        }
    }
    return -1;
}

double detector_uball::get_thres(const double &thetalab, const int &aid, const int &zid)
{
    int ring_id = this->get_ring_id(thetalab);
    return this->get_thres(ring_id, aid, zid);
}

bool detector_uball::cover(const double &thetalab, const double &phi)
{

    int ring_id = this->get_ring_id(thetalab);
    int det_id = this->get_det_id(thetalab, phi);
    return (ring_id != -1 && det_id != -1);
}

bool detector_uball::punch_thr(const double &ekinlab, const double &thetalab, const int &aid, const int &zid)
{
    if (aid > this->MaxA || zid > this->MaxZ)
    {
        return false;
    }
    int ring_id = this->get_ring_id(thetalab);
    return (ring_id != -1) ? ekinlab >= this->ekinlabcut[ring_id][aid][zid][0] / (double)aid : false;
}

bool detector_uball::ready_csi(const double &thetalab, const double &phi)
{
    int ring_id = this->get_ring_id(thetalab);
    int det_id = this->get_det_id(thetalab, phi);
    return (this->csi_hits[ring_id][det_id] == 0) ? 1 : 0;
}

void detector_uball::add_csi_hit(const double &thetalab, const double &phi)
{

    int ring_id = this->get_ring_id(thetalab);
    int det_id = this->get_det_id(thetalab, phi);

    if (ring_id != -1 && det_id != -1)
    {
        this->csi_hits[ring_id][det_id]++;
    }
    if (this->csi_hits[ring_id][det_id] >= 2 && this->opt_csi_mhit == 1)
    {
        std::cout << Form("warning : ring %i, det %i count > 1", ring_id, det_id) << std::endl;
    }
    return;
}

TH2D *detector_uball::get_hit_map()
{

    TH2D *h2 = new TH2D("uball_hitmap", "uball_hitmap", this->NumDet, -0.5, this->NumDet - 0.5, this->NumRing, -0.5, this->NumRing - 0.5);

    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            if (this->theta_range[ring][det][0] != -1 && this->theta_range[ring][det][1] != -1 && this->phi_range[ring][det][0] != -1 && this->phi_range[ring][det][1] != -1)
            {
                h2->SetBinContent(det + 1, ring + 1, 10. * (ring + 1) + (det + 1));
            }
        }
    }
    return h2;
}
