#include "Microball.hh"

Microball::Microball()
{
    OptionChargedParticle = 1;
    OptionMultipleHit = 1;
    OptionKinergy = 1;
    OptinoCoverage = 1;

    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            this->Theta[ring][det][0] = -1;
            this->Phi[ring][det][0] = -1;
            this->Theta[ring][det][1] = -1;
            this->Phi[ring][det][1] = -1;
            this->CsIHits[ring][det] = 0;
        }

        for (int aid = 0; aid < this->MaxA; aid++)
        {
            for (int zid = 0; zid < this->MaxA; zid++)
            {
                this->KinergyThreshold[ring][aid][zid] = 0.;
            }
        }
    }
}

void Microball::ReadGeometry(const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }
    std::ifstream infile(filename.c_str());
    infile.ignore(99, '\n');

    int ring, det;
    double theta1, theta2, phi1, phi2;
    while (infile >> ring)
    {
        infile >> det >> theta1 >> theta2 >> phi1 >> phi2;

        this->Theta[ring][det][0] = theta1;
        this->Theta[ring][det][1] = theta2;
        this->Phi[ring][det][0] = phi1;
        this->Phi[ring][det][1] = phi2;
    }
    infile.close();
}

void Microball::ReadConfiguration(const std::string &reaction, const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }

    std::ifstream stream(filename.c_str());
    stream.ignore(99, '\n');
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream iss(line);
        std::string sys;
        int ring_id, det_id;
        iss >> sys >> ring_id;
        if (sys != reaction)
        {
            continue;
        }

        std::vector<int> det_array;
        while (iss >> det_id)
        {
            det_array.push_back(det_id);
        }

        for (int det = 0; det < this->NumDet; det++)
        {
            if (std::find(det_array.begin(), det_array.end(), det) - det_array.begin() == det_array.size())
            {
                this->RemoveCsI(ring_id, det);
            }
        }
    }
}

void Microball::RemoveCsI(const int &ring, const int &det)
{
    this->Theta[ring][det][0] = -1;
    this->Phi[ring][det][0] = -1;
    this->Theta[ring][det][1] = -1;
    this->Phi[ring][det][1] = -1;
    this->CsIHits[ring][det] = 0;
}

void Microball::ResetCsI()
{
    for (int i = 0; i < this->NumRing; i++)
    {
        for (int j = 0; j < this->NumDet; j++)
        {
            this->CsIHits[i][j] = 0;
        }
    }
}

/**
 * @brief Read threshold kinetic energy (data is in unit of MeV not MeV/A)
 *
 * @param ring
 * @param filename
 */
void Microball::ReadThreshold(const int &ring, const std::string &filename)
{
    if (!fs::exists(filename))
    {
        std::string msg = Form("file does not exist : %s", filename.c_str());
        throw std::invalid_argument(msg.c_str());
    }

    std::ifstream infile(filename.c_str());
    infile.ignore(99, '\n');

    int A, Z;
    double thres;
    while (infile >> A)
    {
        infile >> Z >> thres;
        if (A > this->maxA || Z > this->maxZ)
        {
            continue;
        }
        this->KinergyThreshold[ring][A][Z] = thres;
    }

    TF1 *fcn = new TF1("fcn", "[0]*x*x + [1]*x + [2]");

    std::vector<std::vector<double>> paras(this->MaxZ, std::vector<double>(3));

    for (int zid = 0; zid < this->MaxZ; zid++)
    {
        TGraph *gr = new TGraph();
        int npoints = 0;
        for (int aid = 0; aid < this->MaxA; aid++)
        {
            double thres = this->KinergyThreshold[ring][aid][zid];
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
                this->KinergyThreshold[ring][aid][zid] = 9999.;
            }
            else if (this->KinergyThreshold[ring][aid][zid] == 0.)
            {
                double num_a = 1.0 * aid;
                double fitted_thres = fcn->EvalPar(&num_a, &paras[zid][0]);
                this->KinergyThreshold[ring][aid][zid] = (fitted_thres > 0) ? fitted_thres : 0.;
            }
        }
    }

    infile.close();
    return;
}

void Microball::HiRA_Coordinate()
{
    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            if (this->Phi[ring][det][0] != -1. && this->Phi[ring][det][1] != -1. && this->Theta[ring][det][0] != -1. && this->Theta[ring][det][1] != -1.)
            {
                this->Phi[ring][det][0] += 90.;
                this->Phi[ring][det][1] += 90.;

                if (this->Phi[ring][det][0] >= 360.)
                {
                    this->Phi[ring][det][0] -= 360.;
                }
                if (this->Phi[ring][det][1] > 360.)
                {
                    this->Phi[ring][det][1] -= 360.;
                }
                if (this->Phi[ring][det][0] < 0.)
                {
                    this->Phi[ring][det][0] += 360.;
                }
                if (this->Phi[ring][det][1] <= 0.)
                {
                    this->Phi[ring][det][1] += 360.;
                }
            }
        }
    }
}

void Microball::ResetPhiRange()
{
    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {

            if (this->Phi[ring][det][0] != -1. && this->Phi[ring][det][1] != -1. && this->Theta[ring][det][0] != -1. && this->Theta[ring][det][1] != -1.)
            {
                if (this->Phi[ring][det][0] >= 360.)
                {
                    this->Phi[ring][det][0] -= 360.;
                }
                else if (this->Phi[ring][det][1] > 360.)
                {
                    this->Phi[ring][det][1] -= 360.;
                }
                else if (this->Phi[ring][det][0] < 0.)
                {
                    this->Phi[ring][det][0] += 360.;
                }
                else if (this->Phi[ring][det][1] <= 0.)
                {
                    this->Phi[ring][det][1] += 360.;
                }
            }
        }
    }
}

int Microball::GetRingID(const double &thetalab)
{
    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            if (thetalab >= this->Theta[ring][det][0] && thetalab <= this->Theta[ring][det][1])
            {
                return ring;
            }
        }
    }
    return -1;
}

int Microball::GetDetID(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    if (ring_id == -1)
    {
        return -1;
    }

    for (int det = 0; det < this->NumDet; det++)
    {

        if (this->Phi[ring_id][det][0] < this->Phi[ring_id][det][1])
        {
            if (phi >= this->Phi[ring_id][det][0] && phi <= this->Phi[ring_id][det][1])
            {
                return det;
            }
        }
        else if (this->Phi[ring_id][det][0] > this->Phi[ring_id][det][1])
        {
            if ((phi >= this->Phi[ring_id][det][0] && phi <= 360) || (phi >= 0.0 && phi <= this->Phi[ring_id][det][1]))
            {
                return det;
            }
        }
    }
    return -1;
}

double Microball::GetThreshold(const double &thetalab, const int &aid, const int &zid)
{
    int ring_id = this->GetRingID(thetalab);
    return this->GetThreshold(ring_id, aid, zid);
}

bool Microball::IsChargedParticle(const int &Z)
{
    return (Z > 0 || !OptionChargedParticle) ? 1 : 0;
}

bool Microball::IsCovered(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    int det_id = this->GetDetID(thetalab, phi);
    return (ring_id != -1 && det_id != -1) || !this->OptinoCoverage;
}

bool Microball::IsAccepted(const double &ekinlab, const double &thetalab, const int &aid, const int &zid)
{
    if (!this->OptionKinergy)
    {
        return true;
    }
    int ring_id = this->GetRingID(thetalab);
    return (ring_id != -1 && aid > this->MaxA && zid > this->MaxZ) ? ekinlab >= this->KinergyThreshold[ring_id][aid][zid] / (double)aid : false;
}

bool Microball::IsReadyCsI(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    int det_id = this->GetDetID(thetalab, phi);
    return (this->CsIHits[ring_id][det_id] == 0 || !this->OptionMultipleHit) ? 1 : 0;
}

void Microball::AddCsIHit(const double &thetalab, const double &phi)
{
    int ring_id = this->GetRingID(thetalab);
    int det_id = this->GetDetID(thetalab, phi);

    if (ring_id != -1 && det_id != -1)
    {
        this->CsIHits[ring_id][det_id]++;
    }
    if (this->CsIHits[ring_id][det_id] >= 2 && this->OptionMultipleHit == 1)
    {
        std::cout << Form("warning : ring %i, det %i count > 1", ring_id, det_id) << std::endl;
    }
    return;
}

int Microball::GetCsIHits()
{
    int counter = 0;
    for (int ring = 0; ring < this->NumRing; ring++)
    {
        for (int det = 0; det < this->NumDet; det++)
        {
            counter += this->GetCsIHits(ring, det);
        }
    }
}