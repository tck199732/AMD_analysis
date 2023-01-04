#include "particle.hh"

void particle::autofill(const double &betacms, const double &beam_rapidity)
{

    double nuc_mass = 938.272;
    this->aid = this->nid + this->zid;
    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));
    this->pt = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    this->pcms = TMath::Sqrt(pow(this->pt, 2.) + pow(this->pz, 2.));
    this->ekincms = TMath::Sqrt(pow(this->pcms, 2.) + pow(nuc_mass, 2.)) - nuc_mass;
    this->thetacms = TMath::ATan2(this->pt, this->pz) * TMath::RadToDeg();

    this->pzlab = gamma * (this->pz + betacms * (this->ekincms + nuc_mass));
    this->plab = TMath::Sqrt(pow(this->pt, 2.) + pow(this->pzlab, 2.));
    this->thetalab = TMath::ATan2(this->pt, this->pzlab) * TMath::RadToDeg();
    this->ekinlab = TMath::Sqrt(pow(this->plab, 2.) + pow(nuc_mass, 2.)) - nuc_mass;
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg();

    if (this->coordinate == "hira")
    {
        if (this->phi < 0)
        {
            this->phi += 360.;
        }
        if (this->phi > 360)
        {
            this->phi -= 360.;
        }
    }
    else if (this->coordinate == "uball")
    {
        if (this->phi < -18.)
        {
            this->phi += 360.;
        }
        if (this->phi > 342.)
        {
            this->phi -= 360.;
        }
    }
    this->etrans = TMath::Sqrt(pow(this->pt, 2.) + pow(nuc_mass, 2.)) - nuc_mass;

    this->rapidity_lab = 0.5 * TMath::Log((this->ekinlab + this->pzlab + nuc_mass) / (this->ekinlab - this->pzlab + nuc_mass));
    this->rapidity_lab_normed = this->rapidity_lab / beam_rapidity;

    this->rapidity_cms = 0.5 * TMath::Log((this->ekincms + this->pz + nuc_mass) / (this->ekincms - this->pz + nuc_mass));

    this->pseudo_rapidity_lab = 0.5 * TMath::Log((this->plab + this->pzlab) / (this->plab - this->pzlab));

    this->pseudo_rapidity_cms = 0.5 * TMath::Log((this->pcms + this->pz) / (this->pcms - this->pz));

    this->pid = 99;
    this->name = "";
    if (this->nid == 1 && this->zid == 0)
    {
        this->pid = 0;
        this->name = "n";
    }
    else if (this->nid == 0 && this->zid == 1)
    {
        this->pid = 1;
        this->name = "p";
    }
    else if (this->nid == 1 && this->zid == 1)
    {
        this->pid = 2;
        this->name = "d";
    }
    else if (this->nid == 2 && this->zid == 1)
    {
        this->pid = 3;
        this->name = "t";
    }
    else if (this->nid == 1 && this->zid == 2)
    {
        this->pid = 4;
        this->name = "3He";
    }
    else if (this->nid == 2 && this->zid == 2)
    {
        this->pid = 5;
        this->name = "4He";
    }

    if (this->x != 0. && this->y != 0 && this->z != 0. && this->t != 0)
    {
        this->r = TMath::Sqrt(pow(this->x, 2.) + pow(this->y, 2.) + pow(this->z, 2.));
    }
    else
    {
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
        this->t = 0.;
        this->r = 0.;
    }
}

void particle::set_xyzt(const double &x, const double &y, const double &z, const double &t)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    this->r = TMath::Sqrt(pow(this->x, 2.) + pow(this->y, 2.) + pow(this->z, 2.));
}
