#ifndef POTTS_H
#define POTTS_H

#include <vector>
#include <cmath>
#include "model.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace potts{

template <class Lattice, class RNG01>
class Potts : public Model<Lattice, RNG01>{
public:
  Potts(int q, double h, int L):
    Model<Lattice, RNG01>(L),
    q_(q),
    nsites_(this->lat_.num_sites()),nbonds_(this->lat_.num_bonds()),
    NoverQ_(1.0*nsites_/q_),
    h_(h),
    spins_(nsites_,0),
    nzero_(nsites_),
    ene_(-nbonds_){}
  
  double mag() const{ return nzero_ - NoverQ_;}
  double ene() const{ return ene_ - h_*mag();}

  void update(double beta, RNG01 &rng);

private:
  int q_;
  int nsites_;
  int nbonds_;
  double NoverQ_;

  double h_;
  std::vector<int> spins_;
  int nzero_; // number of s=0 spins
  double ene_;
};

template <class Lattice, class RNG01>
void Potts<Lattice, RNG01>::update(double beta, RNG01 &rng)
{
  for(int site = 0; site < nsites_; ++site){
    const int now  = spins_[site];
    const int cand = q_ * rng();
    int de = 0;
    foreach(int n, this->lat_.neighbors(site)){
      de += now  == spins_[n] ? 1 : 0;
      de -= cand == spins_[n] ? 1 : 0;
    }
    const double p = std::exp(-beta*de);
    if(rng() < p){
      spins_[site] = cand;
      ene_ += de;
      nzero_ -= now == 0 ? 1:0;
      nzero_ += cand == 0 ? 1:0;
    }
  }
}

} // end of namespace potts

#endif
