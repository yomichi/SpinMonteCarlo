#ifndef ISING_H
#define ISING_H

#include <vector>
#include <cmath>
#include "model.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace potts{

template <class Lattice, class RNG01>
class Ising : public Model<Lattice, RNG01>{
public:
  Ising(double h, int L):
    Model<Lattice, RNG01>(L),
    nsites_(this->lat_.num_sites()),nbonds_(this->lat_.num_bonds()),
    h_(h),
    spins_(nsites_,1),
    nup_(nsites_),
    ene_(-nbonds_-h_*nsites_){}
  
  double mag() const{ return 2*nup_ - nsites_; }
  double ene() const{ return ene_;}

  void update(double beta, RNG01 &rng);

private:
  int nsites_;
  int nbonds_;
  double h_;

  std::vector<int> spins_;
  int nup_; // number of up (s=1) spins
  double ene_;
};

template <class Lattice, class RNG01>
void Ising<Lattice, RNG01>::update(double beta, RNG01 &rng)
{
  for(int site = 0; site < nsites_; ++site){
    const int now  = spins_[site];
    double de = 0;
    foreach(int n, this->lat_.neighbors(site)){
      de += now  == spins_[n] ? 2 : -2;
    }
    de += h_*now;
    const double p = std::exp(-beta*de);
    if(rng() < p){
      spins_[site] *= -1;
      ene_ += de;
      nup_ -= now;
    }
  }
}

} // end of namespace potts

#endif
