#ifndef POTTS_H
#define POTTS_H

#include "../../util/squarelattice.hpp"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace potts{

template <class Lattice = util::SquareLattice>
class Potts{
public:
  Potts(int q, int L):
    q_(q),lat_(L),
    nsites_(lat_.num_sites()),nbonds_(lat_.num_bonds()),
    NoverQ_(1.0*nsites_/q_),
    spins_(nsites_,0),
    nzero_(nsites_){}
  
  double mag() const{ return nzero_ - NoverQ_;}
  double ene() const{ return ene_;}

  template <class RNG01>
  void update(double beta, RNG01 &rng);

private:
  int q_;
  Lattice lat_;
  int nsites_;
  int nbonds_;
  double NoverQ_;

  std::vector<int> spins_;
  int nzero_; // number of s=0 spins
  double ene_;
};

template <class Lattice>
template <class RNG01>
void Potts<Lattice>::update(double beta, RNG01 &rng)
{
  for(int site = 0; site < nsites_; ++site){
    const int now  = spins_[site];
    const int cand = q_ * rng();
    int de = 0;
    foreach(int n, lat_.neighbors(site)){
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
