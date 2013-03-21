#ifndef CHAIN_LATTICE_HPP
#define CHAIN_LATTICE_HPP

#include "lattice.hpp"

class chain_lattice : public lattice{
public:
  chain_lattice(int L) :lattice(L,L) {init_impl(L);}
  chain_lattice(std::vector<int> const &L) :lattice(L[0], L[0]) {init_impl(L[0]);}
  virtual void init(std::vector<int> const &L) { init(L[0]); }
  void init(int L)
  {
    reset(L,L,num_neighbor_sites());
    init_impl(L);
  }
  virtual ~chain_lattice(){}

  virtual static int num_neighbor_sites() const{ return 2;}

private:
  void init_impl(int L)
  {
    for(int s=0; s<num_sites(); ++s){
      source_[s] = s;
      target_[s] = (s+1)%L;
      neighbor_sites_[s][0] = (s-1+L)%L;
      neighbor_sites_[s][1] = (s+1)%L;
    }
  }
};

#endif
