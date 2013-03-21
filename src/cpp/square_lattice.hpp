#ifndef SQUARE_LATTICE_HPP
#define SQUARE_LATTICE_HPP

#include "lattice.hpp"

class square_lattice : public lattice{
public:
  square_lattice(int L) :lattice(L*L,2*L*L, num_neighbor_sites()) {init_impl(L,L);}
  square_lattice(int L, int W) :lattice(L*W,2*L*W, num_neighbor_sites()) {init_impl(L,W);}
  square_lattice(std::vector<int> const &L) :lattice(L[0]*L[1], 2*L[0]*L[1],num_neighbor_sites()) {init_impl(L[0],L[1]);}
  virtual void init(std::vector<int> const &L) { init(L[0],L[1]); }
  void init(int L){ init(L,L);}
  void init(int L, int W)
  {
    reset(L*W,2*L*W,num_neighbor_sites());
    init_impl(L,W);
  }
  virtual ~square_lattice(){}

  virtual static int num_neighbor_sites() const{ return 4;}

private:
  void init_impl(int L, int W)
  {
    L_ = L;
    W_ = W;
    for(int s=0; s<num_sites(); ++s){
      const int x = site2x(s);
      const int y = site2y(s);
      neighbor_sites_[s][0] = xy2site(x-1,y);
      neighbor_sites_[s][1] = xy2site(x+1,y);
      neighbor_sites_[s][2] = xy2site(x,y-1);
      neighbor_sites_[s][3] = xy2site(x,y+1);
    }
    for(int b=0; b<num_bonds(); ){
      const int s = b/2;

    }
  }
  int xy2site(int x, int y) const
  {
    while(x < 0) x+=L_;
    x %= L_;
    while(y < 0) y+=W_;
    y %= W_;
    return y*L_+x
  }
  int site2x(int site) const{ return site%L_;}
  int site2y(int site) const{ return site/L_;}

  int L_;
  int W_;
};

#endif
