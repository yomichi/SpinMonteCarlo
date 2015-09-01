#ifndef SQUARE_LATTICE_HPP
#define SQUARE_LATTICE_HPP

#include <vector>
#include <cmath>

namespace wanglandau{
namespace detail{
inline int modulo(int x, int y) { return x - std::floor(1.0*x/y)*y; }
} // end of namespace detail

class SquareLattice{
public:
  SquareLattice(int L);
  int num_sites() const{ return N_;}
  int num_bonds() const{ return 2*N_;}
  int index(int x, int y) const { return detail::modulo(x,L_) + detail::modulo(y,L_)*L_;}
  int x(int index) const { return index % L_;}
  int y(int index) const { return index / L_;}
  std::vector<int> neighbors(int index) const{ return neighbors_[index];}
  std::vector<int> neighbors(int x, int y) const{ return neighbors(index(x,y));}
private:
  int L_;
  int N_;
  std::vector<std::vector<int> > neighbors_;
};

SquareLattice::SquareLattice(int L):L_(L), N_(L*L)
{
  neighbors_.reserve(N_);
  for(int y=0; y<L_; ++y){
    for(int x=0; x<L_; ++x){
      std::vector<int> neighbor(4);
      neighbor[0] = index(x, y-1);
      neighbor[1] = index(x-1, y);
      neighbor[2] = index(x+1, y);
      neighbor[3] = index(x, y+1);
      neighbors_.push_back(neighbor);
    }
  }
}

} // end of namespace wanglandau

#endif 
