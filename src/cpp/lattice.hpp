#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <cassert>
#include <vector>

class lattice{
public:
  virtual init(std::vector<int> sizes) =0;
  int num_sites() const { return nsites_;}
  int num_bonds() const { return nbonds_;}
  int source(int bond) const { assert(0<=bond && bond<nbonds_); return source_[bond];}
  int target(int bond) const { assert(0<=bond && bond<nbonds_); return target_[bond];}
  std::vector<int> const& neighbor_sites(int site) const{ assert(0<=site && site<nsites_); return neighbor_sites_[site];}
  virtual static int num_neighbor_sites() const =0;
protected:
  lattice(int nsites, int nbonds, nneighbors)
     :nsites_(nsites)
     ,nbonds_(nbonds)
     ,source_(nbonds)
     ,target_(nbonds)
     ,neighbor_sites_(nsites, std::vector<int>(nneighbors))
  {}
  void reset(int nsites, int nbonds, int nneighbors)
  {
    nsites_ = nsites;
    nbonds_ = nbonds;
    source_.resize(nbonds);
    target_.resize(nbonds);
    neighbor_sites_.resize(nsites, std::vector<int>(nneighbors));
  }
  virtual ~lattice(){}

  int nsites_;
  int nbonds_;
  std::vector<int> source_;
  std::vector<int> target_;
  std::vector<std::vector<int> > neighbor_sites_;
};

#endif
