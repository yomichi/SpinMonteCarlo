#include "potts.h"
#include <boost/random.hpp>
#include <ctime>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/iterator/counting_iterator.hpp>

namespace {
template <class T>
inline int binary_search(std::vector<T> xs, T x){
  int left = 0;
  int right = xs.size()-1;
  while(right - left > 1){
    int mid = (right+left)/2;
    if(x == xs[mid]){
      return mid;
    }else if(x < xs[mid]){
      right = mid;
    }else{
      left = mid;
    }
  }
  return left;
}
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  using namespace boost;
  using namespace boost::program_options;
  options_description opt("options");
  opt.add_options()
    ("help,h", "show this message")
    ("q,q", value<int>()->default_value(2), "number of state of a spin")
    ("L,L", value<int>()->default_value(10), "length of lattice")
    ("beta-min", value<double>()->default_value(0.1), "minimum beta")
    ("beta-max", value<double>()->default_value(1.0), "maximum beta")
    ("beta-num", value<int>()->default_value(10), "number of beta")
    ("thermalization", value<int>()->default_value(8), "number of exchange to thermalize")
    ("mcs", value<int>()->default_value(64), "number of exchange to measure")
    ("interval,i", value<int>()->default_value(1024), "temperature exchange interval")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);

  if( vm.count("help") ){
    std::cout << opt << std::endl;
    return 0;
  }

  const int q = vm["q"].as<int>();
  const int L = vm["L"].as<int>();
  const double bmin = vm["beta-min"].as<double>();
  const double bmax = vm["beta-max"].as<double>();
  const int nbeta = vm["beta-num"].as<int>();
  const int therm = vm["thermalization"].as<int>();
  const int MCS = vm["mcs"].as<int>();
  const int interval = vm["interval"].as<int>();

  std::vector<double> betas(nbeta);
  for(int i=0; i<nbeta; ++i){
    betas[i] = bmin + i*(bmax-bmin)/(nbeta-1);
  }
  std::vector<int> beta_index(boost::counting_iterator<int>(0),
                              boost::counting_iterator<int>(nbeta));

  std::vector<int> offsets(nprocs+1,0);
  for(int i=0; i<nbeta%nprocs; ++i){
    offsets[i+1] = offsets[i] + nbeta/nprocs+1;
  }
  for(int i = nbeta%nprocs; i<nprocs; ++i){
    offsets[i+1] = offsets[i] + nbeta/nprocs;
  }
  assert(offsets[nprocs] == nbeta);
  const int offset = offsets[rank];

  const int nbeta_local = offsets[rank+1]-offsets[rank];

  std::vector<potts::Potts<potts::SquareLattice> > potts(nbeta_local, potts::Potts<potts::SquareLattice>(q, L));

  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    rnd(boost::mt19937(static_cast<uint32_t>(std::time(0)+rank*42)), boost::uniform_real<>(0.0, 1.0));

  for(int mcs = 0; mcs < therm+mcs; ++mcs){
    for(int lm=0; lm < interval; ++lm){
      for(int i=0; i<nbeta_local; ++i){
        potts[i].update(betas[offset+i], rnd);
      }
    }
    std::vector<double> enes_local(nbeta, 0.0);
    for(int i=0; i<nbeta_local; ++i){
      enes_local[offset+i] = potts[i].ene();
    }
    std::vector<double> enes;
    if(rank==0) enes.resize(nbeta);
    MPI_Reduce(&enes_local[0], &enes[0], nbeta, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
      for(int i=nbeta-1; i>0; --i){ // from low-temperature to high-T
        const int ihigh = beta_index[i];
        const int ilow = beta_index[i-1];
        const double p = std::exp( (betas[ihigh] - betas[ilow])*(enes[ihigh]-enes[ilow]));
        if(rnd() < p){
          std::swap(betas[ihigh], betas[ilow]);
          std::swap(beta_index[i], beta_index[i-1]);
        }
      }
    }
    MPI_Bcast(&betas[0], nbeta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  return 0;
}
