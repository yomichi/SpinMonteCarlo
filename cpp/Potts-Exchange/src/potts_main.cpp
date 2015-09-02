#include "potts.h"
#include "ising.h"
#include "../../util/observable.hpp"
#include "../../util/squarelattice.hpp"
#include <iostream>
#include <boost/random.hpp>
#include <ctime>
#include <mpi.h>
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

int main(int argc, char **argv)
{
  typedef util::SquareLattice Lattice;
  typedef boost::variate_generator<boost::mt19937, boost::uniform_real<> > RNG01;
  typedef boost::shared_ptr<potts::Model<Lattice, RNG01> > ModelPtr;

  MPI_Init(&argc, &argv);

  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RNG01 rnd(boost::mt19937(static_cast<uint32_t>(std::time(0)+rank*42)), boost::uniform_real<>(0.0, 1.0));

  using namespace boost;
  using namespace boost::program_options;
  options_description opt("options");
  opt.add_options()
    ("help", "show this message")
    ("Ising", "simulate Ising model")
    ("q,q", value<int>()->default_value(2), "number of state of a spin")
    ("h,h", value<double>()->default_value(0.0), "magnetic field")
    ("L,L", value<int>()->default_value(10), "length of lattice")
    ("beta-min", value<double>()->default_value(0.1), "minimum beta")
    ("beta-max", value<double>()->default_value(1.0), "maximum beta")
    ("beta-num", value<int>()->default_value(10), "number of beta")
    ("thermalization", value<int>()->default_value(8), "number of exchange to thermalize")
    ("mcs", value<int>()->default_value(64), "number of exchange to measure")
    ("interval,i", value<int>()->default_value(1024), "temperature exchange interval")
    ("no-exchange", "switch off temperature exchange")
    ;

  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);

  if( vm.count("help") ){
    std::cout << opt << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
  }

  const bool bIsing = vm.count("Ising");
  const bool noEx = vm.count("no-exchange");

  const int q = vm["q"].as<int>();
  const int L = vm["L"].as<int>();
  const double h = vm["h"].as<double>();
  const double bmin = vm["beta-min"].as<double>();
  const double bmax = vm["beta-max"].as<double>();
  const int nbeta = vm["beta-num"].as<int>();
  const int therm = vm["thermalization"].as<int>();
  const int MCS = vm["mcs"].as<int>();
  const int interval = vm["interval"].as<int>();

  const double dbeta = (bmax-bmin)/(nbeta-1);

  // the i-th worker runs under beta = betas[i]
  // the beta_index[i]-th worker runs under the i-th lowest beta (bmin + i*dbeta)
  // in other words, the value of the i-th lowest beta is betas[beta_index[i]]
  std::vector<double> betas(nbeta);
  std::vector<int> beta_index(nbeta);
  // inv_beta_index is the inverse of beta_index
  std::vector<int> inv_beta_index(nbeta);
  for(int i=0; i<nbeta; ++i){
    betas[i] = bmin + i*dbeta;
    beta_index[i] = i;
    inv_beta_index[i] = i;
  }

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

  std::vector<ModelPtr> models;
  if(bIsing){
    for(int i=0; i<nbeta_local; ++i){
      models.push_back(ModelPtr(new potts::Ising<Lattice,RNG01>(h,L)));
    }
  }else{
    for(int i=0; i<nbeta_local; ++i){
      models.push_back(ModelPtr(new potts::Potts<Lattice,RNG01>(q,h,L)));
    }
  }

  std::vector<util::Observable> obs_ene(nbeta);
  std::vector<util::Observable> obs_mag(nbeta);

  for(int mcs = 0; mcs < therm+MCS; ++mcs){
    for(int lm=0; lm < interval; ++lm){
      for(int i=0; i<nbeta_local; ++i){
        models[i]->update(betas[offset+i], rnd);
      }
    }
    std::vector<double> enes_local(nbeta);
    std::vector<double> mags_local(nbeta);
    for(int i=0; i<nbeta_local; ++i){
      enes_local[offset+i] = models[i]->ene();
      mags_local[offset+i] = models[i]->mag();
    }
    std::vector<double> enes(rank==0?nbeta:0);
    std::vector<double> mags(rank==0?nbeta:0);
    MPI_Reduce(&enes_local[0], &enes[0], nbeta, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mags_local[0], &mags[0], nbeta, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(noEx){
      if(rank == 0){
        for(int i=0; i<nbeta; ++i){
          obs_ene[i] << enes[beta_index[i]];
          obs_mag[i] << mags[beta_index[i]];
        }
        std::clog << "mcs : " << mcs+1 << "/" << MCS+therm << " done." << std::endl;
      }
    }else{
      if(rank == 0){
        for(int i=0; i<nbeta; ++i){
          obs_ene[i] << enes[beta_index[i]];
          obs_mag[i] << mags[beta_index[i]];
        }

        for(int i=nbeta-1; i>0; --i){ // from low-temperature to high-T
          const int ihigh = beta_index[i];
          const int ilow = beta_index[i-1];
          const double p = std::exp( (betas[ihigh] - betas[ilow])*(enes[ihigh]-enes[ilow]));
          if(rnd() < p){
            std::swap(betas[ihigh], betas[ilow]);
            std::swap(beta_index[i], beta_index[i-1]);
            std::swap(inv_beta_index[beta_index[i]], inv_beta_index[beta_index[i-1]]);
          }
        }
      }
      MPI_Bcast(&betas[0], nbeta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(rank==0){
        if(mcs == therm){
          for(int i=0; i<nbeta; ++i){
            obs_ene[i].reset();
            obs_mag[i].reset();
          }
        }
        std::clog << "mcs : " << mcs+1 << "/" << MCS+therm << " done." << std::endl;
      }
    }
  }

  if(rank==0){
    for(int i=0; i<nbeta; ++i){
      std::cout << betas[beta_index[i]]
         << " " << obs_ene[i].mean() << " " << obs_ene[i].error()
         << " " << obs_mag[i].mean() << " " << obs_mag[i].error()
         << std::endl;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
