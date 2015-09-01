#include "squarelattice.hpp"
#include "wanglandau.h"
#include <boost/random.hpp>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

int main(int argc, char **argv)
{
  using namespace boost;
  using namespace boost::program_options;
  options_description opt("options");
  opt.add_options()
    ("help,h", "show this message")
    ("q,q", value<int>()->default_value(2), "number of state of a spin")
    ("L,L", value<int>()->default_value(10), "length of lattice")
    ("flatness,f", value<double>()->default_value(0.8), "target flatness of histogram")
    ("final,F", value<double>()->default_value(1.0e-8), "minimum updating factor")
    ("interval,i", value<int>()->default_value(32), "interval between checks for flatness")
    ("verbose,v",  "show histogram verbose info");

  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);

  if( vm.count("help") ){
    std::cout << opt << std::endl;
    return 0;
  }

  const int q = vm["q"].as<int>();
  const int L = vm["L"].as<int>();
  const double flatness = vm["flatness"].as<double>();
  const double final_factor = vm["final"].as<double>();
  const int check_interval = vm["interval"].as<int>();
  const bool verbose = vm.count("verbose");

  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    rnd(boost::mt19937(static_cast<uint32_t>(std::time(0))), boost::uniform_real<>(0.0, 1.0));

  wanglandau::SquareLattice lattice(L);

  const int nsites = lattice.num_sites();
  const int nbonds = lattice.num_bonds();

  wanglandau::WangLandau wl(nbonds+1, flatness, final_factor, (q==2?0.4:0.8)*(nbonds+1),
                            std::log(q), nbonds);

  std::vector<int> spins(nsites, 0);
  int npara = nbonds;

  while( !wl.finished()){
    for(int i=0; i<check_interval; ++i){
      for(int s=0; s<nsites; ++s){
        const int site = nsites*rnd();
        const int pspin = spins[site]; // present
        const int cspin = q*rnd();     // candidate
        int cpara = npara;
        foreach(int neighbor, lattice.neighbors(site)){
          cpara -= spins[neighbor] == pspin ? 1 : 0;
          cpara += spins[neighbor] == cspin ? 1 : 0;
        }
        if(rnd() < wl.prob(npara, cpara)){
          spins[site] = cspin;
          npara = cpara;
        }
        wl.visit(npara);
      }
    }
    wl.update(verbose);
  }

  return 0;
}
