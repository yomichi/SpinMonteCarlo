#ifndef SW_H
#define SW_H

#include <alps/parapack/worker.h>

#include "../../util/logging.hpp"

namespace potts{

class Worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> SuperClass;

  util::Logger logger_;

  /*
   * parameters
   */
  int q_;
  double T_, beta_;
  double J_;
  double p_;
  unsigned int nsites_, nbonds_;

  double negspin_;
  double m2_coeff_, m4_coeff_;
  double ene_coeff_, ene_const_;
  double V2_;

  /*
   * configuration (need checkpointing)
   */
  alps::mc_steps mcs_;
  std::vector<int> spins_;

  void init(alps::Parameters const& params);
  int update();

public:
  Worker(alps::Parameters const& params);

  ~Worker(){}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs);

  bool is_thermalized() const {return mcs_.is_thermalized();}
  double progress() const {return mcs_.progress();}

  void load(alps::IDump& dp)
  {
    dp >> mcs_
       >> spins_
       ;
  }
  void save(alps::ODump& dp) const
  {
    dp << mcs_
       << spins_
       ;
  }
  void run(alps::ObservableSet& obs);
};

class Evaluator : public alps::parapack::simple_evaluator{
  public:
  Evaluator(alps::Parameters const& params)
    :params_(params)
  {
  }
  ~Evaluator(){}
  void evaluate(alps::ObservableSet& obs) const;
private:
  alps::Parameters params_;
};

} // end of namespace potts

#endif
