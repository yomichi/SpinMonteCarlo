#include "SW.hpp"

#include "../../util/union_find.hpp"
#include "../../util/pow.hpp"

#include <cmath>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/timer/timer.hpp>

namespace potts{

Worker::Worker(alps::Parameters const& params):
  logger_(std::string("log_")
          + boost::lexical_cast<std::string>(params["TASK_ID"])
          + "_"
          + boost::lexical_cast<std::string>(params["CLONE_ID"])
          ),
  SuperClass(params),
  mcs_(params),
  T_(evaluate("T", params)),
  beta_(1.0/T_),
  J_(params.value_or_default("J", 1.0)),
  nsites_(num_sites()),
  nbonds_(num_bonds()),
  spins_(nsites_, 1),
  V2_(nsites_*nsites_)
{
  init(params);
}

void Worker::init(alps::Parameters const& params)
{
  if(params.defined("MODEL")){
    std::string modelname = params["MODEL"];
    if(modelname ==  std::string("Ising")){
      ene_const_ = nbonds_ * J_;
      q_ = 2;
      J_ *= 2;
    }else if(modelname == std::string("Potts")){
      q_ = evaluate("q", params);
      ene_const_ = 0.0;
    }else{
      std::cerr << modelname << " model is not implemented.";
      exit(127);
    }
  }else{
    std::cerr << "parameter \"MODEL\" is not defined.";
    exit(127);
  }
  negspin_ = 1.0 / (q_-1);
  p_ = -boost::math::expm1(-beta_*J_);
  m2_coeff_ = negspin_;
  m4_coeff_ = (1.0 + negspin_*negspin_*negspin_)/q_;
  ene_coeff_ = J_ / p_;
}

void Worker::init_observables(alps::Parameters const&, alps::ObservableSet& obs)
{
  obs << alps::RealObservable("Time");
  obs << alps::RealObservable("Speed");
  obs << alps::RealObservable("Number of Sites");
  obs << alps::RealObservable("Magnetization^2");
  obs << alps::RealObservable("Magnetization^4");
  obs << alps::RealObservable("|Magnetization|");
  obs << alps::RealObservable("Energy");
  obs << alps::RealObservable("Susceptibility");
  obs << alps::RealObservable("Activated Bonds");
  obs << alps::RealObservable("Activated Bonds^2");
}


void Worker::run(alps::ObservableSet& obs)
{
  boost::timer::cpu_timer tm;
  SCOPED_LOGGING(logger_, "Worker::run()");

  ++mcs_;

  LOGGING(logger_, "clusterize");

  int activated = 0;
  std::vector<union_find::Node> nodes(nsites_);
  foreach(bond_descriptor b, bonds()){
    const int lsite = source(b);
    const int rsite = target(b);
    if( spins_[lsite] == spins_[rsite] && random_01() < p_){
      union_find::unify(nodes, lsite, rsite);
      ++activated;
    }
  }
  const int nc = union_find::clusterize(nodes);

  LOGGING(logger_, "flip spins");

  double ma = 0.0;
  std::vector<int> cluster_size(nc);
  std::vector<int> cluster_spin(nc);
  for(int cl=0; cl<nc; ++cl){
    cluster_spin[cl] = q_ * random_01();
  }
  for(int site=0; site<nsites_; ++site){
    const int id = union_find::cluster_id(nodes, site);
    ++cluster_size[id];
    spins_[site] = cluster_spin[id];
    ma += cluster_spin[id];
  }

  if(!is_thermalized()){
    return;
  }

  LOGGING(logger_, "measure");

  const double invV = 1.0 / nsites_;
  double m2 = 0.0;
  double m4 = 0.0;
  for(int ci=0; ci < nc; ++ci){
    const double cl2 = util::pow2(cluster_size[ci] * invV);
    m4 += m4_coeff_ * cl2*cl2;
    m4 += 6 * cl2 * m2;
    m2 += m2_coeff_ * cl2;
  }

  double ene = ene_const_ - ene_coeff_ * activated;

  obs["Magnetization^2"] << m2;
  obs["Magnetization^4"] << m4;
  obs["|Magnetization|"] << std::abs(ma) * invV;
  obs["Susceptibility"] << m2 * nsites_ *  beta_;
  obs["Number of Sites"] << 1.0 * nsites_;
  obs["Energy"] << ene;
  obs["Activated Bonds"] << 1.0*activated;
  obs["Activated Bonds^2"] << 1.0*activated*activated;

  const double sec = tm.elapsed().wall * 1.0e-9;
  obs["Time"] << sec;
  obs["Speed"] << 1.0/sec;
}

void Evaluator::evaluate(alps::ObservableSet& obs) const
{
  const double T = alps::evaluate("T", params_);
  const double beta = 1.0/T;

  alps::RealObsevaluator m2 = obs["Magnetization^2"];
  alps::RealObsevaluator m4 = obs["Magnetization^4"];
  alps::RealObsevaluator ma = obs["|Magnetization|"];
  alps::RealObsevaluator nsites = obs["Number of Sites"];
  alps::RealObsevaluator n = obs["Activated Bonds"];
  alps::RealObsevaluator n2 = obs["Activated Bonds^2"];

  alps::RealObsevaluator binder = m4 / (m2*m2);
  binder.rename("Binder Ratio");
  obs.addObservable(binder);

  alps::RealObsevaluator csus = (m2 - ma * ma) * nsites * beta;
  csus.rename("Connected Susceptibility");
  obs.addObservable(csus);

  const double J = static_cast<double>(params_.value_or_default("J", 1.0)) * (params_["MODEL"]==std::string("Ising") ? 2.0:1.0);

  const double bj = -beta*J;
  const double coeff = bj / (boost::math::expm1(bj));
  alps::RealObsevaluator spec = (n2 - n*n + std::exp(bj)*n) * coeff * coeff;
  spec.rename("Specific Heat");
  obs.addObservable(spec);
}

} // end of namespace potts

