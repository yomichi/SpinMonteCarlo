#include "wanglandau.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <algorithm>

namespace wanglandau{

WangLandau::WangLandau(int hsize, double flatness,
                       double final_factor, int initial_threshold,
                       double normalization, int normalization_origin)
  :hsize_(hsize), hist_(hsize), g_(hsize),
   flatness_(flatness), factor_(1.0), final_(final_factor),
   threshold_(initial_threshold), stage_(1),
   normalization_(normalization), normalization_origin_(normalization_origin)
{
}

bool WangLandau::check_flat(bool verbose = false)
{
  long hmin = std::numeric_limits<long>::max();
  int count = 0;
  long hsum = 0;
  for(int i=0; i<hsize_; ++i){
    if(hist_[i] == 0) continue;
    ++count;
    hsum += hist_[i];
    hmin = std::min(hmin, hist_[i]);
  }
  const double mean = static_cast<double>(hsum)/count;
  const double fn = hmin / mean;
  if(count < threshold_ || fn < flatness_){
    if(verbose){
      std::cout << "Histogram is not flat. "
               << "[ hmin = " << hmin
               << ", mean = " << mean
               << ", flatness = " << fn
               << ", nonzero bin = " << count
               << ", threshold = " << threshold_
               << ", factor = " << factor_
               << "]" << std::endl;
    }
    return false;
  }else{
    if(verbose){
      std::cout << "Histogram is     flat. "
               << "[ hmin = " << hmin
               << ", mean = " << mean
               << ", flatness = " << fn
               << ", nonzero bin = " << count
               << ", threshold = " << threshold_
               << ", factor = " << factor_
               << "]" << std::endl;
    }
    threshold_ = count;
    return true;
  }
}

void WangLandau::update(bool verbose=false)
{
  bool flat = check_flat(verbose);
  if(flat){

    /*
     * save and reset
     */

    std::string filename("hist-");
    filename += boost::lexical_cast<std::string>(stage_);
    filename += ".dat";
    std::ofstream ofs(filename.c_str());
    ofs << "# $1 : index" << std::endl;
    ofs << "# $2 : log of DoS" << std::endl;
    ofs << "# $3 : population" << std::endl;

    const double offset = normalization_ - g_[normalization_origin_];
    for(int i=0; i<hsize_; ++i){
      g_[i] += offset;
      if(hist_[i] != 0)
        ofs << i << " " << g_[i] << " " << hist_[i] << std::endl;
      hist_[i] = 0;
    }

    std::cout << "# " <<  stage_ << " iteration (factor = " << factor_ << " ) finished." << std::endl;
    factor_ *= 0.5;
    ++stage_;
  }
}

} // end of namespace wanglandau
