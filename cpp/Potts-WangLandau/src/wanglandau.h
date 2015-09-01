#ifndef WANGLANDAU_H
#define WANGLANDAU_H

#include <vector>
#include <cmath>
#include <string>

namespace wanglandau{

class WangLandau{
public:
  WangLandau(int hsize, double flatness,
             double final_factor, int initial_threshold,
             double normalization, int normalization_origin);

  void visit(int index) { ++hist_[index]; g_[index] += factor_;}
  double prob(int now, int candidate) const { return std::exp( g_[now] - g_[candidate]); }
  bool finished() const { return factor_ < final_; }
  void update(bool verbose);
private:
  int hsize_;
  std::vector<long> hist_;
  std::vector<double> g_;
  double flatness_;
  double factor_;
  double final_;
  int threshold_;
  int stage_;
  double normalization_;
  int normalization_origin_;

  bool check_flat(bool verbose);
  void print(std::string const& filename) const;
  void reset();
};

} // end of namespace wanglandau

#endif
