#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include <cmath>

namespace util{

class Observable{
public:
  Observable():num_(0),sum_(0.0),sum2_(0.0){}

  void reset() { num_ = 0; sum_ = sum2_ = 0.0;}
  Observable& add(double v)
  {
    ++num_;
    sum_ += v;
    sum2_ += v*v;
    return *this;
  }
  Observable& operator<<(double v) { return add(v);}
  
  long int num() const { return num_;}
  double mean() const { return num_ > 0 ? sum_ / num_ : 0.0/0.0; }
  double var() const{
    return num_ == 0 ? 0.0/0.0
         : num_ == 1 ? 1.0/0.0
         :             (sum2_ - sum_*sum_/num_)/(num_-1);
  }
  double error() const { return std::sqrt(var());}

private:
  long int num_;
  double sum_;
  double sum2_;
};

} // end of namespace util

#endif
