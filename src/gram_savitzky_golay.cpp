/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <gram_savitzky_golay/gram_savitzky_golay.h>
#include <iostream>

namespace gram_sg
{
// OK
double GramPoly(const int i, const int m, const int k, const int s)
{
  if(k > 0)
  {
    return (4. * k - 2.) / (k * (2. * m - k + 1.)) * (i * GramPoly(i, m, k - 1, s) + s * GramPoly(i, m, k - 1, s - 1))
           - ((k - 1.) * (2. * m + k)) / (k * (2. * m - k + 1.)) * GramPoly(i, m, k - 2, s);
  }
  else
  {
    if(k == 0 && s == 0)
      return 1.;
    else
      return 0.;
  }
}

// OK
double GenFact(const int a, const int b)
{
  double gf = 1.;

  for(int j = (a - b) + 1; j <= a; j++)
  {
    gf *= j;
  }
  return gf;
}

double Weight(const int i, const int t, const int m, const int n, const int s)
{
  double w = 0;
  for(int k = 0; k <= n; ++k)
  {
    w = w
        + (2 * k + 1) * (GenFact(2 * m, k) / GenFact(2 * m + k + 1, k + 1)) * GramPoly(i, m, k, 0)
              * GramPoly(t, m, k, s);
  }
  return w;
}

std::vector<double> ComputeWeights(const int m, const int t, const int n, const int s)
{
  std::vector<double> weights(2 * static_cast<size_t>(m) + 1);
  for(int i = 0; i < 2 * m + 1; ++i)
  {
    weights[static_cast<size_t>(i)] = Weight(i - m, t, m, n, s);
  }
  return weights;
}

void SavitzkyGolayFilter::init()
{
  // Compute weights for the time window 2*m+1, for the t'th least-square
  // point of the s'th derivative
  weights_ = ComputeWeights(static_cast<int>(conf_.m), conf_.t, static_cast<int>(conf_.n), static_cast<int>(conf_.s));
  dt_ = std::pow(conf_.time_step(), conf_.derivation_order());
}

std::ostream & operator<<(std::ostream & os, const SavitzkyGolayFilterConfig & conf)
{
  os << "m                       : " << conf.m << std::endl
     << "Window Size (2*m+1)     : " << 2 * conf.m + 1 << std::endl
     << "n (Order)               :" << conf.n << std::endl
     << "s (Differentiate)       : " << conf.s << std::endl
     << "t: Filter point ([-m,m]): " << conf.t << std::endl;
  return os;
}

} // namespace gram_sg
