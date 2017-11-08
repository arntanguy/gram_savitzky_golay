#include "gram_savitzky_golay/gram_savitzky_golay.h"
#include <iostream>

namespace gram_sg
{
// OK
double GramPoly(const int i, const int m, const int k, const int s)
{
  if (k > 0)
  {
    return (4. * k - 2.) / (k * (2. * m - k + 1.)) * (i * GramPoly(i, m, k - 1, s) + s * GramPoly(i, m, k - 1, s - 1)) - ((k - 1.) * (2. * m + k)) / (k * (2. * m - k + 1.)) * GramPoly(i, m, k - 2, s);
  }
  else
  {
    if (k == 0 && s == 0)
      return 1.;
    else
      return 0.;
  }
}

// OK
double GenFact(const int a, const int b)
{
  double gf = 1.;

  for (int j = (a - b) + 1; j <= a; j++)
  {
    gf *= j;
  }
  return gf;
}

double Weight(const int i, const int t, const int m, const int n, const int s)
{
  double w = 0;
  for (int k = 0; k <= n; ++k)
  {
    w = w + (2 * k + 1) * (GenFact(2 * m, k) / GenFact(2 * m + k + 1, k + 1)) * GramPoly(i, m, k, 0) * GramPoly(t, m, k, s);
  }
  return w;
}

std::vector<double> ComputeWeights(const int m, const int t, const int n, const int s)
{
  std::vector<double> weights(2 * m + 1);
  for (int i = 0; i < 2 * m + 1; ++i)
  {
    weights[i] = Weight(i - m, t, m, n, s);
  }
  return weights;
}

SavitzkyGolayFilter::SavitzkyGolayFilter(const int m, const int t, const int n, const int s) : conf_(m, n, s, t)
{
  // Compute weights for the time window 2*m+1, for the t'th least-square
  // point of the s'th derivative
  weights_ = ComputeWeights(m, t, n, s);
}

double SavitzkyGolayFilter::filter(const std::vector<double>& x)
{
  assert(x.size() == weights_.size());
  double res = 0;
  for (unsigned i = 0; i < x.size(); i++)
  {
    res += weights_[i] * x[i];
  }
  return res;
}
} /* gram_sg */
