#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

namespace gram_sg
{
/*! GRAMPOLY Calculates the Gram Polynomial (s=0) or its sth derivative
  *  evaluated at i, order k, over 2m+1 points
  */
double GramPoly(const int i, const int m, const int k, const int s);

/*! GenFact Calculates the generalized factorial (a)(a-1)(a-2)...(a-b+1)
 *    Detailed explanation goes here
 */
double GenFact(const int a, const int b);

/*!
     * Weight Calculates the weight of the ith data point for the t'th
     * Least-Square point of the s'th derivative, over 2m+1 points, order n
     */
double Weight(const int i, const int t, const int m, const int n, const int s);

/*!
 * Computes the weights for each data point over the window of size 2*m+1,
 * evaluated at time t, with order n and derivative s
 */
std::vector<double> ComputeWeights(const int m, const int t, const int n, const int s);

struct SavitzkyGolayFilterConfig
{
  //! Window size is 2*m+1
  int m;
  //! Polynomial order
  int n;
  //! Derivation order (0 for no derivation)
  int s;
  //! Time at which the filter is applied
  // For real-time, should be t=m
  int t;

  SavitzkyGolayFilterConfig(const int m, const int n, const int s, const int t) : m(m), n(n), s(s), t(t)
  {
  }
};

class SavitzkyGolayFilter
{
 public:
 private:
  const SavitzkyGolayFilterConfig conf_;
  std::vector<double> weights_;

 public:
  SavitzkyGolayFilter(const int m, const int t, const int n, const int s);

  /*!
     * Apply Savitzky-Golay convolution to the data
     * x should have size 2*m+1
     */
  double filter(const std::vector<double>& x);

  std::vector<double> weights() const
  {
    return weights_;
  }
};

} /* gram_sg */
