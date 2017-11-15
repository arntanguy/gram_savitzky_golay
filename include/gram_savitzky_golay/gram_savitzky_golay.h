#pragma once

#include <Eigen/Core>
#include <cassert>
#include <cstddef>
#include <sstream>
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
  //! Time at which the filter is applied
  // For real-time, should be t=m
  int t;
  //! Polynomial order
  int n;
  //! Derivation order (0 for no derivation)
  int s;

  SavitzkyGolayFilterConfig(const int m, const int t, const int n, const int s) : m(m), t(t), n(n), s(s)
  {
  }

  friend std::ostream& operator<<(std::ostream& os, const SavitzkyGolayFilterConfig& conf);
};

class SavitzkyGolayFilter
{
 public:
 private:
  const SavitzkyGolayFilterConfig conf_;
  std::vector<double> weights_;
  void init();

 public:
  SavitzkyGolayFilter(const int m, const int t, const int n, const int s);
  SavitzkyGolayFilter(const SavitzkyGolayFilterConfig& conf);

  /**
   * @brief Apply Savitzky-Golay convolution to the data x should have size 2*m+1
   * As the function only applies a convolution, it runs in O(2m+1), and should
   * be rather fast.
   *
   * @param v Container with the data to be filtered.
   * Should have 2*m+1 elements
   * Type of elements needs to be compatible with multiplication by a scalar,
   * and addition with itself
   * Common types would be std::vector<double>, std::vector<Eigen::VectorXd>, boost::circular_buffer<Eigen::Vector3d>...
   *
   * @param initial_value Initial value for the filter's accumulator. It should be the zero.
   *
   * @return Filtered value according to the precomputed filter weights.
   */
  template <typename ContainerT>
  typename ContainerT::value_type filter(const ContainerT& v, typename ContainerT::value_type&& initial_value) const
  {
    assert(v.size() == weights_.size());
    using T = typename ContainerT::value_type;
    T res = initial_value;
    int i = 0;
    for (const T& value : v)
    {
      res += weights_[i] * value;
      ++i;
    }
    return res;
  }

  std::vector<double> weights() const
  {
    return weights_;
  }
  SavitzkyGolayFilterConfig config() const
  {
    return conf_;
  }
};

} /* gram_sg */
