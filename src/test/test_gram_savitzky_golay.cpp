//Link to Boost
#define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include "gram_savitzky_golay/gram_savitzky_golay.h"

using namespace gram_sg;

BOOST_AUTO_TEST_CASE(TestGorryTables)
{
  // Compare with tables in the paper from Gorry.
  // Convolution weights for quadratic initial-point smoothing:
  // polynomial order = 2, derivative = 0
  std::vector<double> sg7_gram{
      32,
      15,
      3,
      -4,
      -6,
      -3,
      5};

  int m = 3;
  SavitzkyGolayFilter filter(m, -m, 2, 0);

  const auto& filter_weights = filter.weights();
  for (unsigned int i = 0; i < sg7_gram.size(); i++)
  {
    std::cout << "ref: " << sg7_gram[i] << ", computed: " << filter_weights[i] * 42 << std::endl;
    BOOST_REQUIRE_CLOSE(sg7_gram[i], filter_weights[i] * 42, 10e-6);
  }

  // BOOST_CHECK( test_object.is_valid() );
}

BOOST_AUTO_TEST_CASE(TestGorryDerivative)
{
  // Compare with tables in the paper from Gorry.
  // Convolution weights for quadratic initial-point first derivative:
  // polynomial order = 2, derivative = 1
  std::vector<double> sg7_deriv_gram{
      -13,
      -2,
      5,
      8,
      7,
      2,
      -7};

  int m = 3;
  SavitzkyGolayFilter filter(m, -m, 2, 1);

  const auto& filter_weights = filter.weights();
  for (unsigned int i = 0; i < sg7_deriv_gram.size(); i++)
  {
    std::cout << "ref: " << sg7_deriv_gram[i] << ", computed: " << filter_weights[i] * 28 << std::endl;
    BOOST_REQUIRE_CLOSE(sg7_deriv_gram[i], filter_weights[i] * 28, 10e-6);
  }
}
