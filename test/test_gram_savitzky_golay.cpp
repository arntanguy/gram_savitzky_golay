// Copyright 2017-2018 CNRS-UM LIRMM
// Copyright 2017-2018 Arnaud TANGUY <arnaud.tanguy@lirmm.fr>
//
// This file is part of robcalib.
//
// robcalib is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// robcalib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with robcalib.  If not, see <http://www.gnu.org/licenses/>.

//Link to Boost
#define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include <chrono>
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


BOOST_AUTO_TEST_CASE(TestIdentity)
{
  int m = 3;
  SavitzkyGolayFilter filter(m, 0, 2, 0);
  std::vector<double> data = {1, 1, 1, 1, 1, 1, 1};
  double res = filter.filter(data);
  BOOST_REQUIRE_CLOSE(res, 1, 10e-6);
}

BOOST_AUTO_TEST_CASE(TestRealTimeFilter)
{
  // Window size is 2*m+1
  const int m = 3;
  // Polynomial Order
  const int n = 2;
  // Initial Point Smoothing (ie evaluate polynomial at first point in the window)
  // Points are defined in range [-m;m]
  const int t = m;
  // Derivate? 0: no derivation, 1: first derivative...
  SavitzkyGolayFilter filter(m, t, n, 0);

  // Filter some data
  std::vector<double> data = {.1, .7, .9, .7, .8, .5, -.3};
  double result = filter.filter(data);
  double result_ref = -0.22619047619047616;
  BOOST_REQUIRE_CLOSE(result, result_ref, 10e-6);
}

BOOST_AUTO_TEST_CASE(TestRealTimeDerivative)
{
  // Window size is 2*m+1
  const int m = 3;
  // Polynomial Order
  const int n = 2;
  // Initial Point Smoothing (ie evaluate polynomial at first point in the window)
  // Points are defined in range [-m;m]
  const int t = m;


  // Test First Order Derivative
  SavitzkyGolayFilter filter(m, t, n, 1);
  SavitzkyGolayFilter filter_dt(m, t, n, 1, 0.005);
  BOOST_REQUIRE(filter_dt.config().time_step() == 0.005);

  // Filter some data
  std::vector<double> data = {.1, .2, .3, .4, .5, .6, .7};
  double result = filter.filter(data);
  double result_ref = 0.1;
  BOOST_REQUIRE_CLOSE(result, result_ref, 10e-6);

  // Test filtering with timestep=0.005
  data = {.1, .2, .3, .4, .5, .6, .7};
  result = filter_dt.filter(data);
  result_ref = 0.1/filter_dt.config().time_step();
  BOOST_REQUIRE_CLOSE(result, result_ref, 10e-6);

  // Filter some data
  data = {-1, -2, -3, -4, -5, -6, -7};
  result = filter.filter(data);
  result_ref = -1;
  BOOST_REQUIRE_CLOSE(result, result_ref, 10e-6);
  // Test filtering with timestep=0.005
  result = filter_dt.filter(data);
  result_ref = -1./filter_dt.config().time_step();
  BOOST_REQUIRE_CLOSE(result, result_ref, 10e-6);



  // Test Second Order Derivative
  SavitzkyGolayFilter second_order_filter(m, t, n, 2);

  // Filter some data
  data = {.1, .2, .3, .4, .5, .6, .7};
  result = second_order_filter.filter(data);
  BOOST_CHECK_SMALL(result, 10e-6);

  // Filter some data
  data = {-1, -2, -3, -4, -5, -6, -7};
  result = second_order_filter.filter(data);
  BOOST_CHECK_SMALL(result, 10e-6);
}

// Test derivation on a known polynomial function
BOOST_AUTO_TEST_CASE(TestPolynomialDerivative)
{
  // Polynomial is a*x^3 + bx^2 + c*x^1 + d
  double a = 10;
  double b = 2;
  double c = -3;
  double d = -4;
  double timeStep = 0.42;


  // Window size is 2*m+1
  const int m = 50;
  // Polynomial Order
  const int n = 3;
  // Points are defined in range [-m;m]
  // Eval at central point
  const int t = 0;

  SavitzkyGolayFilter filter_order1(m, t, n, 1, timeStep);
  SavitzkyGolayFilter filter_order2(m, t, n, 2, timeStep);
  std::vector<double> data;
  std::vector<double> derivative_order1, derivative_order2;
  data.resize(2*m+1);
  derivative_order1.resize(2*m+1);
  derivative_order2.resize(2*m+1);
  // Generate some data points
  for (unsigned x = 0; x < data.size(); ++x) {
    data[x] = a * std::pow(x, 3) + b * std::pow(x, 2) + c * std::pow(x, 1) + d;
    derivative_order1[x] = (3*a*std::pow(x, 2) + 2*b*std::pow(x,1) + c)/timeStep;
    derivative_order2[x] = (6*a*std::pow(x, 1) + 2*b) / std::pow(timeStep,2);
  }
  const auto result_order1 = filter_order1.filter(data);
  const auto expected_result_order1 = derivative_order1[m];
  const auto result_order2 = filter_order2.filter(data);
  const auto expected_result_order2 = derivative_order2[m];

  BOOST_REQUIRE_CLOSE(result_order1, expected_result_order1, 10e-8);
  BOOST_REQUIRE_CLOSE(result_order2, expected_result_order2, 10e-8);
}

BOOST_AUTO_TEST_CASE(FilterSpeed)
{
    int m = 10000;
    SavitzkyGolayFilter filter(m, 0, 2, 0);
    std::vector<double> data(2*m+1,1.);

    double totalTime = 0;
    int Nsample = 1000;
    for (int i = 0; i < Nsample; ++i) {
      auto start_time = std::chrono::high_resolution_clock::now();
      filter.filter(data);
      auto end_time = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
      totalTime += elapsed.count();
    }
    totalTime /= Nsample;
    std::cout << "Filtering performed in " << totalTime << " (ms)" << std::endl;;

    BOOST_REQUIRE(totalTime < 0.0001);
}
