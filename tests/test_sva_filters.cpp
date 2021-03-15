/*
 * Copyright 2017-2020 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#define BOOST_TEST_MODULE TestSVAFilters
#include <boost/test/unit_test.hpp>
#include <gram_savitzky_golay/sva_filters.h>

#ifndef M_PI
#  include <boost/math/constants/constants.hpp>
#  define M_PI boost::math::constants::pi<double>()
#endif

using namespace gram_sg;

BOOST_AUTO_TEST_CASE(test_transform_filter)
{
  /** Filtering **/
  {
    gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
    TransformFilter<sva::PTransformd> filter(sg_conf);

    Eigen::Matrix3d rot_sg;
    rot_sg = Eigen::AngleAxisd(1.2, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(1.9, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(2.3, Eigen::Vector3d::UnitZ());
    auto X_init = sva::PTransformd{rot_sg, Eigen::Vector3d{0.1, 5, 0.3}};
    filter.reset(X_init);
    auto res_init = filter.filter();
    BOOST_CHECK_MESSAGE(X_init.matrix().isApprox(res_init.matrix()), "\n" << X_init.matrix() << "\n-----\n"
                                                                          << res_init.matrix());
  }

  {
    gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
    TransformFilter<sva::PTransformd> filter(sg_conf);
    Eigen::Matrix3d rot_sg;
    rot_sg = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-M_PI / 2, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ());
    auto X_init = sva::PTransformd{rot_sg, Eigen::Vector3d{0.1, 5, 0.3}};
    filter.reset(X_init);

    auto res_init = filter.filter();
    BOOST_CHECK_MESSAGE(X_init.matrix().isApprox(res_init.matrix()), "\n" << X_init.matrix() << "\n-----\n"
                                                                          << res_init.matrix());
  }
}

BOOST_AUTO_TEST_CASE(test_velocity_filter)
{
  gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
  VelocityFilter<sva::MotionVecd> filter(sg_conf);
  sva::MotionVecd Vinit{{0.1, 0.2, 0.3}, {0.1, 0.20, 0.3}};
  filter.reset(Vinit);
  auto res_init = filter.filter();
  BOOST_CHECK_MESSAGE(Vinit.vector().isApprox(res_init.vector()), "\n" << Vinit.vector() << "\n-----\n"
                                                                       << res_init.vector());
}
