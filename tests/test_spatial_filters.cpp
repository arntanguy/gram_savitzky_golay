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

#define BOOST_TEST_MODULE TestSpatialFilters
#include <boost/test/unit_test.hpp>
#include <gram_savitzky_golay/spatial_filters.h>

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
    TransformFilter filter(sg_conf);

    Eigen::Matrix3d rot_sg;
    rot_sg = Eigen::AngleAxisd(1.2, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(1.9, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(2.3, Eigen::Vector3d::UnitZ());
    Eigen::Affine3d X_init = Eigen::Affine3d::Identity();
    X_init.matrix().block<3, 3>(0, 0) = rot_sg;
    X_init.matrix().block<3, 1>(0, 3) = Eigen::Vector3d{0.1, 5, 0.3};
    filter.reset(X_init);
    const auto & res_init = filter.filter();
    BOOST_CHECK_MESSAGE(X_init.matrix().isApprox(res_init.matrix()), "\n" << X_init.matrix() << "\n-----\n"
                                                                          << res_init.matrix());
  }

  {
    gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
    TransformFilter filter(sg_conf);
    Eigen::Matrix3d rot_sg;
    rot_sg = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(-M_PI / 2, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ());
    Eigen::Affine3d X_init = Eigen::Affine3d::Identity();
    X_init.matrix().block<3, 3>(0, 0) = rot_sg;
    X_init.matrix().block<3, 1>(0, 3) = Eigen::Vector3d{0.1, 5, 0.3};
    filter.reset(X_init);

    const auto & res_init = filter.filter();
    BOOST_CHECK_MESSAGE(X_init.matrix().isApprox(res_init.matrix()), "\n" << X_init.matrix() << "\n-----\n"
                                                                          << res_init.matrix());
  }
}

BOOST_AUTO_TEST_CASE(test_rotation_filter)
{
  Eigen::Matrix3d rot;
  rot = Eigen::AngleAxisd(1.2, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitZ());

  gram_sg::SavitzkyGolayFilterConfig sg_conf(50, 50, 2, 0);
  RotationFilter filter(sg_conf);

  filter.reset(rot);

  const Eigen::Matrix3d res = filter.filter();
  BOOST_CHECK(rot.isApprox(res));
}
