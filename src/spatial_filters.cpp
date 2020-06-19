//  Copyright 2017-2018 CNRS-UM LIRMM
//  Copyright 2017-2018 Arnaud TANGUY <arnaud.tanguy@lirmm.fr>
//
//  This file is part of robcalib.
//
//  robcalib is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  robcalib is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with robcalib.  If not, see <http://www.gnu.org/licenses/>.

#include "gram_savitzky_golay/spatial_filters.h"

#include <boost/circular_buffer.hpp>
#include <Eigen/SVD>

namespace gram_sg
{
RotationFilter::RotationFilter(const gram_sg::SavitzkyGolayFilterConfig & conf)
: sg_conf(conf), sg_filter(conf), buffer(2 * sg_filter.config().m + 1)

{
  reset(Eigen::Matrix3d::Zero());
}

void RotationFilter::reset(const Eigen::Matrix3d & r)
{
  buffer.clear();
  // Initialize to data
  for(size_t i = 0; i < buffer.capacity(); i++)
  {
    buffer.push_back(r);
  }
}

void RotationFilter::reset()
{
  RotationFilter::reset(Eigen::Matrix3d::Zero());
}

void RotationFilter::clear()
{
  buffer.clear();
}

void RotationFilter::add(const Eigen::Matrix3d & r)
{
  buffer.push_back(r);
}
Eigen::Matrix3d RotationFilter::filter() const
{
  // Apply a temporal (savitzky-golay) convolution,
  // followed by an orthogonalization
  const Eigen::Matrix3d & result = sg_filter.filter(buffer);
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(result, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::Matrix3d res = svd.matrixU() * svd.matrixV().transpose();
  return res;
}

TransformFilter::TransformFilter(const gram_sg::SavitzkyGolayFilterConfig & conf) : trans_filter(conf), rot_filter(conf)
{
}

void TransformFilter::reset(const Eigen::Affine3d & T)
{
  trans_filter.reset(T.translation());
  rot_filter.reset(T.rotation());
}

void TransformFilter::reset()
{
  trans_filter.reset();
  rot_filter.reset();
}

void TransformFilter::clear()
{
  trans_filter.clear();
  rot_filter.clear();
}

void TransformFilter::add(const Eigen::Affine3d & T)
{
  trans_filter.add(T.translation());
  rot_filter.add(T.rotation());
}

Eigen::Affine3d TransformFilter::filter() const
{
  const Eigen::Vector3d & trans_res = trans_filter.filter();
  const Eigen::Matrix3d & rot_res = rot_filter.filter();
  Eigen::Matrix4d rot = Eigen::Matrix4d::Identity();
  rot.block<3, 3>(0, 0) = rot_res;
  rot.block<3, 1>(0, 3) = trans_res;
  return Eigen::Affine3d(rot);
}

} // namespace gram_sg
