/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <Eigen/Core>
#include <gram_savitzky_golay/filters/eigen_vector.h>
#include <gram_savitzky_golay/filters/rotation.h>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

namespace gram_sg
{

/**
 * @brief Filters homogeneous transformations
 * Translation and rotation are filtered independently:
 * - translation is filtered using the EigenVectorFilter
 * - rotation is filtered using RotationFilter
 *
 * They are then recombined into a TransformMatrixT
 *
 * \tparm TransformMatrixT Homogeneous matrix 4x4 representing a spatial transformation (translation + rotation). Common
 * types include Eigen::Matrix3d, Eigen::Affine3d Custom matrix types should support eigen-like: translation(),
 * rotation(), Identity(), block() and arithmetic operators
 */
template<typename TransformMatrixT>
struct TransformFilterBase
{
  TransformFilterBase(const gram_sg::SavitzkyGolayFilterConfig & conf) : trans_filter(conf), rot_filter(conf) {}

  void reset(const TransformMatrixT & T)
  {
    trans_filter.reset(T.translation());
    rot_filter.reset(T.rotation());
  }

  void reset()
  {
    trans_filter.reset();
    rot_filter.reset();
  }

  void clear()
  {
    trans_filter.clear();
    rot_filter.clear();
  }

  void add(const TransformMatrixT & T)
  {
    trans_filter.add(T.translation());
    rot_filter.add(T.rotation());
  }

  TransformMatrixT filter() const
  {
    TransformMatrixT res = TransformMatrixT::Identity();
    res.matrix().template block<3, 3>(0, 0) = rot_filter.filter();
    res.matrix().template block<3, 1>(0, 3) = trans_filter.filter();
    return res;
  }

  const gram_sg::SavitzkyGolayFilterConfig & config() const
  {
    return trans_filter.config();
  }

  bool ready() const
  {
    return trans_filter.ready() && rot_filter.ready();
  }

protected:
  EigenVectorFilter<Eigen::Vector3d> trans_filter;
  RotationFilter rot_filter;
};

template<typename TransformMatrixT>
struct TransformFilter : public TransformFilterBase<TransformMatrixT>
{
  using ParentFilter = TransformFilterBase<TransformMatrixT>;
  using ParentFilter::TransformFilterBase;
};

} // namespace gram_sg
