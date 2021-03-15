/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <boost/circular_buffer.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <gram_savitzky_golay/api.h>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

namespace gram_sg
{
using Vector6d = Eigen::Matrix<double, 6, 1>;

template<typename T>
struct GRAM_SAVITZKY_GOLAY_DLLAPI EigenVectorFilter
{
  EigenVectorFilter(const gram_sg::SavitzkyGolayFilterConfig & conf)
  : sg_conf(conf), sg_filter(conf), buffer(2 * sg_filter.config().m + 1)
  {
    reset(T::Zero());
  }

  void reset(const T & data)
  {
    buffer.clear();
    // Initialize to data
    for(size_t i = 0; i < buffer.capacity(); i++)
    {
      buffer.push_back(data);
    }
  }

  void reset()
  {
    reset(T::Zero());
  }

  void clear()
  {
    buffer.clear();
  }

  void add(const T & data)
  {
    buffer.push_back(data);
  }
  T filter() const
  {
    return sg_filter.filter(buffer);
  }

  const gram_sg::SavitzkyGolayFilterConfig & config() const
  {
    return sg_conf;
  }

  bool ready() const
  {
    return buffer.size() == buffer.capacity();
  }

protected:
  /** Filtering **/
  gram_sg::SavitzkyGolayFilterConfig sg_conf;
  gram_sg::SavitzkyGolayFilter sg_filter;
  // Buffers for Savitzky_golay
  boost::circular_buffer<T> buffer;
};

/**
 * Rotation Filter
 * Based on Peter Cork lecture here:
 * https://www.cvl.isy.liu.se/education/graduate/geometry2010/lectures/Lecture7b.pdf
 * Adapted to real time filtering through Savitzky-Golay
 **/
struct GRAM_SAVITZKY_GOLAY_DLLAPI RotationFilter
{
  RotationFilter(const gram_sg::SavitzkyGolayFilterConfig & conf);
  void reset(const Eigen::Matrix3d & r);
  void reset();
  void clear();
  void add(const Eigen::Matrix3d & r);
  Eigen::Matrix3d filter() const;
  bool ready() const
  {
    return buffer.size() == buffer.capacity();
  }

protected:
  /** Filtering **/
  gram_sg::SavitzkyGolayFilterConfig sg_conf;
  gram_sg::SavitzkyGolayFilter sg_filter;
  // Buffers for Savitzky_golay
  boost::circular_buffer<Eigen::Matrix3d> buffer;
};

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

template<typename TransformMatrixT = Eigen::Affine3d>
struct TransformFilter : public TransformFilterBase<TransformMatrixT>
{
  using ParentFilter = TransformFilterBase<TransformMatrixT>;
  using ParentFilter::TransformFilterBase;
};

/**
 * @brief Filters a velocity type
 * The default implementation expects a velocity expressed as a Vector6d (e.g linear and angular velocity or se3), and
 * uses EigenVectorFilter to perform the filtering
 */
template<typename VelocityT>
struct VelocityFilterBase
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  VelocityFilterBase(const gram_sg::SavitzkyGolayFilterConfig & conf) : vfilter(conf) {}

  void reset(const VelocityT & T)
  {
    vfilter.reset(T);
  }
  void reset()
  {
    vfilter.reset();
  }
  void add(const VelocityT & T)
  {
    vfilter.add(T);
  }
  VelocityT filter() const
  {
    return vfilter.filter();
  }

  const gram_sg::SavitzkyGolayFilterConfig & config() const
  {
    return vfilter.config();
  }
  bool ready() const
  {
    return vfilter.ready();
  }

protected:
  EigenVectorFilter<Vector6d> vfilter;
};

template<typename VelocityT>
struct VelocityFilter : public VelocityFilterBase<VelocityT>
{
  using VelocityFilterBase<VelocityT>::VelocityFilterBase;
};

} // namespace gram_sg
