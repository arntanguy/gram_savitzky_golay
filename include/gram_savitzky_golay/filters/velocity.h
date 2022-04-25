/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <gram_savitzky_golay/filters/eigen_vector.h>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

namespace gram_sg
{

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
