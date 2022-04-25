#pragma once

#include <SpaceVecAlg/SpaceVecAlg>
#include <gram_savitzky_golay/filters/velocity.h>

namespace gram_sg
{
template<>
struct VelocityFilter<sva::MotionVecd> : public VelocityFilterBase<Eigen::Vector6d>
{
  using ParentFilter = VelocityFilterBase<Eigen::Vector6d>;
  using ParentFilter::VelocityFilterBase;
  void reset(const sva::MotionVecd & T)
  {
    ParentFilter::reset(T.vector());
  }
  void add(const sva::MotionVecd & T)
  {
    ParentFilter::add(T.vector());
  }
  sva::MotionVecd filter() const
  {
    return ParentFilter::filter();
  }
};

using MotionVecdFilter = VelocityFilter<sva::MotionVecd>;

} // namespace gram_sg
