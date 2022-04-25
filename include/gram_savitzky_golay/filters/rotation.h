/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once
#include <boost/circular_buffer.hpp>
#include <Eigen/Core>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

namespace gram_sg
{

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

} // namespace gram_sg
