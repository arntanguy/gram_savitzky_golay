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
  gram_sg::SavitzkyGolayFilterConfig config() const
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
 * @brief Filters Affine3d
 * The transformations are first converted to their translation and RPY
 * components, and then each component is filtered individually
 * Finally the result is converted back to an Affine3d
 */
struct GRAM_SAVITZKY_GOLAY_DLLAPI TransformFilter
{
  TransformFilter(const gram_sg::SavitzkyGolayFilterConfig & conf);
  void reset(const Eigen::Affine3d & T);
  void reset();
  void clear();
  void add(const Eigen::Affine3d & T);
  Eigen::Affine3d filter() const;
  gram_sg::SavitzkyGolayFilterConfig config() const
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

} // namespace gram_sg
