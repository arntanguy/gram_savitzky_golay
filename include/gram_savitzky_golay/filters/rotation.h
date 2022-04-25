/*
 * Copyright 2017-2018 CNRS-UM LIRMM
 * Copyright 2019-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once
#include <boost/circular_buffer.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <gram_savitzky_golay/gram_savitzky_golay.h>

namespace gram_sg
{

/**
 * Rotation Filter
 * Based on Peter Cork lecture here:
 * https://www.cvl.isy.liu.se/education/graduate/geometry2010/lectures/Lecture7b.pdf
 * Adapted to real time filtering through Savitzky-Golay
 **/
template<typename RotationMatrixT>
struct GRAM_SAVITZKY_GOLAY_DLLAPI RotationFilterBase
{
  RotationFilterBase(const gram_sg::SavitzkyGolayFilterConfig & conf)
  : sg_conf(conf), sg_filter(conf), buffer(2 * sg_filter.config().m + 1)
  {
    reset(RotationMatrixT::Zero());
  }

  void reset(const RotationMatrixT & r)
  {
    buffer.clear();
    // Initialize to data
    for(size_t i = 0; i < buffer.capacity(); i++)
    {
      buffer.push_back(r);
    }
  }

  void reset()
  {
    reset(RotationMatrixT::Zero());
  }

  void clear()
  {
    buffer.clear();
  }

  void add(const RotationMatrixT & r)
  {
    buffer.push_back(r);
  }

  RotationMatrixT filter() const
  {
    // Apply a temporal (savitzky-golay) convolution,
    // followed by an orthogonalization
    const RotationMatrixT & result = sg_filter.filter(buffer);
    Eigen::JacobiSVD<RotationMatrixT> svd(result, Eigen::ComputeFullV | Eigen::ComputeFullU);
    RotationMatrixT res = svd.matrixU() * svd.matrixV().transpose();
    return res;
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
  boost::circular_buffer<RotationMatrixT> buffer;
};

using RotationFilter = RotationFilterBase<Eigen::Matrix3d>;

} // namespace gram_sg
