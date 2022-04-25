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

template<typename T>
struct EigenVectorFilter
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

} // namespace gram_sg
