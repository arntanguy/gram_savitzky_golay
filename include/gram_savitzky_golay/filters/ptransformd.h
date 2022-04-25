#pragma once

#include <SpaceVecAlg/SpaceVecAlg>
#include <gram_savitzky_golay/filters/transform.h>

namespace gram_sg
{

/**
 * @brief Filters PTransform
 * The transformations are first converted to their translation and RPY
 * compenents, and then each component is filtered individually
 * Finally the result is converted back to a PTransform
 */
template<>
struct TransformFilter<sva::PTransformd> : public TransformFilterBase<sva::PTransformd>
{
  using ParentFilter = TransformFilterBase<sva::PTransformd>;
  using ParentFilter::TransformFilterBase;
  sva::PTransformd filter() const
  {
    const Eigen::Vector3d & trans_res = ParentFilter::trans_filter.filter();
    const Eigen::Matrix3d & rot_res = ParentFilter::rot_filter.filter();
    return {rot_res, trans_res};
  }
};

using PTransformdFilter = TransformFilter<sva::PTransformd>;

} // namespace gram_sg
