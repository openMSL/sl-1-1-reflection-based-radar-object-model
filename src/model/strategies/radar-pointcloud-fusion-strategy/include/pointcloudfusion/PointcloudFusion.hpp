//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef POINTCLOUD_FUSION_STRATEGY_HPP
#define POINTCLOUD_FUSION_STRATEGY_HPP

#include <model/include/strategy.hpp>

using namespace osi3;

namespace model
{

class PointcloudFusion : public Strategy
{

    using Strategy::Strategy;
    void apply(SensorData& sensor_data) override;

  public:
  private:
    static void calculate_fused_pointcloud_for_given_lidar_sensors(SensorData& sensor_data);
    static void calculate_fused_pointcloud_for_given_radar_sensors(SensorData& sensor_data);
};
}  // namespace model

#endif  // POINTCLOUD_FUSION_STRATEGY_HPP
