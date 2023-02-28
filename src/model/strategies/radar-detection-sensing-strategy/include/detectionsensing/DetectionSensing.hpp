//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef DETECTION_SENSING_STRATEGY_HPP
#define DETECTION_SENSING_STRATEGY_HPP

#include <model/include/strategy.hpp>

using namespace osi3;

namespace model
{

class DetectionSensing : public Strategy
{

  private:
  public:
    //// Main Functions
    using Strategy::Strategy;
    void apply(SensorData& sensor_data) override;
    typedef std::vector<std::vector<std::vector<std::vector<std::vector<float> > > > > TypeRadarCuboid;
    typedef struct
    {
        std::vector<float> num_peaks;
        std::vector<float> peak_list_idx;
        std::vector<float> peak_list;
        std::vector<float> sub_bin_range;
        std::vector<float> sub_bin_doppler;
        std::vector<float> sub_bin_azimuth;
        std::vector<float> sub_bin_elevation;
        std::vector<float> detection_power;
    } TypeDetectionList;
};
}  // namespace model

#endif  // DETECTION_SENSING_STRATEGY_HPP
