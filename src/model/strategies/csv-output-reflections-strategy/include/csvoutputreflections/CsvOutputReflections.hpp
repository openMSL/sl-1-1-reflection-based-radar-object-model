//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef CSV_OUTPUT_REFLECTIONS_STRATEGY_HPP
#define CSV_OUTPUT_REFLECTIONS_STRATEGY_HPP

#include <string>

#include <model/include/strategy.hpp>

using namespace osi3;

namespace model
{

class CsvOutputReflections : public Strategy
{

    using Strategy::Strategy;

    void apply(SensorData& sensor_data) override;

    std::string file_path_reflections;
    bool first_call = true;

  public:
  private:
    static void write_first_line_to_csv(const std::string& path);
    static void
    write_data_to_csv(const std::string& path, double timestamp, size_t reflection_idx, double azimuth_in_deg, double elevation_in_deg, double distance, double signal_strength);
};

}  // namespace model

#endif  // CSV_OUTPUT_REFLECTIONS_STRATEGY_HPP
