//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef CSV_OUTPUT_DETECTIONS_STRATEGY_HPP
#define CSV_OUTPUT_DETECTIONS_STRATEGY_HPP

#include <model/include/strategy.hpp>
#include <string>

using namespace osi3;

namespace model {

    class CsvOutputDetections : public Strategy {

        using Strategy::Strategy;

        void apply(SensorData &) override;

        std::string file_path_detections;
        bool first_call = true;

    public:

    private:

        void write_first_line_to_CSV(const std::string &path);
        static void write_data_to_CSV(const std::string& path, double timestamp, size_t detection_idx, double azimuth_in_deg, double elevation_in_deg, double distance, double intensity);
    };

}

#endif //CSV_OUTPUT_DETECTIONS_STRATEGY_HPP
