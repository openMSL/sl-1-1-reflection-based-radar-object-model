//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef CSV_OUTPUT_LOGICALDETECTIOONS_STRATEGY_HPP
#define CSV_OUTPUT_LOGICALDETECTIOONS_STRATEGY_HPP

#include <model/include/strategy.hpp>
#include <string>

using namespace osi3;

namespace model {

	class CsvOutputLogicalDetections : public Strategy {

        using Strategy::Strategy;

        void apply(SensorData &) override;

        std::string file_path_logicaldetections;
        bool first_call = true;

    public:

    private:
        static void write_first_line_to_CSV(const std::string &path);
        static void write_data_to_CSV(const std::string& path, double timestamp, size_t detection_idx, float x, float y, float z, float intensity);
    };

}

#endif //CSV_OUTPUT_LOGICALDETECTIOONS_STRATEGY_HPP
