//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#include "csvoutputlogicaldetections/CsvOutputLogicalDetections.hpp"
#include <fstream>
#include <iostream>
#include <ctime>

#ifdef _WIN32
    #include <math.h>
    #include <direct.h>
#else
    #include <sys/stat.h>
#endif

using namespace model;
using namespace osi3;

static bool first_call = true;

void model::CsvOutputLogicalDetections::apply(SensorData &sensor_data) {
    log("Starting .csv output for logical detections");

    if ((sensor_data.sensor_view_size() == 0) || (!sensor_data.has_feature_data())) {
        log("No sensor view or feature data received");
        return;
    }

    if (!sensor_data.sensor_view(0).has_global_ground_truth()) {
        log("No global ground truth received");
        return;
    }

    auto time_nanos = sensor_data.sensor_view(0).global_ground_truth().timestamp().nanos();
    auto time_seconds = sensor_data.sensor_view(0).global_ground_truth().timestamp().seconds();
    double timestamp = (float)time_seconds + (double) time_nanos / 1000000000;

    if (sensor_data.logical_detection_data().logical_detection_size() == 0) {
        log("No logical detections for .csv output at timestamp " + std::to_string(timestamp));
        return;
    }

    /// Write header line of .csv on first call
    if (first_call){
        #include <csvoutputlogicaldetections/set_csv_file_path_logicaldetections.cpp>
        if (sensor_data.sensor_view(0).radar_sensor_view_size() > 0) {
            write_first_line_to_CSV(file_path_logicaldetections);
            first_call = false;
        } else {
            log("RadarSensorView empty at timestamp " + std::to_string(timestamp));
        }
    }

    size_t logical_detection_idx = 0;
    for (const auto &logical_detection: sensor_data.logical_detection_data().logical_detection()) {
        write_data_to_CSV(file_path_logicaldetections, timestamp, logical_detection_idx,
                        float(logical_detection.position().x()), float(logical_detection.position().y()), float(logical_detection.position().z()),
						float(logical_detection.intensity()));
        logical_detection_idx++;
    }
}

void CsvOutputLogicalDetections::write_first_line_to_CSV(const std::string &path) {
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << "timestamp_in_s, detection_id, x_in_m, y_in_m, z_in_m, intensity_in_%" << std::endl;
    my_file.close();
}

void CsvOutputLogicalDetections::write_data_to_CSV(const std::string& path, double timestamp, size_t detection_idx, float x, float y, float z, float intensity) {
    // Debug output of Data
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << timestamp << ", " << detection_idx << ", " << x << ", " << y << ", " << z << ", " << intensity << std::endl;
    my_file.close();
}
