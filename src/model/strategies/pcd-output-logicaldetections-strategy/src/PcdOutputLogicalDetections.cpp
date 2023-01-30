//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#include "pcdoutputlogicaldetections/PcdOutputLogicalDetections.hpp"
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

void model::PcdOutputLogicalDetections::apply(SensorData &sensor_data) {
    log("Starting .pcd output for logical detections");

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
    double timestamp = (double)time_seconds + (double) time_nanos / 1000000000;

    if (sensor_data.logical_detection_data().logical_detection_size() == 0) {
        log("No logical detections for .csv output at timestamp " + std::to_string(timestamp));
        return;
    }

    if (first_call) {
        #include <pcdoutputlogicaldetections/set_pcd_file_path_logicaldetections.cpp>
        first_call = false;
    }
    
    /// Header needed in every file
    std::string filename = "LogicalDetections_";
    filename.append(std::to_string(timestamp));
    filename.append(".pcd");
    #if defined(_WIN32)
        std::string path = path_string + "\\" + filename;
    #else
        std::string path = path_string + "/" + filename;
    #endif
    writePcdHeader(path, sensor_data);

    for (const auto &logical_detection: sensor_data.logical_detection_data().logical_detection()) {
        write2Pcd(path, float(logical_detection.position().x()), float(logical_detection.position().y()), float(logical_detection.position().z()),
                    float(logical_detection.intensity()));
    }
}

void PcdOutputLogicalDetections::writePcdHeader(const std::string& path, const SensorData& sensor_data) {
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << "# .PCD v0.7 - Point Cloud Data file format" << std::endl;
    my_file << "VERSION " << "0.7" << std::endl;
    my_file << "FIELDS " << "x y z intensity" << std::endl;
    my_file << "SIZE " << "4 4 4 4" << std::endl;
    my_file << "TYPE " << "F F F F" << std::endl;
    my_file << "COUNT " << "1 1 1 1" << std::endl;
    my_file << "WIDTH " << sensor_data.logical_detection_data().logical_detection_size() << std::endl;
    my_file << "HEIGHT " << "1" << std::endl;
    my_file << "VIEWPOINT " << "0 0 0 1 0 0 0" << std::endl;
    my_file << "POINTS " << sensor_data.logical_detection_data().logical_detection_size() << std::endl;
    my_file << "DATA " << "ascii" << std::endl;
    my_file.close();
}

void PcdOutputLogicalDetections::write2Pcd(const std::string& path, float x, float y, float z, float intensity) {
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << x << " " << y << " " << z << " " << intensity << std::endl;
    my_file.close();
}
