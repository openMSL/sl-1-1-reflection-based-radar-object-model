//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "csvoutputreflections/CsvOutputReflections.hpp"

#include <ctime>
#include <fstream>
#include <iostream>

#ifdef _WIN32
#include <math.h>

#include <direct.h>
#else
#include <cmath>

#include <sys/stat.h>
#endif

const float speed_of_light = 299792458.0;

using namespace model;
using namespace osi3;

void model::CsvOutputReflections::apply(SensorData& sensor_data)
{
    log("Starting .csv output for reflections");

    if (sensor_data.sensor_view_size() == 0)
    {
        log("No sensor view received");
        return;
    }

    if (!sensor_data.sensor_view(0).has_global_ground_truth())
    {
        log("No global ground truth received");
        return;
    }

    auto time_nanos = sensor_data.sensor_view(0).global_ground_truth().timestamp().nanos();
    auto time_seconds = sensor_data.sensor_view(0).global_ground_truth().timestamp().seconds();
    double timestamp = (float)time_seconds + (double)time_nanos / 1000000000;

    /// Loop over all received sensor views
    for (const auto& sensor_view : sensor_data.sensor_view())
    {

        auto no_of_radar_sensors = sensor_view.radar_sensor_view_size();
        if (profile.sensor_view_configuration.radar_sensor_view_configuration_size() != no_of_radar_sensors)
        {
            alert("Number of RadarSensorViews different to profile/SensorViewConfiguration!");
        }

        /// Loop over all received lidar sensors per sensor view
        if (no_of_radar_sensors > 0)
        {
            if (first_call)
            {
#include <csvoutputreflections/set_csv_file_path_reflections.cpp>
                write_first_line_to_csv(file_path_reflections);
                first_call = false;
            }
            for (const auto& radar_sensor_view : sensor_view.radar_sensor_view())
            {
                /// Run through all rendering results
                size_t radar_reflection_idx = 0;
                for (const auto& reflection : radar_sensor_view.reflection())
                {
                    auto elevation_rad = reflection.source_vertical_angle();
                    auto azimuth_rad = reflection.source_horizontal_angle();
                    auto distance = 0.5 * reflection.time_of_flight() * speed_of_light;
                    auto signal_strength_in_dBm = reflection.signal_strength();

                    write_data_to_csv(
                        file_path_reflections, timestamp, radar_reflection_idx, azimuth_rad * 180 / M_PI, elevation_rad * 180 / M_PI, distance, signal_strength_in_dBm);
                    radar_reflection_idx++;
                }
            }
        }
        else
        {
            log("Lidar or radar sensor view empty at timestamp " + std::to_string(timestamp));
        }
    }
}

void CsvOutputReflections::write_first_line_to_csv(const std::string& path)
{
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << "timestamp_in_s, reflection_id, azimuth_in_deg, elevation_in_deg, distance_in_m, signal_strength_in_dBm" << std::endl;
    my_file.close();
}

void CsvOutputReflections::write_data_to_csv(const std::string& path,
                                             double timestamp,
                                             size_t reflection_idx,
                                             double azimuth_in_deg,
                                             double elevation_in_deg,
                                             double distance,
                                             double signal_strength)
{
    std::fstream my_file;
    my_file.open(path, std::ios::app);
    my_file << timestamp << ", " << reflection_idx << ", " << azimuth_in_deg << ", " << elevation_in_deg << ", " << distance << ", " << signal_strength << std::endl;
    my_file.close();
}
