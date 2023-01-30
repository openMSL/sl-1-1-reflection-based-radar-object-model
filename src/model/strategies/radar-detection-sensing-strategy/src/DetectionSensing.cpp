//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef Speed_of_Light
#define Speed_of_Light 299792458
#endif

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "detectionsensing/DetectionSensing.hpp"
#include "detectionsensing/cfar_functions.hpp"
#include <random>
#include <string>
#include <vector>
#include <utility>
#include <iostream>


#ifdef _WIN32
#include <math.h>
#else

#include <cmath>

#endif


using namespace model;
using namespace osi3;

void DetectionSensing::apply(SensorData &sensor_data) {
    log("Starting Radar Sensor Model");


    if (sensor_data.sensor_view().size() > 0) {
        auto no_of_radar_sensors = sensor_data.sensor_view(0).radar_sensor_view_size();
        log("Number of simulated radar sensors: " + std::to_string(no_of_radar_sensors));
        if (profile.sensor_view_configuration.radar_sensor_view_configuration_size() != no_of_radar_sensors)
            alert("Number of radar sensor view different to profile/SensorViewConfiguration!");

        /// Get timestamp from ground truth
        double timestamp = (double) sensor_data.sensor_view(0).global_ground_truth().timestamp().seconds() +
                           (double) sensor_data.sensor_view(0).global_ground_truth().timestamp().nanos() / 1000000000;
        log("GT timestamp: " + std::to_string(timestamp));

        if (no_of_radar_sensors > 0) {
            /// Loop over all received Sensors
            for (uint64_t sensor_idx = 0; sensor_idx < no_of_radar_sensors; sensor_idx++) {
                /// Pointer to radar sensor view for readability
                auto *radar_sensor_view = &sensor_data.sensor_view(0).radar_sensor_view(sensor_idx);

                auto no_of_reflections = radar_sensor_view->reflection_size();
                log("Number of reflections: " + std::to_string(no_of_reflections));

                /// Start new set of detections
                sensor_data.mutable_feature_data()->add_radar_sensor();
                auto *current_sensor = sensor_data.mutable_feature_data()->mutable_radar_sensor(sensor_idx);
                current_sensor->clear_detection();
                current_sensor->mutable_header()->mutable_mounting_position()->CopyFrom(
                        profile.sensor_view_configuration.radar_sensor_view_configuration(
                                sensor_idx).mounting_position());
                /// Get configuration of current sensor
                auto number_range_bin = profile.detection_sensing_parameters.number_range_bin;
                auto number_doppler_bin = profile.detection_sensing_parameters.number_doppler_bin;
                auto number_azimuth_bin = profile.detection_sensing_parameters.number_azimuth_bin;
                auto number_elevation_bin = profile.detection_sensing_parameters.number_elevation_bin;
                auto range_resolution = profile.detection_sensing_parameters.range_resolution;
                auto doppler_resolution = profile.detection_sensing_parameters.doppler_resolution;
                auto azimuth_resolution = profile.detection_sensing_parameters.azimuth_resolution;
                auto elevation_resolution = profile.detection_sensing_parameters.elevation_resolution;
                auto max_range = profile.sensor_view_configuration.range();
                auto emitter_frequency = profile.sensor_view_configuration.radar_sensor_view_configuration(
                        sensor_idx).emitter_frequency();
                auto emitter_power = profile.detection_sensing_parameters.emitter_strength;
                auto range_window_function = profile.detection_sensing_parameters.range_window_function;
                auto doppler_window_function = profile.detection_sensing_parameters.doppler_window_function;
                auto azimuth_window_function = profile.detection_sensing_parameters.azimuth_window_function;
                auto elevation_window_function = profile.detection_sensing_parameters.elevation_window_function;
                auto window_data_per_bin = profile.detection_sensing_parameters.window_data_per_bin;
                auto bin_affect_range = profile.detection_sensing_parameters.bin_affect_range;
                auto max_number_detections = profile.detection_sensing_parameters.max_number_detections;
                auto power_threshold = profile.detection_sensing_parameters.power_threshold;
                auto num_peaks_max = profile.detection_sensing_parameters.num_peaks_max;
                auto rcs_calibration = profile.detection_sensing_parameters.rcs_calibration;
                /// Initialize additional variables
                std::vector<float> complex_signal_strength(2);
                std::vector<float> complex_signal_strength_sum(2);
                float wavelength = Speed_of_Light / emitter_frequency;

                /// Reserve memory for the radar cuboid to prevent re-allocations at every complex signal strength iteration
                type_radar_cuboid radar_cuboid (number_range_bin, std::vector<std::vector<std::vector<std::vector<float> > > >(number_doppler_bin, std::vector<std::vector<std::vector<float> > >(number_azimuth_bin, std::vector<std::vector<float> >(number_elevation_bin, std::vector<float>(2)))));

                /// Run through all reflections and append valid ones to the radar cuboid
                for (uint64_t reflection_idx = 0; reflection_idx < no_of_reflections; reflection_idx++) {
                    /// Definition of variable for readability
                    float distance = (float)(radar_sensor_view->reflection(reflection_idx).time_of_flight())*static_cast<float>(Speed_of_Light);
                    float signal_strength =
                            pow(10, ((float) (radar_sensor_view->reflection(reflection_idx).signal_strength()) / 10)) *
                            emitter_power;

                    if (distance > 0 && distance <= max_range && signal_strength >= 0 &&
                        std::isfinite(signal_strength)) {
                        /// Definition of variable for readability
                        float relative_velocity = -(float)(radar_sensor_view->reflection(reflection_idx).doppler_shift())*(static_cast<float>(Speed_of_Light))/(2*emitter_frequency);
                        float azimuth_angle = (float) (radar_sensor_view->reflection(
                                reflection_idx).source_horizontal_angle());
                        float elevation_angle = (float) (radar_sensor_view->reflection(
                                reflection_idx).source_vertical_angle());

                        /// Range calculation in bins
                        int range_bin = (int) (distance / range_resolution);
                        float range_bin_decimal = distance / range_resolution - range_bin;
                        int range_window = static_cast<int>(static_cast<float>(size(range_window_function) + 1) / 2 +
                                                            round(range_bin_decimal * window_data_per_bin));

                        /// Radial velocity calculation in bins
                        int doppler_bin = 0;
                        if (relative_velocity > 0) {
                            doppler_bin = static_cast<int>((number_doppler_bin - 1) -
                                                           floor(fmodf(relative_velocity / doppler_resolution,
                                                                       number_doppler_bin)));
                        } else {
                            doppler_bin = static_cast<int>(floor(
                                    fmodf(static_cast<float>(fabs(relative_velocity / doppler_resolution)),
                                          static_cast<float>(number_doppler_bin))));
                        }
                        float doppler_bin_decimal = static_cast<float>(fabs(relative_velocity) / doppler_resolution -
                                                                       floor(fabs(relative_velocity) /
                                                                             doppler_resolution));
                        int doppler_window = static_cast<int>(
                                static_cast<float>(size(doppler_window_function) + 1) / 2 -
                                round(doppler_bin_decimal * window_data_per_bin));

                        /// Azimuth calculation in bins
                        int azimuth_bin = 0;
                        int azimuth_signum;
                        int azimuth_window;
                        float azimuth_bin_decimal = static_cast<float>(sin(azimuth_angle) / azimuth_resolution -
                                                                       floor(sin(azimuth_angle) / azimuth_resolution));
                        if (azimuth_angle > 0) {
                            azimuth_signum = 1;
                            azimuth_bin = static_cast<int>((number_azimuth_bin - 1) - floor(fmodf(
                                    static_cast<float>(fabs(sin(azimuth_angle) / azimuth_resolution)),
                                    (float) number_azimuth_bin)));
                            azimuth_window = static_cast<int>(
                                    static_cast<float>(size(azimuth_window_function) + 1) / 2 -
                                    round(azimuth_bin_decimal * window_data_per_bin));
                        } else {
                            azimuth_signum = -1;
                            azimuth_bin = static_cast<int>(floor(
                                    fmodf(static_cast<float>(fabs(sin(azimuth_angle) / azimuth_resolution)),
                                          number_azimuth_bin)));
                            azimuth_window = static_cast<int>(
                                    static_cast<float>(size(azimuth_window_function) + 1) / 2 +
                                    round(azimuth_bin_decimal * window_data_per_bin));
                        }

                        /// Elevation calculation in bins
                        int elevation_bin = 0;
                        int elevation_signum;
                        int elevation_window;
                        float elevation_bin_decimal = static_cast<float>(sin(elevation_angle) / elevation_resolution -
                                                                         floor(sin(elevation_angle) / elevation_resolution));
                        if (elevation_angle > 0) {
                            elevation_signum = 1;
                            elevation_bin = static_cast<int>((number_elevation_bin - 1) - floor(fmodf(
                                    static_cast<float>(fabs(sin(elevation_angle) / elevation_resolution)),
                                    (float) number_elevation_bin)));
                            elevation_window = static_cast<int>(
                                    static_cast<float>(size(elevation_window_function) + 1) / 2 -
                                    round(elevation_bin_decimal * window_data_per_bin));
                        } else {
                            elevation_signum = -1;
                            elevation_bin = static_cast<int>(floor(
                                    fmodf(static_cast<float>(fabs(sin(elevation_angle) / elevation_resolution)),
                                          number_elevation_bin)));
                            elevation_window = static_cast<int>(
                                    static_cast<float>(size(elevation_window_function) + 1) / 2 +
                                    round(elevation_bin_decimal * window_data_per_bin));
                        }

                        /// Calculation of phase based on travelled distance of ray
                        auto phase = static_cast<float>(std::fmod(2 * distance, wavelength) / wavelength * 2 * M_PI);

                        /// Calculation of complex signal strength based on received signal and phase
                        complex_signal_strength[0] = static_cast<float>(sqrt(signal_strength) *
                                                                        cos(phase));    // real part
                        complex_signal_strength[1] = static_cast<float>(sqrt(signal_strength) *
                                                                        sin(phase));    // imag part

                        /// Calculation of radar cuboid with windowing effects
                        for (int range_index = static_cast<int>(-bin_affect_range) + 1;
                             range_index <= static_cast<int>(bin_affect_range) - 1; range_index++) {
                            /// Range FFT (convolve window function over range bin)
                            if (range_bin + range_index > 0 && range_bin + range_index < number_range_bin) {
                                for (int doppler_index = -static_cast<int>(bin_affect_range) + 1;
                                     doppler_index <= static_cast<int>(bin_affect_range) - 1; doppler_index++) {
                                    /// Doppler FFT (convolve window function over doppler bin)
                                    int doppler_convolve_index = 0;
                                    if (doppler_bin + doppler_index < 0) {
                                        doppler_convolve_index = 1;
                                    } else if (doppler_bin + doppler_index >= number_doppler_bin) {
                                        doppler_convolve_index = -1;
                                    }

                                    for (int azimuth_index = -number_azimuth_bin / 2 + 1; azimuth_index <= number_azimuth_bin / 2; azimuth_index++) {
                                        /// Azimuth FFT (convolve window function over azimuth bin)
                                        int azimuth_convolve_index = 0;
                                        if (azimuth_bin + azimuth_index < 0) {
                                            azimuth_convolve_index = 1;
                                        } else if (azimuth_bin + azimuth_index >= number_azimuth_bin) {
                                            azimuth_convolve_index = -1;
                                        }
                                        for (int elevation_index = -number_elevation_bin / 2 + 1; elevation_index <= number_elevation_bin / 2; elevation_index++) {
                                            /// Elevation FFT (convolve window function over elevation bin)
                                            int elevation_convolve_index = 0;
                                            if (elevation_bin + elevation_index < 0) {
                                                elevation_convolve_index = 1;
                                            } else if (elevation_bin + elevation_index >= number_elevation_bin) {
                                                elevation_convolve_index = -1;
                                            }
                                            /// Calculation of windowing factor as a FFT equivalent representation
                                            float windowing_factor =
                                                    range_window_function[(range_window -
                                                                          range_index * window_data_per_bin)%(size(range_window_function))] *
                                                    doppler_window_function[(doppler_window +
                                                                            doppler_index * window_data_per_bin)%(size(doppler_window_function))] *
                                                    azimuth_window_function[(azimuth_window +
                                                                            azimuth_signum * azimuth_index *
                                                                            window_data_per_bin)%(size(azimuth_window_function)-1)] *
                                                    elevation_window_function[(elevation_window +
                                                                              elevation_signum * elevation_index *
                                                                              window_data_per_bin)%(size(elevation_window_function)-1)];

                                            /// Calculation of complex signal strength
                                            complex_signal_strength_sum[0] = radar_cuboid[range_bin + range_index]
                                                                             [doppler_bin + doppler_index +
                                                                              doppler_convolve_index * number_doppler_bin]
                                                                             [azimuth_bin + azimuth_index +
                                                                              azimuth_convolve_index *
                                                                              number_azimuth_bin]
                                                                              [elevation_bin + elevation_index +
                                                                              elevation_convolve_index * number_elevation_bin][0] +
                                                                             complex_signal_strength[0] * windowing_factor;
                                            complex_signal_strength_sum[1] = radar_cuboid[range_bin + range_index]
                                                                             [doppler_bin + doppler_index +
                                                                              doppler_convolve_index * number_doppler_bin]
                                                                             [azimuth_bin + azimuth_index +
                                                                              azimuth_convolve_index *
                                                                              number_azimuth_bin]
                                                                             [elevation_bin + elevation_index +
                                                                              elevation_convolve_index * number_elevation_bin][1] +
                                                                             complex_signal_strength[1] * windowing_factor;
                                            radar_cuboid[range_bin + range_index][doppler_bin + doppler_index +
                                                                                  doppler_convolve_index *
                                                                                  number_doppler_bin]
                                                                                  [azimuth_bin + azimuth_index + azimuth_convolve_index *
                                                                                  number_azimuth_bin]
                                                                                  [elevation_bin + elevation_index +
                                                                                  elevation_convolve_index * number_elevation_bin][0] =
                                                                                          complex_signal_strength_sum[0];
                                            radar_cuboid[range_bin + range_index][doppler_bin + doppler_index +
                                                                                  doppler_convolve_index *
                                                                                  number_doppler_bin]
                                                                                  [azimuth_bin + azimuth_index +
                                                                                  azimuth_convolve_index * number_azimuth_bin]
                                                                                  [elevation_bin + elevation_index +
                                                                                  elevation_convolve_index * number_elevation_bin][1] =
                                                                                          complex_signal_strength_sum[1];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                /// CFAR Peakdetection
                for (size_t r = 0; r < number_range_bin; ++r) {
                    float* doppler_vec = new float[number_doppler_bin];
                    memset(doppler_vec, 0, sizeof(doppler_vec));
                    for (size_t d = 0; d < number_doppler_bin; ++d) {
                        for (size_t b = 0; b < number_azimuth_bin; ++b) {
                            for (size_t e = 0; e < number_elevation_bin; ++e) {
                                float power = get_cube_power_at(radar_cuboid, r, d, b, e);
                                /// Look for Beam bin that contains the maximum power in each Doppler bin
                                if (power > power_threshold && power > doppler_vec[d]) {
                                    doppler_vec[d] = power;
                                }
                            }
                        }
                    }
                    std::vector<int> doppler_idx;
                    for (size_t p = 0; p < num_peaks_max; ++p) {
                        float current_doppler_max = 0;
                        int current_doppler_max_idx = 0;
                        for (size_t d = 0; d < number_doppler_bin; ++d) {
                            if (doppler_vec[d] > 0 && doppler_vec[d] > current_doppler_max) {
                                current_doppler_max = doppler_vec[d];
                                current_doppler_max_idx = d;
                            }
                        }
                        if (current_doppler_max == 0) {
                            p = num_peaks_max;
                            continue;
                        }
                        if (current_doppler_max_idx == 0) {
                            doppler_vec[number_doppler_bin - 1] = 0;
                            doppler_vec[current_doppler_max_idx] = 0;
                            doppler_vec[current_doppler_max_idx + 1] = 0;
                        } else if (current_doppler_max_idx == number_doppler_bin - 1) {
                            doppler_vec[current_doppler_max_idx - 1] = 0;
                            doppler_vec[current_doppler_max_idx] = 0;
                            doppler_vec[0] = 0;
                        } else {
                            doppler_vec[current_doppler_max_idx - 1] = 0;
                            doppler_vec[current_doppler_max_idx] = 0;
                            doppler_vec[current_doppler_max_idx + 1] = 0;
                        }
                        doppler_idx.push_back(current_doppler_max_idx);
                    }
                    for (int dpl: doppler_idx) {
                        for (size_t azi = 0; azi < number_azimuth_bin; ++azi) {
                            for(size_t ele = 0; ele < number_elevation_bin; ++ele) {
                                raw_detection rawDetection = get_detections_by_spectral_interpolation(radar_cuboid, r, dpl,
                                                                                                      azi, ele, profile,
                                                                                                      sensor_idx);

                                if (!std::isnan(rawDetection.sub_bin_range)) {
                                    /// Add new detection and fill it
                                    auto detection = current_sensor->add_detection();
                                    detection->mutable_position()->set_distance(
                                            rawDetection.sub_bin_range * range_resolution);
                                    if (static_cast<float>(rawDetection.sub_bin_azimuth) /
                                        static_cast<float>(number_azimuth_bin) < 0.5) {
                                        detection->mutable_position()->set_azimuth(
                                                rawDetection.sub_bin_azimuth * azimuth_resolution);
                                    } else {
                                        detection->mutable_position()->set_azimuth(
                                                (static_cast<float>(number_azimuth_bin) - rawDetection.sub_bin_azimuth) *
                                                azimuth_resolution * -1);
                                    }
                                    if (static_cast<float>(rawDetection.sub_bin_elevation) /
                                        static_cast<float>(number_elevation_bin) < 0.5) {
                                        detection->mutable_position()->set_elevation(
                                                rawDetection.sub_bin_elevation * elevation_resolution);
                                    } else {
                                        detection->mutable_position()->set_elevation(
                                                (static_cast<float>(number_elevation_bin) - rawDetection.sub_bin_elevation) *
                                                        elevation_resolution * -1);
                                    }
                                    detection->set_radial_velocity(rawDetection.sub_bin_doppler * doppler_resolution);
                                    detection->set_rcs(10 * log10(rawDetection.raw_detection_power *
                                                                  pow(detection->mutable_position()->distance(), 4)) +
                                                       rcs_calibration);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
