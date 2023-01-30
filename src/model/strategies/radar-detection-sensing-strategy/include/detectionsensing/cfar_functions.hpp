//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef REFLECTIONBASEDRADARMODEL_CFAR_FUNCTIONS_HPP
#define REFLECTIONBASEDRADARMODEL_CFAR_FUNCTIONS_HPP

/** --------------------------------------------------------------------------------------------------------------------
 *  @file    cfar_functions.hpp
 *  @author  Markus Klingenstein
 *  @date    12.10.2019
 *  @version 0.1
 *
 * This header contains the functions to perform target detection and spectral estimation of the detection position.
 * The main function get_detections_by_spectral_interpolation controls the the algorithm. It performs the classification of a specified test-cell
 * inside the given radar cube. First, a check for local maximum of the test-cell is performed. Then the classification
 * is done by a two staged ordered statistics constant false alarm rate (OS-CFAR) algorithm which can be adjusted by the
 * parameters in the parameter header files. Thereby a classification in the beam-dimension is performed, followed by a
 * classification in range-dimension. If the test-cell is classified as a target, a spectral estimation based on
 * predefined lookup-tables is performed in all three dimensions. The return struct contains the interpolated position-
 * values and the power-value inside the test-cell.
 * The lookup-tables are calculated in the separate matlab-file. Also the CFAR-Factor for a given combination of window
 * lenght, reference value position and false alarm rate is calculated in matlab. The following values have shown to
 * have a good performance with the current state of the simulation plugin:
 *
 * const __uint8_t    range_CFAR_window_size = 32;
 * const __uint8_t    range_CFAR_ref_pos = 24;
 * const float        range_CFAR_factor = 16.45;
 * const __uint8_t    beam_CFAR_window_size = 16;
 * const __uint8_t    beam_CFAR_ref_pos = 12;
 * const float        beam_CFAR_factor = 26.97;
 *
 * The reference window size and the reference value position should be considered as ratio: ref_pos / window_size
 * By setting the ratio bigger then 1/2 the false alarms at clutter-edges can be reduced. Setting the ratio to big will
 * result in in a detection loss if several target reflections are close together in the spectrum.
 * The corresponding matlab- and python-files can be found in the repository "math_740_19_radar_cfar". The given python-
 * scripts allow the display and analysis of the output from the CFAR-algorithm in this header.
 *
 * ------------------------------------------------------------------------------------------------------------------ */

#pragma once

#include <vector>
#include <valarray>
#include <random>

    // Data structure containing the results of the CFAR-functions. The values are returned as NAN if no target is
    // detected.
    typedef struct {
        float sub_bin_range;
        float sub_bin_doppler;
        float sub_bin_azimuth;
        float sub_bin_elevation;
        float raw_detection_power;
    } raw_detection;

    using radar_cuboid = std::vector<std::vector<std::vector<std::vector<std::vector<float>>>>>;

    raw_detection get_detections_by_spectral_interpolation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, const model::Profile &profile, uint64_t i);
    //raw_detection get_detections_by_spectral_interpolation(radar_cuboid &cube);

    bool check_local_max(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_range_bin, int number_doppler_bin,
                         int number_azimuth_bin, int number_elevation_bin);
    bool check_local_max_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_range_bin);
    bool check_local_max_doppler(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_doppler_bin);
    bool check_local_max_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_azimuth_bin);
    bool check_local_max_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_elevation_bin);//new
    float get_CFAR_threshold_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int range_CFAR_window_size,
                                   int number_range_bin, float noise_value, int range_CFAR_ref_pos,
                                   float range_CFAR_factor);
    float get_CFAR_threshold_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_azimuth_bin,
                                  int azimuth_CFAR_window_size, int azimuth_CFAR_ref_pos, float noise_value,
                                  float azimuth_CFAR_factor);
    float get_CFAR_threshold_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_elevation_bin,
                                       int elevation_CFAR_window_size, int elevation_CFAR_ref_pos, float noise_value,
                                       float elevation_CFAR_factor);//new
    float
    correct_windowing_bias_interpolation_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_range_bin,
                                               std::vector<float> lookup_range_interpolation, int lookup_size_interpolation,
                                               float lookup_stepsize_interpolation);
    float correct_windowing_bias_interpolation_doppler(radar_cuboid &cube, size_t r, size_t d, size_t b,size_t e,
                                                       int number_doppler_bin,
                                                       std::vector<float> lookup_doppler_interpolation, int lookup_size_interpolation,
                                                       float lookup_stepsize_interpolation);
    float
    correct_windowing_bias_interpolation_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_azimuth_bin,
                                              std::vector<float> lookup_azimuth_interpolation, int lookup_size_interpolation,
                                              float lookup_stepsize_interpolation);
    float
    correct_windowing_bias_interpolation_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_elevation_bin,
                                          std::vector<float> lookup_elevation_interpolation, int lookup_size_interpolation,
                                          float lookup_stepsize_interpolation);//new
    float get_cube_power_at(radar_cuboid &cube, size_t range, size_t doppler, size_t beam, size_t elevation);

#endif //REFLECTIONBASEDRADARMODEL_CFAR_FUNCTIONS_HPP