//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

//#pragma once

#include "detectionsensing/cfar_functions.hpp"

#include <iostream>

#include "detectionsensing/DetectionSensing.hpp"

/**
 * Main control function. Performs peak detection and subbin Interpolation for given cell under test
 *
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @return
 *      preCluster struct
 */

RawDetection get_detections_by_spectral_interpolation(radar_cuboid& cube, int range, int doppler, int azimuth, int elevation, const model::Profile& profile)
{
    RawDetection raw_detection = {NAN, NAN, NAN, NAN};
    float cut_power = get_cube_power_at(cube, range, doppler, azimuth, elevation);
    if (check_local_max(cube,
                        range,
                        doppler,
                        azimuth,
                        elevation,
                        cut_power,
                        profile.detection_sensing_parameters.number_range_bin,
                        profile.detection_sensing_parameters.number_doppler_bin,
                        profile.detection_sensing_parameters.number_azimuth_bin,
                        profile.detection_sensing_parameters.number_elevation_bin))
    {
        if (get_cfar_threshold_beam(cube,
                                    range,
                                    doppler,
                                    azimuth,
                                    elevation,
                                    profile.detection_sensing_parameters.number_azimuth_bin,
                                    profile.detection_sensing_parameters.azimuth_CFAR_window_size,
                                    profile.detection_sensing_parameters.azimuth_CFAR_ref_pos,
                                    profile.detection_sensing_parameters.noise_value,
                                    profile.detection_sensing_parameters.azimuth_CFAR_factor) < cut_power)
        {
            if (get_cfar_threshold_range(cube,
                                         range,
                                         doppler,
                                         azimuth,
                                         elevation,
                                         profile.detection_sensing_parameters.range_CFAR_window_size,
                                         profile.detection_sensing_parameters.number_range_bin,
                                         profile.detection_sensing_parameters.noise_value,
                                         profile.detection_sensing_parameters.range_CFAR_ref_pos,
                                         profile.detection_sensing_parameters.range_CFAR_factor) < cut_power)
            {
                if (get_cfar_threshold_elevation(cube,
                                                 range,
                                                 doppler,
                                                 azimuth,
                                                 elevation,
                                                 profile.detection_sensing_parameters.number_elevation_bin,
                                                 profile.detection_sensing_parameters.elevation_CFAR_window_size,
                                                 profile.detection_sensing_parameters.elevation_CFAR_ref_pos,
                                                 profile.detection_sensing_parameters.noise_value,
                                                 profile.detection_sensing_parameters.elevation_CFAR_factor) < cut_power)
                {

                    // if test cell is local maximum and declared as a detection by the OS-CFAR algorithms in azimuth and
                    // range dimension the spectral estimation is performed and the results are stored
                    raw_detection.sub_bin_range =
                        correct_windowing_bias_interpolation_range(cube,
                                                                   range,
                                                                   doppler,
                                                                   azimuth,
                                                                   elevation,
                                                                   profile.detection_sensing_parameters.number_range_bin,
                                                                   static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_range_interpolation),
                                                                   profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                   profile.detection_sensing_parameters.lookup_stepsize_interpolation) +
                        static_cast<float>(range);
                    raw_detection.sub_bin_doppler =
                        correct_windowing_bias_interpolation_doppler(cube,
                                                                     range,
                                                                     doppler,
                                                                     azimuth,
                                                                     elevation,
                                                                     profile.detection_sensing_parameters.number_doppler_bin,
                                                                     static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_doppler_interpolation),
                                                                     profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                     profile.detection_sensing_parameters.lookup_stepsize_interpolation) +
                        static_cast<float>(doppler);
                    raw_detection.sub_bin_azimuth =
                        correct_windowing_bias_interpolation_beam(cube,
                                                                  range,
                                                                  doppler,
                                                                  azimuth,
                                                                  elevation,
                                                                  profile.detection_sensing_parameters.number_azimuth_bin,
                                                                  static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_azimuth_interpolation),
                                                                  profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                  profile.detection_sensing_parameters.lookup_stepsize_interpolation) +
                        static_cast<float>(azimuth);
                    raw_detection.sub_bin_elevation =
                        correct_windowing_bias_interpolation_elevation(cube,
                                                                       range,
                                                                       doppler,
                                                                       azimuth,
                                                                       elevation,
                                                                       profile.detection_sensing_parameters.number_elevation_bin,
                                                                       static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_elevation_interpolation),
                                                                       profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                       profile.detection_sensing_parameters.lookup_stepsize_interpolation) +
                        static_cast<float>(elevation);

                    raw_detection.raw_detection_power = get_cube_power_at(cube, range, doppler, azimuth, elevation);
                }
            }
        }
    }
    return raw_detection;
}

/**
 * Check whether the given cell is a local maximum (peak) in all three dimensions
 *
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in all dimensions
 */

bool check_local_max(radar_cuboid& cube,
                     int range,
                     int doppler,
                     int azimuth,
                     int elevation,
                     float cut_power,
                     int number_range_bin,
                     int number_doppler_bin,
                     int number_azimuth_bin,
                     int number_elevation_bin)
{
    bool r_max;
    bool d_max;
    bool b_max;
    bool e_max;
    r_max = check_local_max_range(cube, range, doppler, azimuth, elevation, cut_power, number_range_bin);
    d_max = check_local_max_doppler(cube, range, doppler, azimuth, elevation, cut_power, number_doppler_bin);
    b_max = check_local_max_beam(cube, range, doppler, azimuth, elevation, cut_power, number_azimuth_bin);
    e_max = check_local_max_elevation(cube, range, doppler, azimuth, elevation, cut_power, number_elevation_bin);
    return r_max && d_max && b_max && e_max;
}

/**
 * Check whether the given cell is a local maximum (peak) in range dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. If this is the case, only the direct neigbour-cell is considered in the check
 * for maximum.
 *
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in range dimension
 */
bool check_local_max_range(radar_cuboid& cube, int range, int doppler, int azimuth, int elevation, float cut_power, int number_range_bin)
{
    bool output;
    if (range == 0)
    {
        output = (cut_power >= get_cube_power_at(cube, range + 1, doppler, azimuth, elevation));
    }
    else if (range == -1)
    {
        output = (cut_power >= get_cube_power_at(cube, range - 1, doppler, azimuth, elevation));
    }
    else if (range == number_range_bin - 1)
    {
        output = (cut_power >= get_cube_power_at(cube, range - 1, doppler, azimuth, elevation));
    }
    else
    {
        output = ((get_cube_power_at(cube, range - 1, doppler, azimuth, elevation) < cut_power) && (cut_power > get_cube_power_at(cube, range + 1, doppler, azimuth, elevation)));
    }
    return output;
}

/**
 * Check whether the given cell is a local maximum (peak) in doppler dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 *
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in doppler dimension
 */
bool check_local_max_doppler(radar_cuboid& cube, int range, int doppler, int azimuth, int elevation, float cut_power, int number_doppler_bin)
{
    bool output;
    if (doppler == number_doppler_bin - 1)
    {
        output = ((get_cube_power_at(cube, range, 0, azimuth, elevation) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler - 1, azimuth, elevation)));
    }
    else if (doppler == 0)
    {
        output = ((get_cube_power_at(cube, range, number_doppler_bin - 1, azimuth, elevation) < cut_power) &&
                  (cut_power > get_cube_power_at(cube, range, doppler + 1, azimuth, elevation)));
    }
    else
    {
        output = ((get_cube_power_at(cube, range, doppler - 1, azimuth, elevation) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler + 1, azimuth, elevation)));
    }
    return output;
}

/**
 * Check whether the given cell is a local maximum (peak) in beam dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 *
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in beam dimension
 */
bool check_local_max_beam(radar_cuboid& cube, int range, int doppler, int azimuth, int elevation, float cut_power, int number_azimuth_bin)
{
    bool output;
    if (azimuth == number_azimuth_bin - 1)
    {
        output = ((get_cube_power_at(cube, range, doppler, 0, elevation) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler, azimuth - 1, elevation)));
    }
    else if (azimuth == 0)
    {
        output = ((get_cube_power_at(cube, range, doppler, number_azimuth_bin - 1, elevation) < cut_power) &&
                  (cut_power > get_cube_power_at(cube, range, doppler, azimuth + 1, elevation)));
    }
    else
    {
        output = ((get_cube_power_at(cube, range, doppler, azimuth - 1, elevation) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler, azimuth + 1, elevation)));
    }
    return output;
}

/*!
 * Check whether the given cell is a local maximum (peak) in elevation dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 * @param cube      The original radar cube
 * @param range         Range position
 * @param doppler         Doppler position
 * @param azimuth         Beam position
 * @param elevation         Elevation position
 * @param cut_power Power value of cell under test
 * @param number_elevation_bin
 * @return
 *      True if local max (peak) in beam dimension
 */
bool check_local_max_elevation(radar_cuboid& cube, int range, int doppler, int azimuth, int elevation, float cut_power, int number_elevation_bin)
{
    bool output;
    if (elevation == number_elevation_bin - 1)
    {
        output = ((get_cube_power_at(cube, range, doppler, azimuth, 0) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler, azimuth, elevation - 1)));
    }
    else if (elevation == 0)
    {
        output = ((get_cube_power_at(cube, range, doppler, azimuth, number_elevation_bin - 1) < cut_power) &&
                  (cut_power > get_cube_power_at(cube, range, doppler, azimuth, elevation + 1)));
    }
    else
    {
        output = ((get_cube_power_at(cube, range, doppler, azimuth, elevation - 1) < cut_power) && (cut_power > get_cube_power_at(cube, range, doppler, azimuth, elevation + 1)));
    }
    return output;
}

/**
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in range dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Threshold value for given cell
 */
float get_cfar_threshold_range(radar_cuboid& cube,
                               int range,
                               int doppler,
                               int azimuth,
                               int elevation,
                               int range_cfar_window_size,
                               int number_range_bin,
                               float noise_value,
                               int range_cfar_ref_pos,
                               float range_cfar_factor)
{
    auto half_window = static_cast<int8_t>(range_cfar_window_size / 2);
    std::vector<float> cfar_window;
    float ref_value;
    // Set the reference window and sort it. Thereby check if a symmetric positioning of the reference window around
    // the CUT is possible. If this is not the case, a non-symmetric reference window is used:
    if (range < half_window)
    {
        for (int cnt = 0; cnt < range_cfar_window_size + 1; ++cnt)
        {
            if (range != cnt)
            {
                cfar_window.push_back(get_cube_power_at(cube, cnt, doppler, azimuth, elevation));
            }
        }
    }
    else if (range + half_window > number_range_bin - 1)
    {
        for (int cnt = number_range_bin - 1 - range_cfar_window_size; cnt < number_range_bin; ++cnt)
        {
            if (range != cnt)
            {
                cfar_window.push_back(get_cube_power_at(cube, cnt, doppler, azimuth, elevation));
            }
        }
    }
    else
    {
        for (int cnt = range - half_window; cnt < range + half_window + 1; ++cnt)
        {
            if (range != cnt)
            {
                cfar_window.push_back(get_cube_power_at(cube, cnt, doppler, azimuth, elevation));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(cfar_window.begin(), cfar_window.end());
    ref_value = cfar_window[range_cfar_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0)
    {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1 / noise_value);
        ref_value = static_cast<float>(distribution(generator));
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * range_cfar_factor;
}

/**
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in beam dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Threshold value for given cell
 */
float get_cfar_threshold_beam(radar_cuboid& cube,
                              int range,
                              int doppler,
                              int azimuth,
                              int elevation,
                              int number_azimuth_bin,
                              int azimuth_cfar_window_size,
                              int azimuth_cfar_ref_pos,
                              float noise_value,
                              float azimuth_cfar_factor)
{
    auto half_window = static_cast<int>(azimuth_cfar_window_size / 2);
    std::vector<float> cfar_window;
    int window_size_buffer;
    float ref_value;
    // Declaration of the reference-window. If no symmetric window around the test-cell is possible, the reference
    // window is expanded on the other side of the spectrum, allowing a symmetric reference window in all cases:
    if (azimuth < half_window)
    {
        for (int cnt = 0; cnt < azimuth + half_window + 1; ++cnt)
        {
            if (cnt != azimuth)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, cnt, elevation));
            }
        }
        window_size_buffer = static_cast<int>(cfar_window.size());
        for (int cnt = number_azimuth_bin - azimuth_cfar_window_size + window_size_buffer; cnt < number_azimuth_bin; ++cnt)
        {
            cfar_window.push_back(get_cube_power_at(cube, range, doppler, cnt, elevation));
        }
    }
    else if (azimuth + half_window > number_azimuth_bin - 1)
    {
        for (int cnt = azimuth - half_window; cnt < number_azimuth_bin; ++cnt)
        {
            if (cnt != azimuth)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, cnt, elevation));
            }
        }
        window_size_buffer = static_cast<int>(cfar_window.size());
        for (int cnt = 0; cnt < azimuth_cfar_window_size - window_size_buffer; ++cnt)
        {
            cfar_window.push_back(get_cube_power_at(cube, range, doppler, cnt, elevation));
        }
    }
    else
    {
        for (int cnt = azimuth - half_window; cnt < azimuth + half_window + 1; ++cnt)
        {
            if (cnt != azimuth)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, cnt, elevation));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(cfar_window.begin(), cfar_window.end());
    ref_value = cfar_window[azimuth_cfar_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0)
    {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1 / noise_value);
        ref_value = static_cast<float>(distribution(generator));
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * azimuth_cfar_factor;
}

/*!
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in elevation dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @param number_elevation_bin
 * @param elevation_cfar_window_size
 * @param elevation_cfar_ref_pos
 * @param noise_value
 * @param elevation_cfar_factor
 * @return
 */
float get_cfar_threshold_elevation(radar_cuboid& cube,
                                   int range,
                                   int doppler,
                                   int azimuth,
                                   int elevation,
                                   int number_elevation_bin,
                                   int elevation_cfar_window_size,
                                   int elevation_cfar_ref_pos,
                                   float noise_value,
                                   float elevation_cfar_factor)
{
    int half_window = elevation_cfar_window_size / 2;
    std::vector<float> cfar_window;
    int window_size_buffer;
    float ref_value;
    // Declaration of the reference-window. If no symmetric window around the test-cell is possible, the reference
    // window is expanded on the other side of the spectrum, allowing a symmetric reference window in all cases:
    if (elevation < half_window)
    {
        for (int cnt = 0; cnt < elevation + half_window + 1; ++cnt)
        {
            if (cnt != elevation)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, azimuth, cnt));
                std::cout << "e < half_window, cnt: " + std::to_string(cnt) << std::endl;
            }
        }
        window_size_buffer = static_cast<int>(cfar_window.size());
        for (int cnt = number_elevation_bin - elevation_cfar_window_size + window_size_buffer; cnt < number_elevation_bin; ++cnt)
        {
            cfar_window.push_back(get_cube_power_at(cube, range, doppler, azimuth, cnt));
        }
    }
    else if (elevation + half_window > number_elevation_bin - 1)
    {
        for (int cnt = elevation - half_window; cnt < number_elevation_bin; ++cnt)
        {
            if (cnt != elevation)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, azimuth, cnt));
            }
        }
        window_size_buffer = static_cast<int>(cfar_window.size());
        for (int cnt = 0; cnt < elevation_cfar_window_size - window_size_buffer; ++cnt)
        {
            cfar_window.push_back(get_cube_power_at(cube, range, doppler, azimuth, cnt));
        }
    }
    else
    {
        for (int cnt = elevation - half_window; cnt < elevation + half_window + 1; ++cnt)
        {
            if (cnt != elevation)
            {
                cfar_window.push_back(get_cube_power_at(cube, range, doppler, azimuth, cnt));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(cfar_window.begin(), cfar_window.end());
    ref_value = cfar_window[elevation_cfar_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0)
    {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1 / noise_value);
        ref_value = static_cast<float>(distribution(generator));
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * elevation_cfar_factor;
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in range dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * range dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Decimal part of the peak position in range dimension
 */
float correct_windowing_bias_interpolation_range(radar_cuboid& cube,
                                                 int range,
                                                 int doppler,
                                                 int azimuth,
                                                 int elevation,
                                                 int number_range_bin,
                                                 std::vector<float> lookup_range_interpolation,
                                                 int lookup_size_interpolation,
                                                 float lookup_stepsize_interpolation)
{
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t is_next;
    // Get maximum neighbour of the CUT in range-dimension, thereby check if CUT is at beginning or end of spectrum.
    // If at beginning or end only direct neighbour is considered.
    if (range == 0)
    {
        max_neighbor_pow = get_cube_power_at(cube, range + 1, doppler, azimuth, elevation);
        is_next = 1;
    }
    else if (range == number_range_bin - 1)
    {
        max_neighbor_pow = get_cube_power_at(cube, range - 1, doppler, azimuth, elevation);
        is_next = -1;
    }
    else
    {
        next_pow = get_cube_power_at(cube, range + 1, doppler, azimuth, elevation);
        prev_pow = get_cube_power_at(cube, range - 1, doppler, azimuth, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, range, doppler, azimuth, elevation) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (!check_local_max_range(cube, range, doppler, azimuth, elevation, get_cube_power_at(cube, range, doppler, azimuth, elevation), number_range_bin))
    {
        sub_bin_correction = 0;
    }
    else if (ref_value < lookup_range_interpolation[lookup_size_interpolation - 1])
    {
        sub_bin_correction = 0.5;
    }
    else if (ref_value > lookup_range_interpolation[0])
    {
        sub_bin_correction = 0;
    }
    else
    {
        // search the position of the reference value inside the lookup-table
        for (auto& lu_step : lookup_range_interpolation)
        {
            if (ref_value > lu_step)
            {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1) +
                             ((lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt) - lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1)) /
                              (lookup_range_interpolation[lu_step_cnt] - lookup_range_interpolation[lu_step_cnt - 1])) *
                                 (ref_value - lookup_range_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * static_cast<float>(is_next);  // return signed decimal bias estimation
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in doppler dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * doppler dimension. The function corrects the bias induced by the signal-processing of the radar under consideration
 * of the used window functions. The path to the used lookup-table is specified in the parameter files. For details on
 * the creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Decimal part of the peak position in doppler dimension
 */
float correct_windowing_bias_interpolation_doppler(radar_cuboid& cube,
                                                   int range,
                                                   int doppler,
                                                   int azimuth,
                                                   int elevation,
                                                   int number_doppler_bin,
                                                   std::vector<float> lookup_doppler_interpolation,
                                                   int lookup_size_interpolation,
                                                   float lookup_stepsize_interpolation)
{
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t is_next;
    // get maximum neighbour of the CUT in doppler-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (doppler == 0)
    {
        next_pow = get_cube_power_at(cube, range, doppler + 1, azimuth, elevation);
        prev_pow = get_cube_power_at(cube, range, number_doppler_bin - 1, azimuth, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else if (doppler == number_doppler_bin - 1)
    {
        next_pow = get_cube_power_at(cube, range, 0, azimuth, elevation);
        prev_pow = get_cube_power_at(cube, range, doppler - 1, azimuth, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else
    {
        next_pow = get_cube_power_at(cube, range, doppler + 1, azimuth, elevation);
        prev_pow = get_cube_power_at(cube, range, doppler - 1, azimuth, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, range, doppler, azimuth, elevation) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_doppler_interpolation[lookup_size_interpolation - 1])
    {
        sub_bin_correction = 0.5;
    }
    else if (ref_value > lookup_doppler_interpolation[0])
    {
        sub_bin_correction = 0;
    }
    else
    {
        // search the position of the reference value inside the lookup-table
        for (auto& lu_step : lookup_doppler_interpolation)
        {
            if (ref_value > lu_step)
            {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1) +
                             ((lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt) - lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1)) /
                              (lookup_doppler_interpolation[lu_step_cnt] - lookup_doppler_interpolation[lu_step_cnt - 1])) *
                                 (ref_value - lookup_doppler_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * static_cast<float>(is_next);  // return signed decimal bias estimation
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in beam dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * beam dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Decimal part of the peak position in beam dimension
 */
float correct_windowing_bias_interpolation_beam(radar_cuboid& cube,
                                                int range,
                                                int doppler,
                                                int azimuth,
                                                int elevation,
                                                int number_azimuth_bin,
                                                std::vector<float> lookup_azimuth_interpolation,
                                                int lookup_size_interpolation,
                                                float lookup_stepsize_interpolation)
{
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t is_next;
    // get maximum neighbour of the CUT in beam-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (azimuth == 0)
    {
        next_pow = get_cube_power_at(cube, range, doppler, azimuth + 1, elevation);
        prev_pow = get_cube_power_at(cube, range, doppler, number_azimuth_bin - 1, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else if (azimuth == number_azimuth_bin - 1)
    {
        next_pow = get_cube_power_at(cube, range, doppler, 0, elevation);
        prev_pow = get_cube_power_at(cube, range, doppler, azimuth - 1, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else
    {
        next_pow = get_cube_power_at(cube, range, doppler, azimuth + 1, elevation);
        prev_pow = get_cube_power_at(cube, range, doppler, azimuth - 1, elevation);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, range, doppler, azimuth, elevation) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_azimuth_interpolation[lookup_size_interpolation - 1])
    {
        sub_bin_correction = 0.5;
    }
    else if (ref_value > lookup_azimuth_interpolation[0])
    {
        sub_bin_correction = 0;
    }
    else
    {
        // search the position of the reference value inside the lookup-table
        for (auto& lu_step : lookup_azimuth_interpolation)
        {
            if (ref_value > lu_step)
            {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1) +
                             ((lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt) - lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1)) /
                              (lookup_azimuth_interpolation[lu_step_cnt] - lookup_azimuth_interpolation[lu_step_cnt - 1])) *
                                 (ref_value - lookup_azimuth_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * static_cast<float>(is_next);  // return signed decimal bias estimation
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in elevation dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * beam dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param range         Range peak position
 * @param doppler         Doppler peak position
 * @param azimuth         Beam peak position
 * @param elevation         Elevation peak position
 * @return
 *      Decimal part of the peak position in beam dimension
 */
float correct_windowing_bias_interpolation_elevation(radar_cuboid& cube,
                                                     int range,
                                                     int doppler,
                                                     int azimuth,
                                                     int elevation,
                                                     int number_elevation_bin,
                                                     std::vector<float> lookup_elevation_interpolation,
                                                     int lookup_size_interpolation,
                                                     float lookup_stepsize_interpolation)
{
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t is_next;
    // get maximum neighbour of the CUT in beam-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (elevation == 0)
    {
        next_pow = get_cube_power_at(cube, range, doppler, azimuth, elevation + 1);
        prev_pow = get_cube_power_at(cube, range, doppler, azimuth, number_elevation_bin - 1);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else if (elevation == number_elevation_bin - 1)
    {
        next_pow = get_cube_power_at(cube, range, doppler, azimuth, 0);
        prev_pow = get_cube_power_at(cube, range, doppler, azimuth, elevation - 1);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    else
    {
        next_pow = get_cube_power_at(cube, range, doppler, azimuth, elevation + 1);
        prev_pow = get_cube_power_at(cube, range, doppler, azimuth, elevation - 1);
        if (prev_pow > next_pow)
        {
            max_neighbor_pow = prev_pow;
            is_next = -1;
        }
        else
        {
            max_neighbor_pow = next_pow;
            is_next = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, range, doppler, azimuth, elevation) / max_neighbor_pow;

    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_elevation_interpolation[lookup_size_interpolation - 1])
    {
        sub_bin_correction = 0.5;
    }
    else if (ref_value > lookup_elevation_interpolation[0])
    {
        sub_bin_correction = 0;
    }
    else
    {
        // search the position of the reference value inside the lookup-table
        for (auto& lu_step : lookup_elevation_interpolation)
        {
            if (ref_value > lu_step)
            {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1) +
                             ((lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt) - lookup_stepsize_interpolation * static_cast<float>(lu_step_cnt - 1)) /
                              (lookup_elevation_interpolation[lu_step_cnt] - lookup_elevation_interpolation[lu_step_cnt - 1])) *
                                 (ref_value - lookup_elevation_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * static_cast<float>(is_next);  // return signed decimal bias estimation
}

/**
 * Accesses the power at the (range, doppler, beam) position.
 *
 * @param cube      The original radar cube
 * @param range     Range position
 * @param doppler   Doppler position
 * @param beam      Beam position
 * @param elevation Elevation position
 * @return
 *      Power at the (range, doppler, beam) position.
 */
float get_cube_power_at(radar_cuboid& cube, const int range, const int doppler, const int beam, const int elevation)
{
    return static_cast<float>(pow(cube[range][doppler][beam][elevation][0], 2) + pow(cube[range][doppler][beam][elevation][1], 2));
}
