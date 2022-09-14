//
// Created by lukas on 04/08/2021.
//

//#pragma once

#include "detectionsensing/DetectionSensing.hpp"
#include "detectionsensing/cfar_functions.hpp"

#include <iostream>

/**
 * Main control function. Performs peak detection and subbin Interpolation for given cell under test
 *
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @return
 *      preCluster struct
 */

raw_detection get_detections_by_spectral_interpolation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, const model::Profile &profile, uint64_t sensor_idx) {
    raw_detection RawDetection = {NAN, NAN, NAN, NAN};
    float cut_power = get_cube_power_at(cube, r, d, b, e);
    if (check_local_max(cube, r, d, b, e, cut_power,
                        profile.detection_sensing_parameters.number_range_bin,
                        profile.detection_sensing_parameters.number_doppler_bin,
                        profile.detection_sensing_parameters.number_azimuth_bin,
                        profile.detection_sensing_parameters.number_elevation_bin)) {
        if (get_CFAR_threshold_beam(cube, r, d, b, e,
                                    profile.detection_sensing_parameters.number_azimuth_bin,
                                    profile.detection_sensing_parameters.azimuth_CFAR_window_size,
                                    profile.detection_sensing_parameters.azimuth_CFAR_ref_pos,
                                    profile.detection_sensing_parameters.noise_value,
                                    profile.detection_sensing_parameters.azimuth_CFAR_factor) < cut_power) {
            if (get_CFAR_threshold_range(cube, r, d, b, e,
                                         profile.detection_sensing_parameters.range_CFAR_window_size,
                                         profile.detection_sensing_parameters.number_range_bin,
                                         profile.detection_sensing_parameters.noise_value,
                                         profile.detection_sensing_parameters.range_CFAR_ref_pos,
                                         profile.detection_sensing_parameters.range_CFAR_factor) < cut_power) {
                if(get_CFAR_threshold_elevation(cube, r, d, b, e,
                                                profile.detection_sensing_parameters.number_elevation_bin,
                                                profile.detection_sensing_parameters.elevation_CFAR_window_size,
                                                profile.detection_sensing_parameters.elevation_CFAR_ref_pos,
                                                profile.detection_sensing_parameters.noise_value,
                                                profile.detection_sensing_parameters.elevation_CFAR_factor) < cut_power){

                    // if test cell is local maximum and declared as a detection by the OS-CFAR algorithms in azimuth and
                    // range dimension the spectral estimation is performed and the results are stored
                    RawDetection.sub_bin_range =
                            correct_windowing_bias_interpolation_range(cube, r, d, b, e,
                                                                       profile.detection_sensing_parameters.number_range_bin,
                                                                       static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_range_interpolation),
                                                                       profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                       profile.detection_sensing_parameters.lookup_stepsize_interpolation) + r;
                    RawDetection.sub_bin_doppler =
                            correct_windowing_bias_interpolation_doppler(cube, r, d, b, e,
                                                                         profile.detection_sensing_parameters.number_doppler_bin,
                                                                         static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_doppler_interpolation),
                                                                         profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                         profile.detection_sensing_parameters.lookup_stepsize_interpolation) + d;
                    RawDetection.sub_bin_azimuth =
                            correct_windowing_bias_interpolation_beam(cube, r, d, b, e,
                                                                      profile.detection_sensing_parameters.number_azimuth_bin,
                                                                      static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_azimuth_interpolation),
                                                                      profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                      profile.detection_sensing_parameters.lookup_stepsize_interpolation) + b;
                    RawDetection.sub_bin_elevation =
                            correct_windowing_bias_interpolation_elevation(cube, r, d, b, e,
                                                                           profile.detection_sensing_parameters.number_elevation_bin,
                                                                           static_cast<std::vector<float>>(profile.detection_sensing_parameters.lookup_elevation_interpolation),
                                                                           profile.detection_sensing_parameters.lookup_size_interpolation,
                                                                           profile.detection_sensing_parameters.lookup_stepsize_interpolation) + e;




                    RawDetection.raw_detection_power = get_cube_power_at(cube, r, d, b, e);
                }
            }
        }

    }
    return RawDetection;
}

/**
 * Check whether the given cell is a local maximum (peak) in all three dimensions
 *
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in all dimensions
 */

bool check_local_max(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_range_bin,int number_doppler_bin,
                     int number_azimuth_bin, int number_elevation_bin) {
    bool r_max, d_max, b_max, e_max;
    r_max = check_local_max_range(cube, r, d, b, e, cut_power, number_range_bin);
    d_max = check_local_max_doppler(cube, r, d, b, e, cut_power, number_doppler_bin);
    b_max = check_local_max_beam(cube, r, d, b, e, cut_power, number_azimuth_bin);
    e_max = check_local_max_elevation(cube, r, d, b, e, cut_power, number_elevation_bin);
    return r_max && d_max && b_max && e_max ;
}

/**
 * Check whether the given cell is a local maximum (peak) in range dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. If this is the case, only the direct neigbour-cell is considered in the check
 * for maximum.
 *
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in range dimension
 */
bool check_local_max_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_range_bin) {
    if (r == 0) {
        return (cut_power >= get_cube_power_at(cube, r + 1, d, b, e));
    } else if (r ==  - 1) {
        return (cut_power >= get_cube_power_at(cube, r - 1, d, b, e));
    } else if (r == number_range_bin-1) {
        return (cut_power >= get_cube_power_at(cube, r - 1, d, b, e));
    } else {
        return ((get_cube_power_at(cube, r - 1, d, b, e) < cut_power) && (cut_power >
        get_cube_power_at(cube, r + 1, d, b, e)));
    }
}

/**
 * Check whether the given cell is a local maximum (peak) in doppler dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 *
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in doppler dimension
 */
bool check_local_max_doppler(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_doppler_bin) {
    if (d == number_doppler_bin - 1) {
        return ((get_cube_power_at(cube, r, 0, b, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d - 1, b, e)));
    } else if (d == 0) {
        return ((get_cube_power_at(cube, r, number_doppler_bin - 1, b, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d + 1, b, e)));
    } else {
        return ((get_cube_power_at(cube, r, d - 1, b, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d + 1, b, e)));
    }
}

/**
 * Check whether the given cell is a local maximum (peak) in beam dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 *
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @param cut_power Power value of cell under test
 * @return
 *      True if local max (peak) in beam dimension
 */
bool check_local_max_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_azimuth_bin) {
    if (b == number_azimuth_bin - 1) {
        return ((get_cube_power_at(cube, r, d, 0, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d, b - 1, e)));
    } else if (b == 0) {
        return ((get_cube_power_at(cube, r, d, number_azimuth_bin - 1, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d, b + 1, e)));
    } else {
        return ((get_cube_power_at(cube, r, d, b - 1, e) < cut_power) &&
        (cut_power > get_cube_power_at(cube, r, d, b + 1, e)));
    }
}



/*!
 * Check whether the given cell is a local maximum (peak) in elevation dimension. Thereby a check is performed if the test-
 * cell is at the edge of the spectrum. In that case the power-value at the other side of the spectrum is used for the
 * maximum check.
 * @param cube      The original radar cube
 * @param r         Range position
 * @param d         Doppler position
 * @param b         Beam position
 * @param e         Elevation position
 * @param cut_power Power value of cell under test
 * @param number_elevation_bin
 * @return
 *      True if local max (peak) in beam dimension
 */
bool check_local_max_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, float cut_power, int number_elevation_bin){
    if (e == number_elevation_bin - 1) {
        return ((get_cube_power_at(cube, r, d, b, 0) < cut_power) &&
                (cut_power > get_cube_power_at(cube, r, d, b, e - 1)));
    } else if (e == 0) {
        return ((get_cube_power_at(cube, r, d, b, number_elevation_bin - 1) < cut_power) &&
                (cut_power > get_cube_power_at(cube, r, d, b, e + 1)));
    } else {
        return ((get_cube_power_at(cube, r, d, b, e - 1) < cut_power) &&
                (cut_power > get_cube_power_at(cube, r, d, b, e + 1)));
    }
}

/**
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in range dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Threshold value for given cell
 */
float get_CFAR_threshold_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int range_CFAR_window_size,
                               int number_range_bin, float noise_value, int range_CFAR_ref_pos,
                               float range_CFAR_factor) {
    int8_t half_window = range_CFAR_window_size / 2;
    std::vector<float> CFAR_window;
    float ref_value;
    // Set the reference window and sort it. Thereby check if a symmetric positioning of the reference window around
    // the CUT is possible. If this is not the case, a non symmetric reference window is used:
    if (r < half_window) {
        for (size_t cnt = 0; cnt < range_CFAR_window_size + 1; ++cnt) {
            if (r != cnt) {
                CFAR_window.push_back(get_cube_power_at(cube, cnt, d, b, e));
            }
        }
    } else if (r + half_window > number_range_bin - 1) {
        for (size_t cnt = number_range_bin - 1 - range_CFAR_window_size; cnt < number_range_bin; ++cnt) {
            if (r != cnt) {
                CFAR_window.push_back(get_cube_power_at(cube, cnt, d, b, e));
            }
        }
    } else {
        for (size_t cnt = r - half_window; cnt < r + half_window + 1; ++cnt) {
            if (r != cnt) {
                CFAR_window.push_back(get_cube_power_at(cube, cnt, d, b, e));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(CFAR_window.begin(), CFAR_window.end());
    ref_value = CFAR_window[range_CFAR_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0) {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1/noise_value);
        ref_value = distribution(generator);
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * range_CFAR_factor;
}

/**
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in beam dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Threshold value for given cell
 */
float get_CFAR_threshold_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_azimuth_bin,
                              int azimuth_CFAR_window_size, int azimuth_CFAR_ref_pos, float noise_value,
                              float azimuth_CFAR_factor) {
    int8_t half_window = azimuth_CFAR_window_size / 2;
    std::vector<float> CFAR_window;
    int8_t window_size_buffer = 0;
    float ref_value;
    // Declaration of the reference-window. If no symmetric window around the test-cell is possible, the reference-
    // window is expanded on the other side of the spectrum, allowing a symmetric reference window in all cases:
    if (b < half_window) {
        for (size_t cnt = 0; cnt < b + half_window + 1; ++cnt) {
            if (cnt != b) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, cnt, e));
            }
        }
        window_size_buffer = CFAR_window.size();
        for (size_t cnt = number_azimuth_bin - azimuth_CFAR_window_size + window_size_buffer; cnt < number_azimuth_bin; ++cnt) {
            CFAR_window.push_back(get_cube_power_at(cube, r, d, cnt, e));
        }
    } else if (b + half_window > number_azimuth_bin - 1) {
        for (size_t cnt = b - half_window; cnt < number_azimuth_bin; ++cnt) {
            if (cnt != b) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, cnt, e));
            }
        }
        window_size_buffer = CFAR_window.size();
        for (size_t cnt = 0; cnt < azimuth_CFAR_window_size - window_size_buffer; ++cnt) {
            CFAR_window.push_back(get_cube_power_at(cube, r, d, cnt, e));
        }
    } else {
        for (size_t cnt = b - half_window; cnt < b + half_window + 1; ++cnt) {
            if (cnt != b) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, cnt, e));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(CFAR_window.begin(), CFAR_window.end());
    ref_value = CFAR_window[azimuth_CFAR_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0) {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1/noise_value);
        ref_value = distribution(generator);
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * azimuth_CFAR_factor;
}


/*!
 * Perform the ordered-statistics constant-false-alarm-rate (OS-CFAR) algorithm in elevation dimension. The reference
 * window is set around the cell-under-test (specified by r,d,b,e) and sorted. The reference value at the specified
 * position is used for the threshold calculation. The CFAR parameters are loaded from the parameter-files. Here the
 * size of the reference window, the position of the used reference value and the CFAR-constant are declared.
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @param number_elevation_bin
 * @param elevation_CFAR_window_size
 * @param elevation_CFAR_ref_pos
 * @param noise_value
 * @param elevation_CFAR_factor
 * @return
 */
float get_CFAR_threshold_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_elevation_bin,
                                   int elevation_CFAR_window_size, int elevation_CFAR_ref_pos, float noise_value,
                                   float elevation_CFAR_factor){
    int8_t half_window = elevation_CFAR_window_size / 2;
    std::vector<float> CFAR_window;
    int8_t window_size_buffer = 0;
    float ref_value;
    // Declaration of the reference-window. If no symmetric window around the test-cell is possible, the reference-
    // window is expanded on the other side of the spectrum, allowing a symmetric reference window in all cases:
    if (e < half_window) {
        for (size_t cnt = 0; cnt < e + half_window + 1; ++cnt) {
            if (cnt != e) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, b, cnt));
                std::cout<<"e < half_window, cnt: " + std::to_string(cnt)<<std::endl;
            }
        }
        window_size_buffer = CFAR_window.size();
        for (size_t cnt = number_elevation_bin - elevation_CFAR_window_size + window_size_buffer; cnt < number_elevation_bin; ++cnt) {
            CFAR_window.push_back(get_cube_power_at(cube, r, d, b, cnt));
        }
    } else if (e + half_window > number_elevation_bin - 1) {
        for (size_t cnt = e - half_window; cnt < number_elevation_bin; ++cnt) {
            if (cnt != e) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, b, cnt));
            }
        }
        window_size_buffer = CFAR_window.size();
        for (size_t cnt = 0; cnt < elevation_CFAR_window_size - window_size_buffer; ++cnt) {
            CFAR_window.push_back(get_cube_power_at(cube, r, d, b, cnt));
        }
    } else {
        for (size_t cnt = e - half_window; cnt < e + half_window + 1; ++cnt) {
            if (cnt != e) {
                CFAR_window.push_back(get_cube_power_at(cube, r, d, b, cnt));
            }
        }
    }
    // sort the reference window and get the reference value at the specified position of the ordered statistic
    std::sort(CFAR_window.begin(), CFAR_window.end());
    ref_value = CFAR_window[elevation_CFAR_ref_pos - 1];
    // if the reference value is zero, a random noise value is generated:
    if (ref_value == 0) {
        std::random_device generator;
        std::exponential_distribution<double> distribution(1/noise_value);
        ref_value = distribution(generator);
    }
    // return the threshold value by multiplying the reference value with the CFAR-constant:
    return ref_value * elevation_CFAR_factor;
}


/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in range dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * range dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Decimal part of the peak position in range dimension
 */
float correct_windowing_bias_interpolation_range(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_range_bin,
                                                 std::vector<float> lookup_range_interpolation, int lookup_size_interpolation,
                                                 float lookup_stepsize_interpolation) {
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t isNext;
    // Get maximum neighbour of the CUT in range-dimension, thereby check if CUT is at beginning or end of spectrum.
    // If at beginning or end only direct neighbour is considered.
    if (r == 0) {
        max_neighbor_pow = get_cube_power_at(cube, r + 1, d, b, e);
        isNext = 1;
    } else if (r == number_range_bin - 1) {
        max_neighbor_pow = get_cube_power_at(cube, r - 1, d, b, e);
        isNext = -1;
    } else {
        next_pow = get_cube_power_at(cube, r + 1, d, b, e);
        prev_pow = get_cube_power_at(cube, r - 1, d, b, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, r, d, b, e) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (!check_local_max_range(cube, r, d, b, e, get_cube_power_at(cube, r, d, b, e),  number_range_bin)) {
        sub_bin_correction = 0;
    } else if (ref_value < lookup_range_interpolation[lookup_size_interpolation - 1]) {
        sub_bin_correction = 0.5;
    } else if (ref_value > lookup_range_interpolation[0]) {
        sub_bin_correction = 0;
    } else {
        // search the position of the reference value inside the lookup-table
        for (auto &lu_step : lookup_range_interpolation) {
            if (ref_value > lu_step) {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * (lu_step_cnt - 1) +
                ((lookup_stepsize_interpolation * lu_step_cnt -
                lookup_stepsize_interpolation * (lu_step_cnt - 1)) /
                (lookup_range_interpolation[lu_step_cnt] - lookup_range_interpolation[lu_step_cnt - 1])) *
                (ref_value - lookup_range_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * isNext;  // return signed decimal bias estimation
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in doppler dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * doppler dimension. The function corrects the bias induced by the signal-processing of the radar under consideration
 * of the used window functions. The path to the used lookup-table is specified in the parameter files. For details on
 * the creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Decimal part of the peak position in doppler dimension
*/
float
correct_windowing_bias_interpolation_doppler(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_doppler_bin,
                                             std::vector<float> lookup_doppler_interpolation, int lookup_size_interpolation,
                                             float lookup_stepsize_interpolation) {
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t isNext;
    // get maximum neighbour of the CUT in doppler-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (d == 0) {
        next_pow = get_cube_power_at(cube, r, d + 1, b, e);
        prev_pow = get_cube_power_at(cube, r, number_doppler_bin - 1, b, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else if (d == number_doppler_bin - 1) {
        next_pow = get_cube_power_at(cube, r, 0, b, e);
        prev_pow = get_cube_power_at(cube, r, d - 1, b, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else {
        next_pow = get_cube_power_at(cube, r, d + 1, b, e);
        prev_pow = get_cube_power_at(cube, r, d - 1, b, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, r, d, b, e) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_doppler_interpolation[lookup_size_interpolation - 1]) {
        sub_bin_correction = 0.5;
    } else if (ref_value > lookup_doppler_interpolation[0]) {
        sub_bin_correction = 0;
    } else {
        // search the position of the reference value inside the lookup-table
        for (auto &lu_step : lookup_doppler_interpolation) {
            if (ref_value > lu_step) {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * (lu_step_cnt - 1) +
                ((lookup_stepsize_interpolation * lu_step_cnt -
                lookup_stepsize_interpolation * (lu_step_cnt - 1)) /
                (lookup_doppler_interpolation[lu_step_cnt] -
                lookup_doppler_interpolation[lu_step_cnt - 1])) *
                (ref_value - lookup_doppler_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * isNext;  // return signed decimal bias estimation
}

/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in beam dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * beam dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Decimal part of the peak position in beam dimension
 */
float
correct_windowing_bias_interpolation_beam(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_azimuth_bin,
                                          std::vector<float> lookup_azimuth_interpolation, int lookup_size_interpolation,
                                          float lookup_stepsize_interpolation) {
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t isNext;
    // get maximum neighbour of the CUT in beam-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (b == 0) {
        next_pow = get_cube_power_at(cube, r, d, b + 1, e);
        prev_pow = get_cube_power_at(cube, r, d, number_azimuth_bin - 1, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else if (b == number_azimuth_bin - 1) {
        next_pow = get_cube_power_at(cube, r, d, 0, e);
        prev_pow = get_cube_power_at(cube, r, d, b - 1, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else {
        next_pow = get_cube_power_at(cube, r, d, b + 1, e);
        prev_pow = get_cube_power_at(cube, r, d, b - 1, e);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, r, d, b, e) / max_neighbor_pow;
    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_azimuth_interpolation[lookup_size_interpolation - 1]) {
        sub_bin_correction = 0.5;
    } else if (ref_value > lookup_azimuth_interpolation[0]) {
        sub_bin_correction = 0;
    } else {
        // search the position of the reference value inside the lookup-table
        for (auto &lu_step : lookup_azimuth_interpolation) {
            if (ref_value > lu_step) {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * (lu_step_cnt - 1) +
                ((lookup_stepsize_interpolation * lu_step_cnt -
                lookup_stepsize_interpolation * (lu_step_cnt - 1)) /
                (lookup_azimuth_interpolation[lu_step_cnt] -
                lookup_azimuth_interpolation[lu_step_cnt - 1])) *
                (ref_value - lookup_azimuth_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * isNext;  // return signed decimal bias estimation
}


/**
 * Calculates decimal bias estimate using the pre-calculated lookup-table in elevation dimension.
 * Estimation is based on the peak value of the cell-under-test (CUT) at position r, d, b, e, and the maximum neighbour in
 * beam dimension. The function corrects the bias induced by the signal-processing of the radar under consideration of
 * the used window functions. The path to the used lookup-table is specified in the parameter files. For details on the
 * creation of the lookup-table see the corresponding matlab file.
 *
 * @param cube      The original radar cube
 * @param r         Range peak position
 * @param d         Doppler peak position
 * @param b         Beam peak position
 * @param e         Elevation peak position
 * @return
 *      Decimal part of the peak position in beam dimension
 */
float
correct_windowing_bias_interpolation_elevation(radar_cuboid &cube, size_t r, size_t d, size_t b, size_t e, int number_elevation_bin,
                                               std::vector<float> lookup_elevation_interpolation, int lookup_size_interpolation,
                                               float lookup_stepsize_interpolation){
    float ref_value;
    float max_neighbor_pow;
    float next_pow;
    float prev_pow;
    float sub_bin_correction;
    uint8_t lu_step_cnt = 0;
    int8_t isNext;
    // get maximum neighbour of the CUT in beam-dimension, thereby check if CUT is at beginning or end of
    // spectrum. If this is the case, check for maximum neighbour on the other side of the spectrum.
    if (e == 0) {
        next_pow = get_cube_power_at(cube, r, d, b, e + 1);
        prev_pow = get_cube_power_at(cube, r, d, b, number_elevation_bin - 1);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else if (e == number_elevation_bin - 1) {
        next_pow = get_cube_power_at(cube, r, d, b, 0);
        prev_pow = get_cube_power_at(cube, r, d, b, e - 1);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    } else {
        next_pow = get_cube_power_at(cube, r, d, b, e + 1);
        prev_pow = get_cube_power_at(cube, r, d, b, e - 1);
        if (prev_pow > next_pow) {
            max_neighbor_pow = prev_pow;
            isNext = -1;
        } else {
            max_neighbor_pow = next_pow;
            isNext = 1;
        }
    }
    // calculate reference value for bias estimation (CUT / max neighbour):
    ref_value = get_cube_power_at(cube, r, d, b, e) / max_neighbor_pow;

    // check if reference value is in bounds of lookup-table, if not set fixed bias estimation:
    if (ref_value < lookup_elevation_interpolation[lookup_size_interpolation - 1]) {
        sub_bin_correction = 0.5;
    } else if (ref_value > lookup_elevation_interpolation[0]) {
        sub_bin_correction = 0;
    } else {
        // search the position of the reference value inside the lookup-table
        for (auto &lu_step : lookup_elevation_interpolation) {
            if (ref_value > lu_step) {
                break;
            }
            ++lu_step_cnt;
        }
        // estimate the bias by interpolating between the two positions next to the reference value inside the
        // lookup table and the step-size of the lookup-table
        sub_bin_correction = lookup_stepsize_interpolation * (lu_step_cnt - 1) +
                             ((lookup_stepsize_interpolation * lu_step_cnt -
                               lookup_stepsize_interpolation * (lu_step_cnt - 1)) /
                              (lookup_elevation_interpolation[lu_step_cnt] -
                               lookup_elevation_interpolation[lu_step_cnt - 1])) *
                             (ref_value - lookup_elevation_interpolation[lu_step_cnt - 1]);
    }
    return sub_bin_correction * isNext;  // return signed decimal bias estimation

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
float get_cube_power_at(radar_cuboid &cube, const size_t range, const size_t doppler, const size_t beam,  const size_t elevation) {
    return static_cast<float>(pow(cube[range][doppler][beam][elevation][0], 2) + pow(cube[range][doppler][beam][elevation][1], 2));
}
