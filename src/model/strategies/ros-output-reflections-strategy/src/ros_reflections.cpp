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

#include "ros_reflections/ros_reflections.hpp"
#include <sensor_msgs/PointCloud2.h>
#include <pcl_ros/point_cloud.h>
#include <string>
#include <utility>

#ifdef _WIN32
#include <math.h>
#else
#include <cmath>
#endif

using namespace model;
using namespace osi3;

void ros_reflections::apply(SensorData &sensor_data) {
    log("Starting ROS output for reflections");

    ros::param::set("/use_sim_time", true);

    if ((sensor_data.sensor_view_size()==0) || (!sensor_data.has_feature_data())) {
        log("No sensor view or feature data received");
        return;
    }

    if (!sensor_data.sensor_view(0).has_global_ground_truth()) {
        log("No global ground truth received");
        return;
    }

    if (sensor_data.sensor_view(0).radar_sensor_view().size() > 0) {
        for (int sensor_no = 0; sensor_no < sensor_data.sensor_view(0).radar_sensor_view().size(); sensor_no++) {
            worker_pcl = std::make_unique<reflections::WorkerPCL>("reflections_" + std::to_string(sensor_no), "sensor_" + std::to_string(sensor_no));
            worker_pcl->injectRadar(sensor_data, sensor_no, log);
        }
    } else {
        log("No radar sensor view");
        return;
    }
}

ros_reflections::ros_reflections(const Profile &profile, const Log &log, const Alert &alert): model::Strategy(profile, log, alert) {
    auto remapping = std::map<std::string, std::string>();
    remapping.emplace("__master", "http://localhost:11311");
    ros::init(remapping, "sensor_model_fmu");
}

reflections::WorkerPCL::WorkerPCL(const std::string& topic, std::string frame_id) : publisher(node.advertise<pcl::PointCloud<pcl::PointXYZ>>(topic, 1)) , frame_id(std::move(frame_id)) {}

void reflections::WorkerPCL::injectRadar(SensorData &sensor_data, int sensor_no, const Log &log) {
    //loop over all sensors and rendering results

    const auto &radar_sensor_view = sensor_data.sensor_view(0).radar_sensor_view(sensor_no);
    auto no_of_rendering_results = radar_sensor_view.reflection_size();

    auto sim_seconds = uint32_t(sensor_data.sensor_view(0).global_ground_truth().timestamp().seconds());
    auto sim_nanos = uint32_t(sensor_data.sensor_view(0).global_ground_truth().timestamp().nanos());
    auto sim_time = ros::Time(sim_seconds, sim_nanos);
    auto sim_time_behind = ros::Time(sim_seconds, sim_nanos + 500);

    if (no_of_rendering_results == 0) {
        auto timestamp = (double) sim_seconds + (double) sim_nanos / 1000000000;
        log("No reflection from sensor " + std::to_string(sensor_no) + " for ROS output at timestamp " +
            std::to_string(timestamp));
        return;
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tmp(new pcl::PointCloud<pcl::PointXYZ>());
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>());

    pcl_conversions::toPCL(sim_time_behind, cloud_tmp->header.stamp);
    cloud_tmp->header.frame_id = frame_id;
    /// Loop over all received sensor views
    for (const auto &reflection: radar_sensor_view.reflection()) {
        auto ray_vertical_angle_rad = reflection.source_vertical_angle();
        auto ray_horizontal_angle_rad = reflection.source_horizontal_angle();
        auto distance = 0.5 * reflection.time_of_flight() * Speed_of_Light;

        double x = distance * cos(ray_horizontal_angle_rad) * cos(ray_vertical_angle_rad);
        double y = distance * sin(ray_horizontal_angle_rad) * cos(ray_vertical_angle_rad);
        double z = distance * sin(ray_vertical_angle_rad);
        cloud_tmp->points.emplace_back(pcl::PointXYZ(float(x), float(y), float(z)));
    }
    pcl::copyPointCloud(*cloud_tmp, *cloud);
    for (int reflection_idx = 0; reflection_idx < radar_sensor_view.reflection().size(); reflection_idx++) {
        auto signal_strength_in_dB = radar_sensor_view.reflection(
                reflection_idx).signal_strength();
        cloud->points[reflection_idx].intensity = float(signal_strength_in_dB);
    }
    publisher.publish(cloud);
}
