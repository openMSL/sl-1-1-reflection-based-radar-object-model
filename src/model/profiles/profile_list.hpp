//
// Copyright 2022 Technical University of Darmstadt - FZD
// SPDX-License-Identifier: MPL-2.0
//

#ifndef OSMPSENSORFRAMEWORK_PROFILE_LIST_HPP
#define OSMPSENSORFRAMEWORK_PROFILE_LIST_HPP

#include <string>

#include <model/profiles/profile.hpp>

/* TODO add further profiles here */
#include <model/profiles/profile_ReflectionRadar.hpp>

bool CFrameworkPackaging::try_load_profile(const std::string& name)
{
    /* TODO add further profile generators here */
    if (name == "ReflectionRadar")
    {
        profile = model::profile::ReflectionRadar::generate();
        return true;
    }
    return false;
}

#endif  // OSMPSENSORFRAMEWORK_PROFILE_LIST_HPP
