#pragma once


#ifndef INCLUDE_GUARD
#define INCLUDE_GUARD


#define REMFracRes_NAME ""
#define REMFracRes_VER  "1.1.1.3"
#define REMFracRes_VER_MAJOR "1"
#define REMFracRes_VER_MINOR "1"
#define REMFracRes_VER_PATCH "1"
#include <string>


/* Absolute path to test source directory */
namespace remfracres{
static constexpr auto data_path = "/home/pascal/CIG/Workspaces/cpp_workspaces/openGeode/Data/FracRes/";
bool  test_fractures_boolean();
}
 
#endif // INCLUDE_GUARD
