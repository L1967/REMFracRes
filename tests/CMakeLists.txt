# Copyright (c) 2019 - 2022 Geode-solutions
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.14)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Release)
set(ENV{PATH_HOME_WORK} /home/pascal/CIG/Workspaces/cpp_workspaces/OpenGeode_New) 
set(ENV{OpenGeode_DIR} $ENV{PATH_HOME_WORK}/OpenGeode-12.2.7)
set(ENV{OpenGeode-IO_DIR} $ENV{PATH_HOME_WORK}/OpenGeode-IO-5.22.6)
set(ENV{OpenGeode-Geosciences_DIR} $ENV{PATH_HOME_WORK}/OpenGeode-Geosciences-5.7.5)
set(ENV{OpenGeode-GeosciencesIO_DIR} $ENV{PATH_HOME_WORK}/OpenGeode-GeosciencesIO-3.12.0)

set(CMAKE_PREFIX_PATH $ENV{OpenGeode_DIR}/build/opengeode/install $ENV{OpenGeode-IO_DIR}/build/opengeode-io/install $ENV{OpenGeode-Geosciences_DIR}/build $ENV{OpenGeode-GeosciencesIO_DIR}/build)
    
if(NOT TARGET REMFracRes::mylib)
    project(REMFracRes version "1.0.0.1")
    find_package(OpenGeode REQUIRED CONFIG)
    find_package(OpenGeode-Geosciences REQUIRED CONFIG)
    find_package(OpenGeode-IO REQUIRED CONFIG)
    find_package(OpenGeode-GeosciencesIO REQUIRED CONFIG)   
    find_package(REMFracRes REQUIRED CONFIG)
    enable_testing()
endif()

if(NOT TARGET REMFracRes::fractureslib)
    project(REMFracRes version "1.0.0.1")
    find_package(OpenGeode REQUIRED CONFIG)
    find_package(OpenGeode-Geosciences REQUIRED CONFIG)
    find_package(OpenGeode-IO REQUIRED CONFIG)
    find_package(OpenGeode-GeosciencesIO REQUIRED CONFIG)  
    find_package(REMFracRes REQUIRED CONFIG)
    enable_testing()
endif()

add_subdirectory(mylib)
add_subdirectory(fractureslib)

