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

cmake_minimum_required(VERSION 3.0)
# Define the project
project(REMFracRes  VERSION "1.1.1.3")
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Debug)


set(CMAKE_PREFIX_PATH /home/pascal/DevPro/Total_Flowers/install-dir)

option(REMFRACRES_WITH_TESTS "Compile test projects" ON)
option(REMFRACRES_WITH_PYTHON "Compile Python bindings" OFF)

# Get MyModule geode dependencies
find_package(OpenGeode REQUIRED)
find_package(OpenGeode-IO REQUIRED)
find_package(OpenGeode-Geosciences REQUIRED)
find_package(OpenGeode-GeosciencesIO REQUIRED)

############# DATA Directory for test and Version in config.h.in
set(DATA_DIRECTORY /home/pascal/CIG/Workspaces/cpp_workspaces/openGeode/Data/FracRes)
configure_file("${CMAKE_SOURCE_DIR}/src/remfracres/REMFracResConfig.h.in" "${CMAKE_SOURCE_DIR}/include/remfracres/REMFracRes.h")



#---------------------------------------------------------------------------------------------

# Compile options
add_compile_options(-v -std=c++11 )

#c ++ Pattern 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

add_geode_library(
    NAME lib
    FOLDER remfracres
    SOURCES
    	"FracResUnit.cpp"
    	"RemFracResFileManagement.cpp"
        "FracResBeddingAperture.cpp"
        "FracResBeddingProperties.cpp"
        "FracResFracture.cpp"
        "fractureset.cpp"
        "FaultArraySet.cpp"
        "FracResHypothesis.cpp"
        "FracResStage.cpp"
        "FracResScenario.cpp"
        "RemFracResSim.cpp"
        "../primitive_mesh/FracResBeddingLocalisation.cpp"
        "../primitive_mesh/boite.cpp"
        "../primitive_mesh/faille.cpp"
        "../primitive_mesh/fonction.cpp"
        "../primitive_mesh/cylindre.cpp"
        "../primitive_mesh/plan.cpp"
        "../primitive_mesh/point.cpp"
        "../primitive_mesh/polyedre.cpp"
        "../primitive_mesh/predicates.cpp"
        "../primitive_mesh/primitive.cpp"
        "../primitive_mesh/primitive_derive.cpp"
        "../primitive_mesh/contrainte_geo.cpp"
        "../primitive_mesh/fracture_set_geo.cpp"
        "../primitive_mesh/fault_array_set_geo.cpp"
        "../primitive_mesh/fault_array.cpp"
        "../numeric/model.cpp"
    PUBLIC_HEADERS
        "FracResUnit.h"
        "RemFracResFileManagement.h"
        "FracResBeddingAperture.h"
        "FracResBeddingProperties.h"
        "FracResFracture.h"
        "fractureset.h"
        "FaultArraySet.h"
        "FracResHypothesis.h"
        "FracResStage.h"
        "FracResScenario.h"
        "RemFracResSim.h"
        "../primitive_mesh/varglob.h"
        "../primitive_mesh/FracResBeddingLocalisation.h"
        "../primitive_mesh/bounding_box.h"
        "../primitive_mesh/oriented_bounding_box.h"
        "../primitive_mesh/boundingbox.h"
        "../primitive_mesh/boite.h"
        "../primitive_mesh/tetra.h"
        "../primitive_mesh/faille.h"
        "../primitive_mesh/fonction.h"
        "../primitive_mesh/cylindre.h"
        "../primitive_mesh/plan.h"
        "../primitive_mesh/point.h"
        "../primitive_mesh/polyedre.h"
        "../primitive_mesh/predicates.h"
        "../primitive_mesh/primitive.h"
        "../primitive_mesh/primitive_derive.h"
        "../primitive_mesh/contrainte_geo.h"
        "../primitive_mesh/fracture_set_geo.h"
        "../primitive_mesh/fault_arrays_set_geo.h"
        "../primitive_mesh/fault_array.h"
        "../remfracres/fracResFileKeyWord.h"
        "../numeric/model.h"
    PRIVATE_DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        OpenGeode::mesh
        OpenGeode::model
        OpenGeode-IO::mesh
        OpenGeode-IO::model
        OpenGeode-Geosciences::explicit
        OpenGeode-Geosciences::implicit
        OpenGeode-GeosciencesIO::mesh
        OpenGeode-GeosciencesIO::model
)


add_geode_binary(
    SOURCE "src/remfracres/RemFracResSim.cpp"
    DEPENDENCIES
        OpenGeode::basic
        OpenGeode::geometry
        OpenGeode::mesh
        OpenGeode::model
        OpenGeode-IO::mesh
        OpenGeode-IO::model
        OpenGeode-Geosciences::explicit
        OpenGeode-Geosciences::implicit
        OpenGeode-GeosciencesIO::mesh
        OpenGeode-GeosciencesIO::model
        lib
)
set_target_properties(RemFracResSim PROPERTIES VERSION ${PROJECT_VERSION})

############ OTHER WAY FOR CREATION EXECUTABLE AND LIBRARY AND SWIG ENCAPSULATION STAND ALONE NOT USED GEODE CMAKE COMMAND
file(GLOB SOURCESCPP1 ${CMAKE_SOURCE_DIR}/src/remfracres/*.cpp)
file(GLOB_RECURSE GEODE_LIB ${CMAKE_PREFIX_PATH}/lib64/*.so)
file(GLOB_RECURSE GEODE_INCLUDE ${CMAKE_PREFIX_PATH}/include/*.h)
file(GLOB SOURCESPRIMITIVECPP1 ${CMAKE_SOURCE_DIR}/src/primitive_mesh/*.cpp)
file(GLOB SOURCESNUMERICCPP1 ${CMAKE_SOURCE_DIR}/src/numeric/*.cpp)
file(GLOB SOURCESFRACRESCPP1 ${CMAKE_SOURCE_DIR}/src/remfracres/FracResUnit.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/RemFracResFileManagement.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/FracResBeddingAperture.cpp   ${CMAKE_SOURCE_DIR}/src/remfracres/FracResBeddingProperties.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/FracResFracture.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/FracResScenario.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/FracResHypothesis.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/FracResStage.cpp)

file(GLOB_RECURSE SRCCPP ${CMAKE_PREFIX_PATH}/src/*.cpp)
# ---------------------------------------------executable creation need include directories

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/remfracres)
include_directories(${CMAKE_SOURCE_DIR}/include/numeric)
include_directories(${CMAKE_SOURCE_DIR}/include/primitive_mesh)

# GEODE INCLUDE
include_directories(${CMAKE_PREFIX_PATH}/include)
include_directories(${CMAKE_PREFIX_PATH}/include/absl)
include_directories(${CMAKE_PREFIX_PATH}/include/geode)
include_directories(${CMAKE_PREFIX_PATH}/include/bitsery)
include_directories(${CMAKE_PREFIX_PATH}/include/geode/geosciences)
include_directories(${CMAKE_PREFIX_PATH}/include/ghc)
include_directories(${CMAKE_PREFIX_PATH}/include/async++)

#add_executable(REMFracResApp ${SOURCESNUMERICCPP1} ${SOURCESPRIMITIVECPP1} ${SOURCESFRACRESCPP1} ${CMAKE_SOURCE_DIR}/src/remfracres/RemFracResSim.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/fractureset.cpp  ${CMAKE_SOURCE_DIR}/src/remfracres/FaultArraySet.cpp)
#TARGET_LINK_LIBRARIES(REMFracResApp ${GEODE_LIB})
#set_target_properties(REMFracResApp PROPERTIES VERSION ${PROJECT_VERSION})

# ---------------------------------------------library creation If NEEDED
#add_library(REMFracResApp_lib SHARED ${GEODE_INCLUDE} ${SOURCESFRACRESCPP1} ${CMAKE_SOURCE_DIR}/src/remfracres/RemFracResSim.cpp ${CMAKE_SOURCE_DIR}/src/remfracres/fractureset.cpp  ${CMAKE_SOURCE_DIR}/src/remfracres/FaultArraySet.cpp ${SOURCESNUMERICCPP1}  ${SOURCESPRIMITIVECPP1})
#set_target_properties(REMFracResApp_lib PROPERTIES VERSION ${PROJECT_VERSION})
#TARGET_LINK_LIBRARIES(REMFracResApp_lib ${GEODE_LIB})

# ---------------------------------------------swig creation
set(ENV{JAVA_HOME} /home/pascal/CIG/devkit/jdk8) 
# Get Package Java
find_package(Java 1.8 REQUIRED COMPONENTS Runtime Development)
find_package(JNI REQUIRED)
include(UseJava)

#Get package SWIG
find_package(SWIG REQUIRED)
include(UseSWIG)
include(${SWIG_USE_FILE})
#include directories
include_directories(${JNI_INCLUDE_DIRS})
include_directories($ENV{JAVA_HOME}/include)
include_directories($ENV{JAVA_HOME}/include/linux)

#--------------swig using
# Get devkitjdk dependency

# get Java tools
set(Java_JAVA_EXECUTABLE /home/pascal/CIG/devkit/jdk8/bin/java) 
MESSAGE ( STATUS "  JAVA_VERSION:           " ${Java_VERSION})
MESSAGE ( STATUS "  JAVA_RUNTIME:           " ${Java_JAVA_EXECUTABLE})
MESSAGE ( STATUS "  JAVA_COMPILE:           " ${Java_JAVAC_EXECUTABLE})
MESSAGE ( STATUS "  JAVA_ARCHIVE:           " ${Java_JAR_EXECUTABLE})
MESSAGE ( STATUS "  JAVA_JAVA_H:            " ${Java_JAVAH_EXECUTABLE})
MESSAGE ( STATUS "  CMAKE_CURRENT_BINARY_DIR:            " ${CMAKE_CURRENT_BINARY_DIR})
#MESSAGE ( FATAL_ERROR "  LIBRARY_OUTPUT_DIRECTORY:            " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#Set SWIG parameters
# Specify an output directory name , Used to place the generated source file . For example, output exampleJAVA_wrap.cxx Location of files, etc 
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/java/src/REMFracRes)
set(SWIG_OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/srcWrap)

# package
set(CMAKE_SWIG_FLAGS -package REMFracRes)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_SOURCE_DIR}/swig/REMFracRes.i PROPERTIES CPLUSPLUS ON)

# To start generating swig library
swig_add_library(REMFracResJava_lib LANGUAGE java SOURCES ${CMAKE_SOURCE_DIR}/swig/REMFracRes.i ${SRCCPP})
swig_link_libraries(REMFracResJava_lib ${GEODE_LIB} ${JNI_LIBRARIES})


# custom command
set(GRAPHSIM_SWIG_JAVA_CLASS_DIR ${CMAKE_SOURCE_DIR}/java/src) 
set(REMFRACRES_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/java/src) 
file(GLOB SWIG_JAVA_FILES )
add_custom_command(
     TARGET REMFracResJava_lib
     POST_BUILD
     COMMAND cmake -E echo "Compiling Java files..."
     COMMAND cmake -E echo "Compiling Java  and class files..." 
     COMMAND ${Java_JAVAC_EXECUTABLE} -d ${REMFRACRES_BINARY_DIR} -sourcepath ${CMAKE_SWIG_OUTDIR} ${REMFRACRES_BINARY_DIR}/REMFracRes/*.java
     COMMAND cmake -E echo "making Jar package..."
     COMMAND cmake -E make_directory ${REMFRACRES_BINARY_DIR}/REMFracRes/jar
     COMMAND ${Java_JAR_EXECUTABLE} -cvf REMFracRes.jar ${REMFRACRES_BINARY_DIR}/REMFracRes/*.class
     
    )  
#      



# For convenience we copy the dynamic library to the current build folder

#set(LIBFRACRESSIM_DIR /home/pascal/CIG/Workspaces/workspace-git/SmgLibrariesREMFracRes/REMFracRes/install/LinuxOP) 
#set(LIBFRACRESSIM_LIB ${LIBFRACRESSIM_DIR}/libREMFracRes.so) 
#file(MAKE_DIRECTORY ${LIBFRACRESSIM_DIR})
#file(GLOB SOURCESJAVACLASS ${REMFRACRES_BINARY_DIR}/*.class)



# Configure version
set(REMFRACRES_VER_MAJOR 0)
set(REMFRACRES_VER_MINOR 25)
set(REMFRACRES_VER_PATCH 3)

#------------------------------------------------------------------------------------------------
# Configure CPacks
if(WIN32)
    set(CPACK_GENERATOR "ZIP")
else()
    set(CPACK_GENERATOR "TGZ")
endif()

# This must always be last!
include(CPack)