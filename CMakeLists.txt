cmake_minimum_required(VERSION 3.10)
project(ReadHipoFile)

set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS)

include(${ROOT_USE_FILE})

# Enable debugging information
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

add_executable(ReadHipoFile reboot.cpp EventProcessor.cpp dcFiducialCuts.hpp variablecalc.cpp linspace.hpp histograms.hpp)
# add_executable(ReadHipoFile main.cpp dcFiducialCuts.hpp variablecalc.cpp linspace.hpp histograms.hpp)


target_link_libraries(ReadHipoFile ${ROOT_LIBRARIES})
target_link_libraries(ReadHipoFile "${CMAKE_CURRENT_SOURCE_DIR}/lib/libhipo4.a")
target_link_libraries(ReadHipoFile "${CMAKE_CURRENT_SOURCE_DIR}/lib/liblz4.a")

# Remove the following line:
# target_include_directories(ReadHipoFile PRIVATE /home/matamoros/Downloads/root/include)

# Specify target-specific include directories
target_include_directories(ReadHipoFile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lz4/)
target_include_directories(ReadHipoFile PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/hipo4/)


