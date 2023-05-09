cmake_minimum_required(VERSION 3.10)
project(ReadHipoFile)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS)

include(${ROOT_USE_FILE})

add_executable(ReadHipoFile main.cpp dcFiducialCuts.hpp variablecalc.cpp linspace.hpp histograms.hpp)

target_link_libraries(ReadHipoFile "${CMAKE_CURRENT_SOURCE_DIR}/lib/libhipo4.a")
target_link_libraries(ReadHipoFile "${CMAKE_CURRENT_SOURCE_DIR}/lib/liblz4.a")
target_link_libraries(ReadHipoFile ${ROOT_LIBRARIES})

target_include_directories(ReadHipoFile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lz4/)
target_include_directories(ReadHipoFile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/hipo4/)
