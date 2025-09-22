cmake_minimum_required(VERSION 3.10)
project(ReadHipoFile)
set(CMAKE_CXX_STANDARD 17)

# -------------------- ROOT --------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Matrix Graf Gpad Tree)
include(${ROOT_USE_FILE})

# Debug flags (keep yours)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# -------------------- RooUnfold --------------------
# Hints cover both an installed prefix and a build-tree config.
find_package(RooUnfold QUIET CONFIG
  HINTS
    ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install
    ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build
)

# -------------------- Executable & sources --------------------
add_executable(ReadHipoFile
  main.cpp
  Event.h  EventReader.cpp  Particle.cpp Event.cpp  EventReader.h Particle.h
  EventProcessor.cpp EventProcessor.h
  constants.h
  CutSet.h CutSet.cpp
  Monitoring.cpp Monitoring.h
  Ratio.h Ratio.cpp
  Dpt.cpp
  cratio.cpp cratio.h
  sratio.cpp sratio.h
  c2ratio.cpp c2ratio.h
  FilePathGenerator.cpp FilePathGenerator.h
  Monunfold.cpp Monunfold.h
)

# -------------------- Include paths for your local libs --------------------
target_include_directories(ReadHipoFile PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../lz4/
  ${CMAKE_CURRENT_SOURCE_DIR}/../hipo4/
)

# -------------------- Link (ROOT + local static libs) --------------------
target_link_libraries(ReadHipoFile
  PRIVATE
    ${ROOT_LIBRARIES}
    "${CMAKE_CURRENT_SOURCE_DIR}/../lib/libhipo4.a"
    "${CMAKE_CURRENT_SOURCE_DIR}/../lib/liblz4.a"
)

# -------------------- Link RooUnfold (handle all package variants) --------------------
if (RooUnfold_FOUND)
  message(STATUS "RooUnfold found: using CMake package")

  # Case 1: the package exports an imported target
  if (TARGET RooUnfold::RooUnfold)
    target_link_libraries(ReadHipoFile PRIVATE RooUnfold::RooUnfold)

  elseif (TARGET RooUnfold)
    target_link_libraries(ReadHipoFile PRIVATE RooUnfold)

  else()
    # Case 2: package sets variables, but not a target
    set(_ru_linked FALSE)

    # Try absolute library paths first
    foreach(_cand
      "${RooUnfold_LIBRARIES}"
      "${RooUnfold_LIBRARY}"
      "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib/libRooUnfold.so"
      "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib/libRooUnfold.a"
      "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.so"
      "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.a"
    )
      if (_cand AND EXISTS "${_cand}")
        message(STATUS "Linking RooUnfold using: ${_cand}")
        target_link_libraries(ReadHipoFile PRIVATE "${_cand}")
        get_filename_component(_ru_dir "${_cand}" DIRECTORY)
        # Set rpath so runtime can locate lib
        set_property(TARGET ReadHipoFile APPEND PROPERTY BUILD_RPATH "${_ru_dir}")
        set_property(TARGET ReadHipoFile APPEND PROPERTY INSTALL_RPATH "${_ru_dir}")
        set(_ru_linked TRUE)
        break()
      endif()
    endforeach()

    # Case 3: only a bare library name is available; add a link directory + rpath and link by name
    if (NOT _ru_linked)
      if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib")
        link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib)
        set_property(TARGET ReadHipoFile APPEND PROPERTY BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib")
        target_link_libraries(ReadHipoFile PRIVATE RooUnfold)
        set(_ru_linked TRUE)
      elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build")
        link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build)
        set_property(TARGET ReadHipoFile APPEND PROPERTY BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build")
        target_link_libraries(ReadHipoFile PRIVATE RooUnfold)
        set(_ru_linked TRUE)
      endif()
    endif()

    if (NOT _ru_linked)
      message(FATAL_ERROR "RooUnfold was found but no usable library path/target was provided.")
    endif()

    # Add include dirs if exposed by the package
    if (DEFINED RooUnfold_INCLUDE_DIRS)
      target_include_directories(ReadHipoFile PRIVATE ${RooUnfold_INCLUDE_DIRS})
    endif()
  endif()

else()
  message(WARNING "RooUnfold CMake package not found. Falling back to manual include/lib paths.")

  # Headers
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/include")
    target_include_directories(ReadHipoFile PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/include)
  else()
    target_include_directories(ReadHipoFile PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/src)
  endif()

  # Library (absolute path)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib/libRooUnfold.so")
    target_link_libraries(ReadHipoFile PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib/libRooUnfold.so)
    set_property(TARGET ReadHipoFile PROPERTIES
      BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib"
      INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/install/lib")
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.so")
    target_link_libraries(ReadHipoFile PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.so)
    set_property(TARGET ReadHipoFile PROPERTIES
      BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build")
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.a")
    target_link_libraries(ReadHipoFile PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/../RooUnfold/build/libRooUnfold.a)
  else()
    message(FATAL_ERROR "Could not find libRooUnfold.{so,a}. Build or install RooUnfold and rerun cmake.")
  endif()
endif()

