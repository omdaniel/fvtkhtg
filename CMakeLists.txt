cmake_minimum_required(VERSION 3.24)

project(fvtkhtg
  VERSION 0.1.0
  DESCRIPTION "Pure Fortran 2018 vtkHyperTreeGrid Construction and Validation Library"
  LANGUAGES Fortran
)

# Tests and examples default ON only when fvtkhtg is the top-level project,
# so add_subdirectory / FetchContent consumers get just the library.
option(FVTKHTG_BUILD_TESTING "Build pure-Fortran unit tests" ${PROJECT_IS_TOP_LEVEL})
option(FVTKHTG_BUILD_ORACLE_TESTS "Build C++ VTK v9.7+ differential oracle tests" OFF)
option(FVTKHTG_BUILD_EXAMPLES "Build runnable fvtkhtg example programs" ${PROJECT_IS_TOP_LEVEL})
option(FVTKHTG_BUILD_VTKHDF_BRIDGE
  "Build the optional fvtkhtg_vtkhdf library that writes grids with fvtkhdf (requires fVTKHDF)" OFF)
option(FVTKHTG_INSTALL "Generate install rules and a CMake package config" ${PROJECT_IS_TOP_LEVEL})

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompilerWarnings)

# Watch canonical directories so Ninja re-runs CMake automatically during incremental development
file(GLOB _fvtkhtg_incremental_watch CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/test/unit/test_*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/test/oracle/CMakeLists.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/examples/*.f90"
)

set(FVTKHTG_CANONICAL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_kinds_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_bit_array_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_cell_data_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_breadth_first_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_tree_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_grid_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_cursor_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_cursor_geometry_s.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_validate_m.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg_grid_export_s.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg.f90"
)

set(FVTKHTG_ACTIVE_SOURCES "")
foreach(_src IN LISTS FVTKHTG_CANONICAL_SOURCES)
  if(EXISTS "${_src}")
    list(APPEND FVTKHTG_ACTIVE_SOURCES "${_src}")
  else()
    message(STATUS "fvtkhtg incremental build: pending module ${_src}")
  endif()
endforeach()

if(NOT FVTKHTG_ACTIVE_SOURCES)
  message(FATAL_ERROR "No Fortran source files found for target fvtkhtg in ${CMAKE_CURRENT_SOURCE_DIR}/src")
endif()

add_library(fvtkhtg STATIC ${FVTKHTG_ACTIVE_SOURCES})
add_library(fvtkhtg::fvtkhtg ALIAS fvtkhtg)
set_target_properties(fvtkhtg PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
fvtkhtg_target_compile_flags(fvtkhtg)
target_include_directories(fvtkhtg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fvtkhtg>
)

if(FVTKHTG_BUILD_TESTING OR FVTKHTG_BUILD_ORACLE_TESTS OR
   (FVTKHTG_BUILD_EXAMPLES AND BUILD_TESTING))
  enable_testing()
endif()

function(fvtkhtg_add_unit_test test_name test_src required_mod_src)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_src}" AND
     EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${required_mod_src}")
    add_executable(${test_name} "${CMAKE_CURRENT_SOURCE_DIR}/${test_src}")
    set_target_properties(${test_name} PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_${test_name}"
    )
    fvtkhtg_target_compile_flags(${test_name})
    target_link_libraries(${test_name} PRIVATE fvtkhtg)
    add_test(NAME ${test_name} COMMAND ${test_name})
  endif()
endfunction()

if(FVTKHTG_BUILD_TESTING)
  fvtkhtg_add_unit_test(
    test_foundations
    "test/unit/test_foundations.f90"
    "src/fvtkhtg_cell_data_m.f90"
  )
  fvtkhtg_add_unit_test(
    test_tree
    "test/unit/test_tree.f90"
    "src/fvtkhtg_tree_m.f90"
  )
  fvtkhtg_add_unit_test(
    test_grid_cursor_validate
    "test/unit/test_grid_cursor_validate.f90"
    "src/fvtkhtg.f90"
  )
  fvtkhtg_add_unit_test(
    test_export
    "test/unit/test_export.f90"
    "src/fvtkhtg_grid_export_s.f90"
  )
  fvtkhtg_add_unit_test(
    test_geometry_neighbors
    "test/unit/test_geometry_neighbors.f90"
    "src/fvtkhtg_cursor_geometry_s.f90"
  )
  fvtkhtg_add_unit_test(
    test_coarsen
    "test/unit/test_coarsen.f90"
    "src/fvtkhtg_cursor_geometry_s.f90"
  )
endif()

function(fvtkhtg_add_example example_name example_src)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${example_src}" AND
     EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/fvtkhtg.f90")
    set(_example_target "example_${example_name}")
    add_executable(${_example_target} "${CMAKE_CURRENT_SOURCE_DIR}/${example_src}")
    set_target_properties(${_example_target} PROPERTIES
      OUTPUT_NAME "${example_name}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/examples"
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_${_example_target}"
    )
    fvtkhtg_target_compile_flags(${_example_target})
    target_link_libraries(${_example_target} PRIVATE fvtkhtg)
    if((DEFINED BUILD_TESTING AND BUILD_TESTING) OR
       (NOT DEFINED BUILD_TESTING AND FVTKHTG_BUILD_TESTING))
      add_test(NAME ${_example_target} COMMAND $<TARGET_FILE:${_example_target}>)
      set_tests_properties(${_example_target} PROPERTIES LABELS "examples")
    endif()
  endif()
endfunction()

if(FVTKHTG_BUILD_EXAMPLES)
  fvtkhtg_add_example("01_basic_quadtree" "examples/01_basic_quadtree.f90")
  fvtkhtg_add_example("02_octree_amr_celldata" "examples/02_octree_amr_celldata.f90")
  fvtkhtg_add_example("03_mask_puremask_interface" "examples/03_mask_puremask_interface.f90")
  fvtkhtg_add_example("04_bf_descriptors_and_views" "examples/04_bf_descriptors_and_views.f90")
  fvtkhtg_add_example("05_breadth_first_export" "examples/05_breadth_first_export.f90")
  fvtkhtg_add_example("06_adaptive_refinement" "examples/06_adaptive_refinement.f90")
endif()

# Optional bridge to the fvtkhdf HyperTreeGrid writer. The core library never
# depends on fvtkhdf; only this separate target does.
if(FVTKHTG_BUILD_VTKHDF_BRIDGE)
  # fVTKHDF's package finds HDF5, whose own package needs the C language
  # enabled (for FindThreads). The bridge itself is pure Fortran.
  enable_language(C)
  find_package(fVTKHDF 0.8 CONFIG REQUIRED)
  add_library(fvtkhtg_vtkhdf STATIC "${CMAKE_CURRENT_SOURCE_DIR}/bridge/vtkhdf/fvtkhtg_vtkhdf.f90")
  add_library(fvtkhtg::fvtkhtg_vtkhdf ALIAS fvtkhtg_vtkhdf)
  set_target_properties(fvtkhtg_vtkhdf PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_vtkhdf"
  )
  fvtkhtg_target_compile_flags(fvtkhtg_vtkhdf)
  target_include_directories(fvtkhtg_vtkhdf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_vtkhdf>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fvtkhtg>
  )
  target_link_libraries(fvtkhtg_vtkhdf PUBLIC fvtkhtg fVTKHDF::fvtkhdf)

  # Round trip: write files through the bridge, read them back with VTK
  if(FVTKHTG_BUILD_TESTING)
    add_subdirectory(test/vtkhdf)
  endif()
endif()

if(FVTKHTG_INSTALL)
  include(CMakePackageConfigHelpers)
  set(_fvtkhtg_cmake_dir "${CMAKE_INSTALL_LIBDIR}/cmake/fvtkhtg")

  install(TARGETS fvtkhtg EXPORT fvtkhtgTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  # Fortran .mod files are compiler-specific: consumers must use the same
  # compiler (fvtkhtgConfig.cmake warns otherwise).
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/fvtkhtg"
    FILES_MATCHING PATTERN "*.mod"
  )
  install(EXPORT fvtkhtgTargets
    FILE fvtkhtgTargets.cmake
    NAMESPACE fvtkhtg::
    DESTINATION "${_fvtkhtg_cmake_dir}"
  )

  if(FVTKHTG_BUILD_VTKHDF_BRIDGE)
    install(TARGETS fvtkhtg_vtkhdf EXPORT fvtkhtgVtkhdfTargets
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_vtkhdf/"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/fvtkhtg"
      FILES_MATCHING PATTERN "*.mod"
    )
    install(EXPORT fvtkhtgVtkhdfTargets
      FILE fvtkhtgVtkhdfTargets.cmake
      NAMESPACE fvtkhtg::
      DESTINATION "${_fvtkhtg_cmake_dir}"
    )
  endif()

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fvtkhtgConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/fvtkhtgConfig.cmake"
    INSTALL_DESTINATION "${_fvtkhtg_cmake_dir}"
  )
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/fvtkhtgConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
  )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/fvtkhtgConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/fvtkhtgConfigVersion.cmake"
    DESTINATION "${_fvtkhtg_cmake_dir}"
  )
endif()

if(FVTKHTG_BUILD_ORACLE_TESTS)
  enable_language(CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  include(VtkOracle)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/oracle/CMakeLists.txt")
    add_subdirectory(test/oracle)
  else()
    message(STATUS "fvtkhtg oracle mode enabled: awaiting test/oracle/CMakeLists.txt")
  endif()
endif()
