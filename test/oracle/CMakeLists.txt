# C++ VTK v9.7+ Differential & LLDB Lockstep Test Suite (`test/oracle/`)
# Built strictly when `-DFVTKHTG_BUILD_ORACLE_TESTS=ON`; `libfvtkhtg` in `src/`
# remains 100% pure Fortran 2018 with zero C/C++ bindings.

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/src/fvtkhtg.f90")
  message(STATUS "fvtkhtg oracle suite: awaiting ${CMAKE_SOURCE_DIR}/src/fvtkhtg.f90")
  return()
endif()

add_library(fvtkhtg_test_c_bridge STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/fvtkhtg_test_c_bridge.f90"
)
set_target_properties(fvtkhtg_test_c_bridge PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_bridge"
  POSITION_INDEPENDENT_CODE ON
)
target_compile_options(fvtkhtg_test_c_bridge PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-g>
  $<$<COMPILE_LANGUAGE:Fortran>:-O0>
)
fvtkhtg_target_compile_flags(fvtkhtg_test_c_bridge)
target_link_libraries(fvtkhtg_test_c_bridge PUBLIC fvtkhtg)

function(fvtkhtg_add_oracle_executable target_name cxx_src)
  add_executable(${target_name} "${CMAKE_CURRENT_SOURCE_DIR}/${cxx_src}")
  target_include_directories(${target_name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_options(${target_name} PRIVATE
    -g
    -O0
    -fno-inline
    -Wall
    -Wextra
    $<$<BOOL:${FVTKHTG_WARNINGS_AS_ERRORS}>:-Werror>
  )
  target_link_libraries(${target_name} PRIVATE
    fvtkhtg_test_c_bridge
    VTK::CommonCore
    VTK::CommonDataModel
  )
  set_target_properties(${target_name} PROPERTIES
    BUILD_RPATH "${FVTKHTG_VTK_ORACLE_INSTALL_DIR}/lib"
  )
  add_test(NAME ${target_name} COMMAND ${target_name})
  set_tests_properties(${target_name} PROPERTIES
    LABELS "oracle"
  )
endfunction()

fvtkhtg_add_oracle_executable(
  test_oracle_dimensions_coords
  "test_oracle_dimensions_coords.cxx"
)
fvtkhtg_add_oracle_executable(
  test_oracle_amr_indexing
  "test_oracle_amr_indexing.cxx"
)
fvtkhtg_add_oracle_executable(
  test_oracle_mask_celldata
  "test_oracle_mask_celldata.cxx"
)
fvtkhtg_add_oracle_executable(
  test_oracle_bf_descriptors
  "test_oracle_bf_descriptors.cxx"
)
fvtkhtg_add_oracle_executable(
  test_oracle_geometry_neighbors
  "test_oracle_geometry_neighbors.cxx"
)
fvtkhtg_add_oracle_executable(
  test_lockstep_driver
  "test_lockstep_driver.cxx"
)

if(EXISTS "${CMAKE_SOURCE_DIR}/tools/lldb/verify_lockstep.py")
  if(NOT FVTKHTG_UV_EXECUTABLE)
    message(FATAL_ERROR "uv executable is required to run tools/lldb/verify_lockstep.py (Python 3.14)")
  endif()
  add_test(
    NAME test_lldb_lockstep
    COMMAND "${FVTKHTG_UV_EXECUTABLE}" run python
      "${CMAKE_SOURCE_DIR}/tools/lldb/verify_lockstep.py"
      "$<TARGET_FILE:test_lockstep_driver>"
      --fortran-archive "$<TARGET_FILE:fvtkhtg>"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  )
  set_tests_properties(test_lldb_lockstep PROPERTIES
    TIMEOUT 180
    LABELS "lldb;oracle"
  )
endif()
