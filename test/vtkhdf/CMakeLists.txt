# Round trip through the fvtkhdf writer and VTK's own VTKHDF reader
# (built only with -DFVTKHTG_BUILD_VTKHDF_BRIDGE=ON).
#
#   vtkhdf_roundtrip_write  builds grids with fvtkhtg, writes them through the
#                           fvtkhtg_vtkhdf bridge, and records what every cell
#                           should be in a plain-text .expect file
#   vtkhdf_roundtrip_read   reads the files with vtkHDFReader (Python `vtk`
#                           9.7 or later, plus numpy-free standard library)
#                           and compares cell by cell; skipped (exit code 77)
#                           when that Python module is not available

if(fVTKHDF_MPI_SUPPORT)
  message(STATUS
    "fvtkhtg: skipping the VTKHDF round-trip test; it needs a serial fVTKHDF build")
  return()
endif()

add_executable(test_vtkhdf_roundtrip_write
  "${CMAKE_CURRENT_SOURCE_DIR}/test_vtkhdf_roundtrip_write.f90"
)
set_target_properties(test_vtkhdf_roundtrip_write PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
fvtkhtg_target_compile_flags(test_vtkhdf_roundtrip_write)
target_link_libraries(test_vtkhdf_roundtrip_write PRIVATE fvtkhtg_vtkhdf)

set(_fvtkhtg_roundtrip_dir "${CMAKE_CURRENT_BINARY_DIR}/roundtrip")
file(MAKE_DIRECTORY "${_fvtkhtg_roundtrip_dir}")

add_test(NAME vtkhdf_roundtrip_write
  COMMAND test_vtkhdf_roundtrip_write
  WORKING_DIRECTORY "${_fvtkhtg_roundtrip_dir}"
)
set_tests_properties(vtkhdf_roundtrip_write PROPERTIES
  FIXTURES_SETUP vtkhdf_roundtrip_files
  LABELS "vtkhdf"
)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME vtkhdf_roundtrip_read
    COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/check_vtkhdf_roundtrip.py" "${_fvtkhtg_roundtrip_dir}"
  )
  set_tests_properties(vtkhdf_roundtrip_read PROPERTIES
    FIXTURES_REQUIRED vtkhdf_roundtrip_files
    SKIP_RETURN_CODE 77
    LABELS "vtkhdf"
  )
else()
  message(STATUS "fvtkhtg: Python 3 not found; the VTKHDF files will be written but not read back")
endif()
