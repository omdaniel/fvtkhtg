# Minimal Debug VTK v9.7+ C++ Oracle locator (gated strictly behind FVTKHTG_BUILD_ORACLE_TESTS)

set(FVTKHTG_VTK_ORACLE_GIT_URL "https://gitlab.kitware.com/vtk/vtk.git")
set(FVTKHTG_VTK_ORACLE_TAG "v9.7.2")
set(FVTKHTG_VTK_ORACLE_ROOT "${CMAKE_SOURCE_DIR}/.cache/vtk-oracle")
set(FVTKHTG_VTK_ORACLE_SRC_DIR "${FVTKHTG_VTK_ORACLE_ROOT}/src")
set(FVTKHTG_VTK_ORACLE_BUILD_DIR "${FVTKHTG_VTK_ORACLE_ROOT}/build")
set(FVTKHTG_VTK_ORACLE_INSTALL_DIR "${FVTKHTG_VTK_ORACLE_ROOT}/install")

if(NOT EXISTS "${FVTKHTG_VTK_ORACLE_INSTALL_DIR}")
  message(FATAL_ERROR
    "VTK ${FVTKHTG_VTK_ORACLE_TAG} oracle install not found at:\n"
    "  ${FVTKHTG_VTK_ORACLE_INSTALL_DIR}\n"
    "Build the minimal Debug oracle (CommonCore + CommonDataModel) from ${FVTKHTG_VTK_ORACLE_GIT_URL} @ ${FVTKHTG_VTK_ORACLE_TAG} "
    "or configure with -DFVTKHTG_BUILD_ORACLE_TESTS=OFF for pure-Fortran builds."
  )
endif()

find_package(VTK REQUIRED
  COMPONENTS
    CommonCore
    CommonDataModel
  PATHS "${FVTKHTG_VTK_ORACLE_INSTALL_DIR}"
  NO_DEFAULT_PATH
)

if(VTK_VERSION VERSION_LESS "9.7.0")
  message(FATAL_ERROR
    "Found VTK oracle version ${VTK_VERSION} at ${FVTKHTG_VTK_ORACLE_INSTALL_DIR}, "
    "but VTK >= 9.7.0 (${FVTKHTG_VTK_ORACLE_TAG}) is strictly required."
  )
endif()

message(STATUS "Found VTK C++ Oracle ${VTK_VERSION} at ${FVTKHTG_VTK_ORACLE_INSTALL_DIR}")

find_program(FVTKHTG_UV_EXECUTABLE
  NAMES uv
  HINTS "$ENV{HOME}/.local/bin" "/opt/homebrew/bin"
)
if(FVTKHTG_UV_EXECUTABLE)
  message(STATUS "Found uv at ${FVTKHTG_UV_EXECUTABLE}")
endif()

