# Strict Fortran 2018 compiler flags and diagnostics for fvtkhtg targets.
#
# Flags are applied per target (PRIVATE) by fvtkhtg_target_compile_flags(), so
# a parent project that pulls fvtkhtg in with add_subdirectory or
# FetchContent does not inherit -std=f2018, -Werror or -fcheck=all.

option(FVTKHTG_WARNINGS_AS_ERRORS "Treat Fortran compiler warnings as errors" OFF)

function(fvtkhtg_target_compile_flags target)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:-std=f2018>
      $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
      $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
      $<$<COMPILE_LANGUAGE:Fortran>:-Wimplicit-interface>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${FVTKHTG_WARNINGS_AS_ERRORS}>>:-Werror>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=all>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(LLVMFlang|Flang)$")
    target_compile_options(${target} PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:-std=f2018>
      $<$<COMPILE_LANGUAGE:Fortran>:-pedantic>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${FVTKHTG_WARNINGS_AS_ERRORS}>>:-Werror>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    )
  endif()
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endfunction()
