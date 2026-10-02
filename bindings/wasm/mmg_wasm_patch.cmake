# Patch MMG's Fortran-header macro to run a host-compiled genheader (see MvbMeshWasm.cmake).
# Usage: cmake -DSRC=<mmg source dir> -P mmg_wasm_patch.cmake
# ExternalProject re-runs the patch step whenever its stamps are reset (a refetch, a reconfigure),
# so the file is restored to the checked-out v5.8.0 text first: the patch always applies to the
# same input, and an upstream change still stops the build below.
execute_process(COMMAND git -C "${SRC}" checkout -- cmake/modules/macros.cmake RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "mmg_wasm_patch: cannot restore ${SRC}/cmake/modules/macros.cmake from git")
endif()
file(READ "${SRC}/cmake/modules/macros.cmake" _m)
string(FIND "${_m}" "COMMAND genheader " _at)
if(_at EQUAL -1)
  message(FATAL_ERROR "mmg_wasm_patch: 'COMMAND genheader ' not found in macros.cmake -- MMG changed, update the patch")
endif()
string(REPLACE "COMMAND genheader " "COMMAND \${NATIVE_GENHEADER} " _m "${_m}")
string(REPLACE "DEPENDS genheader " "DEPENDS " _m "${_m}")
file(WRITE "${SRC}/cmake/modules/macros.cmake" "${_m}")
