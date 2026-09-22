cmake_minimum_required(VERSION 3.18)
project(mvbpp_wasm VERSION 0.1.0 LANGUAGES CXX)

# Must be invoked with the Emscripten toolchain:
#   emcmake cmake -S . -B build-wasm
#   cmake --build build-wasm
if(NOT EMSCRIPTEN)
    message(FATAL_ERROR "This CMakeLists.txt must be built with emcmake.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
include(ExternalProject)

# ── OpenCASCADE 7.9.3 — compiled for WASM ────────────────────────────────────
set(OCCT_PREFIX "${CMAKE_BINARY_DIR}/occt-install")

ExternalProject_Add(
    occt_wasm
    URL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_9_3.tar.gz
    CMAKE_ARGS
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_LIBRARY_TYPE=Static
        -DBUILD_MODULE_Draw=OFF
        -DBUILD_MODULE_DETools=OFF
        -DBUILD_DOC_Overview=OFF
        -DBUILD_PERF_TESTS=OFF
        -DBUILD_TESTING=OFF
        -DBUILD_SAMPLES_MFC=OFF
        -DBUILD_SAMPLES_QT=OFF
        -DBUILD_Inspector=OFF
        -DBUILD_RESOURCES=OFF
        -DUSE_FREETYPE=OFF
        -DUSE_FREEIMAGE=OFF
        -DUSE_TBB=OFF
        -DUSE_OPENVR=OFF
        -DCMAKE_INSTALL_PREFIX=${OCCT_PREFIX}
        -DINSTALL_DIR_LIB=lib
        -DINSTALL_DIR_INCLUDE=include/opencascade
    BUILD_BYPRODUCTS
        ${OCCT_PREFIX}/lib/libTKernel.a
        ${OCCT_PREFIX}/lib/libTKMath.a
        ${OCCT_PREFIX}/lib/libTKG2d.a
        ${OCCT_PREFIX}/lib/libTKG3d.a
        ${OCCT_PREFIX}/lib/libTKGeomBase.a
        ${OCCT_PREFIX}/lib/libTKBRep.a
        ${OCCT_PREFIX}/lib/libTKGeomAlgo.a
        ${OCCT_PREFIX}/lib/libTKTopAlgo.a
        ${OCCT_PREFIX}/lib/libTKPrim.a
        ${OCCT_PREFIX}/lib/libTKBO.a
        ${OCCT_PREFIX}/lib/libTKBool.a
        ${OCCT_PREFIX}/lib/libTKMesh.a
        ${OCCT_PREFIX}/lib/libTKShHealing.a
        ${OCCT_PREFIX}/lib/libTKCDF.a
        ${OCCT_PREFIX}/lib/libTKLCAF.a
        ${OCCT_PREFIX}/lib/libTKXSBase.a
        ${OCCT_PREFIX}/lib/libTKXCAF.a
        ${OCCT_PREFIX}/lib/libTKCAF.a
        ${OCCT_PREFIX}/lib/libTKVCAF.a
        ${OCCT_PREFIX}/lib/libTKBin.a
        ${OCCT_PREFIX}/lib/libTKXml.a
        ${OCCT_PREFIX}/lib/libTKStd.a
        ${OCCT_PREFIX}/lib/libTKOffset.a
        ${OCCT_PREFIX}/lib/libTKRWMesh.a
        ${OCCT_PREFIX}/lib/libTKDESTEP.a
        ${OCCT_PREFIX}/lib/libTKDEIGES.a
        ${OCCT_PREFIX}/lib/libTKDESTL.a
    USES_TERMINAL_BUILD TRUE
)

add_library(occt INTERFACE)
add_dependencies(occt occt_wasm)
target_include_directories(occt INTERFACE ${OCCT_PREFIX}/include/opencascade)
target_link_directories(occt INTERFACE ${OCCT_PREFIX}/lib)

set(OCCT_ESSENTIAL_LIBS
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKPrim TKBO TKBool TKMesh TKShHealing
    TKCDF TKLCAF TKXSBase TKXCAF TKCAF TKVCAF TKBin TKXml TKStd
    TKFillet TKOffset TKRWMesh TKDESTEP TKDEIGES TKDESTL
)
target_link_libraries(occt INTERFACE ${OCCT_ESSENTIAL_LIBS})

# ── nlohmann/json ─────────────────────────────────────────────────────────────
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json/
    GIT_TAG        v3.11.3
)
FetchContent_MakeAvailable(json)

# ── MKF ───────────────────────────────────────────────────────────────────────
set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS"   CACHE PATH "" FORCE)
set(MAS_DIR       "${CMAKE_BINARY_DIR}/_deps/mas-src" CACHE PATH "" FORCE)
set(INCLUDE_MKF_TESTS OFF CACHE BOOL "" FORCE)
set(INCLUDE_PYMKF     OFF CACHE BOOL "" FORCE)
set(ENABLE_MATPLOTPP  OFF CACHE BOOL "" FORCE)
set(ENABLE_NGSPICE    OFF CACHE BOOL "" FORCE)

# ABT #852 — do NOT embed the finished-core inventory in the drawing engine.
#
# MKF's cmrc `data` archive is force-linked whole (see --whole-archive below), so every
# catalog listed in MKF's DATA_LIST lands in the wasm DATA section whether it is read or
# not. mvbpp DRAWS geometry; it never looks up a core part number: nothing in MVB++ or in
# its bindings calls load_cores / find_core_by_name / get_core_shape_names(manufacturer),
# and no adviser or cross-referencer (the only MKF callers of those) is reachable from a
# draw entry point. The two inventory files are 9.48 MB of the engine — 19% of it — read
# by nobody:
#     MAS/data/cores.ndjson        8.20 MB
#     MAS/data/cores_stock.ndjson  1.28 MB
# Measured on the 2026-08-22 engine: 49.75 MB -> 40.27 MB raw (-19.0%), 10.52 -> 10.14 MB
# gzipped, with byte-identical drawCore/drawTurns/drawMagnetic output on cm37 and the
# tests/mas_complete_fixtures battery.
#
# The catalogs mvbpp DOES read stay embedded: core_shapes (Core::resolve_shape),
# core_materials (Core::resolve_material), wires, wire_materials, bobbins and
# insulation_materials all resolve name references in the MAS files the app hands us.
set(EMBED_MAS_CORES OFF CACHE BOOL "" FORCE)

set(_MAS_LOCAL "https://github.com/OpenMagnetics/MAS.git")
set(_MKF_LOCAL "https://github.com/OpenMagnetics/MKF.git")

FetchContent_Declare(
    MAS
    GIT_REPOSITORY "${_MAS_LOCAL}"
    # Track MAS main (kept in lockstep with the MKF main pin below).
    GIT_TAG        main
    GIT_SHALLOW    FALSE
)
FetchContent_MakeAvailable(MAS)
set(MAS_DIR "${mas_SOURCE_DIR}" CACHE PATH "" FORCE)

set(_MAS_SRC "${CMAKE_BINARY_DIR}/_deps/mas-src")
set(_MKF_SRC "${CMAKE_BINARY_DIR}/_deps/mkf-src")
set(_MKF_BIN "${CMAKE_BINARY_DIR}/_deps/mkf-build")

FetchContent_Declare(
    MKF
    GIT_REPOSITORY "${_MKF_LOCAL}"
    # Track MKF main (kept in lockstep with the MAS main pin above).
    GIT_TAG        main
    GIT_SHALLOW    FALSE
    GIT_SUBMODULES ""
)
FetchContent_GetProperties(MKF)
if(NOT mkf_POPULATED)
    FetchContent_Populate(MKF)
endif()
if(IS_DIRECTORY "${_MKF_SRC}/MAS" AND NOT IS_SYMLINK "${_MKF_SRC}/MAS")
    file(REMOVE_RECURSE "${_MKF_SRC}/MAS")
endif()
if(NOT EXISTS "${_MKF_SRC}/MAS")
    file(CREATE_LINK "${_MAS_SRC}" "${_MKF_SRC}/MAS" SYMBOLIC)
endif()
# Initialize submodules MKF needs (cci_coords for codegen, CAS for CAS.hpp,
# PEAS for the PEAS-integrated schemas). MKF main no longer has an EAS submodule,
# so listing it would fail the whole command and leave CAS uninitialised.
if(IS_DIRECTORY "${_MKF_SRC}/cci_coords" AND NOT EXISTS "${_MKF_SRC}/cci_coords/coordinates")
    execute_process(
        COMMAND git submodule update --init cci_coords CAS PEAS
        WORKING_DIRECTORY "${_MKF_SRC}"
    )
endif()
if(NOT TARGET MKF)
    # MKF hardcodes add_library(MKF SHARED ...) — patch it to STATIC.
    file(READ "${_MKF_SRC}/CMakeLists.txt" _mkf_cmake)
    string(REPLACE "add_library(MKF SHARED" "add_library(MKF STATIC" _mkf_cmake "${_mkf_cmake}")
    file(WRITE "${_MKF_SRC}/CMakeLists.txt" "${_mkf_cmake}")
    include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann")
    add_subdirectory("${_MKF_SRC}" "${_MKF_BIN}" EXCLUDE_FROM_ALL)
endif()

# ── MVB++ sources (shared with native build) ──────────────────────────────────
set(REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")

add_library(mvb++ STATIC
    ${REPO_ROOT}/src/MagneticBuilder.cpp
    ${REPO_ROOT}/src/SectionBuilder.cpp
    ${REPO_ROOT}/src/SectionDrawing.cpp
    ${REPO_ROOT}/src/ProjectionDrawing.cpp
    ${REPO_ROOT}/src/SpacerBuilder.cpp
    ${REPO_ROOT}/src/FR4Builder.cpp
    # ABT #1169's accessory builders: the fourth, fifth and sixth sources to reach the native
    # build without this list. The guard below caught them at CONFIGURE this time, which is why
    # the web engine has been stuck at 2026-08-23 rather than linking against a stale half.
    ${REPO_ROOT}/src/BaseBuilder.cpp
    ${REPO_ROOT}/src/PinBuilder.cpp
    ${REPO_ROOT}/src/ShuntBuilder.cpp
    ${REPO_ROOT}/src/StepExporter.cpp
    ${REPO_ROOT}/src/Symmetry.cpp
    ${REPO_ROOT}/src/Utils.cpp
    ${REPO_ROOT}/src/TurnBuilder.cpp
    ${REPO_ROOT}/src/ConductorBuilder.cpp
    ${REPO_ROOT}/src/TerminalFillet.cpp
    ${REPO_ROOT}/src/WireAssembler.cpp
    ${REPO_ROOT}/src/BobbinBuilder.cpp
    ${REPO_ROOT}/src/shapes/ShapeBuilder.cpp
    ${REPO_ROOT}/src/shapes/ShapeDrum.cpp
    ${REPO_ROOT}/src/shapes/ShapeDrumRing.cpp
    ${REPO_ROOT}/src/shapes/ShapeMolded.cpp
    ${REPO_ROOT}/src/shapes/ShapeDrumSemishielded.cpp
    ${REPO_ROOT}/src/shapes/ShapeE.cpp
    ${REPO_ROOT}/src/shapes/ShapeEp.cpp
    ${REPO_ROOT}/src/shapes/ShapeEpx.cpp
    ${REPO_ROOT}/src/shapes/ShapeLp.cpp
    ${REPO_ROOT}/src/shapes/ShapeEc.cpp
    ${REPO_ROOT}/src/shapes/ShapeEr.cpp
    ${REPO_ROOT}/src/shapes/ShapeP.cpp
    ${REPO_ROOT}/src/shapes/ShapePm.cpp
    ${REPO_ROOT}/src/shapes/ShapePQ.cpp
    ${REPO_ROOT}/src/shapes/ShapeRM.cpp
    ${REPO_ROOT}/src/shapes/ShapeT.cpp
    ${REPO_ROOT}/src/shapes/ShapeU.cpp
    ${REPO_ROOT}/src/shapes/ShapeUr.cpp
    ${REPO_ROOT}/src/shapes/ShapeUt.cpp
    ${REPO_ROOT}/src/shapes/Factory.cpp
)

# This source list DUPLICATES the native build's, and it has now drifted THREE times -- missing
# ConductorBuilder.cpp (the wasm was unbuildable from main since June), then ShapeDrum.cpp
# (wasm-ld: "undefined symbol: vtable for mvb::shapes::ShapeDrum"), then TerminalFillet.cpp
# (undefined smoothLayerLinks / filletTerminalCorners). Every one surfaced only at LINK time,
# long after the source landed.
#
# The third one got through a guard written to prevent exactly this, because that guard globbed
# only src/shapes/*.cpp -- and TerminalFillet.cpp lives in src/. A guard narrower than the list
# it protects gives false confidence: the configure passed, so the omission looked impossible.
# So the glob is now RECURSIVE over every source the native build compiles.
get_target_property(_mvb_wasm_srcs mvb++ SOURCES)
file(GLOB_RECURSE _mvb_all_srcs "${REPO_ROOT}/src/*.cpp")
foreach(_s ${_mvb_all_srcs})
    if(NOT "${_s}" IN_LIST _mvb_wasm_srcs)
        message(FATAL_ERROR
            "bindings/wasm/CMakeLists.txt does not list ${_s}. Every src/**/*.cpp must appear in "
            "the wasm mvb++ source list, or the link fails with an undefined symbol or vtable. "
            "If a source is deliberately native-only, exclude it here explicitly rather than by "
            "omission, so the intent is visible.")
    endif()
endforeach()


target_include_directories(mvb++ PUBLIC
    ${REPO_ROOT}/include
    ${MAS_DIRECTORY}
    ${json_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/_deps/spline-src/src
    ${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include/magic_enum
    ${CMAKE_BINARY_DIR}/_deps/svg-src/src
)
target_link_libraries(mvb++ PUBLIC occt MKF nlohmann_json::nlohmann_json)

# ── WASM module ───────────────────────────────────────────────────────────────
add_executable(mvbpp_wasm mvbpp_wasm.cpp)

# cmrc embeds resources via static constructors in libdata.a.  When linking
# static libraries the linker drops unreferenced object files, so the resource
# initializers never run.  Force the entire archive so all cmrc files are kept.
target_link_libraries(mvbpp_wasm PRIVATE mvb++)
target_link_options(mvbpp_wasm PRIVATE
    "LINKER:--whole-archive"
    "$<TARGET_FILE:data>"
    "LINKER:--no-whole-archive"
)

# Enable C++ exception catching only in mvbpp_wasm.cpp so the guard wrappers
# at each entry point can catch std::exception thrown from MKF/OCCT (which
# already throw using the default legacy EH ABI). Scoped to this object so
# we avoid the perf hit of recompiling OCCT/MKF/mvb++ with -fexceptions.
target_compile_options(mvbpp_wasm PRIVATE -fexceptions)
target_link_options(mvbpp_wasm PRIVATE -fexceptions)

set_target_properties(mvbpp_wasm PROPERTIES
    OUTPUT_NAME "mvbpp"
    SUFFIX      ".js"
)

# Emscripten link flags
target_link_options(mvbpp_wasm PRIVATE
    -lembind
    "SHELL:-s MODULARIZE=1"
    "SHELL:-s EXPORT_NAME=createMvbpp"
    "SHELL:-s ALLOW_MEMORY_GROWTH=1"
    "SHELL:-s INITIAL_MEMORY=256MB"
    "SHELL:-s STACK_SIZE=8MB"
    # Emscripten virtual FS for /tmp writes
    "SHELL:-s FORCE_FILESYSTEM=1"
    "SHELL:-s EXPORTED_RUNTIME_METHODS=['FS']"
    # WASM output alongside the JS glue
    "SHELL:-s WASM=1"
)
