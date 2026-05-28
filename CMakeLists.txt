cmake_minimum_required(VERSION 3.18)
project(MVB++ VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

include(FetchContent)
include(ExternalProject)

# ============================================================================
# OpenCASCADE 7.9.3 — build from source, minimal configuration
# ============================================================================
set(OCCT_PREFIX "${CMAKE_BINARY_DIR}/occt-install")

ExternalProject_Add(
    occt_external
    URL https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V7_9_3.tar.gz
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
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
add_dependencies(occt occt_external)
target_include_directories(occt INTERFACE ${OCCT_PREFIX}/include/opencascade)
target_link_directories(occt INTERFACE ${OCCT_PREFIX}/lib)

set(OCCT_ESSENTIAL_LIBS
    TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKPrim TKBO TKBool TKMesh TKShHealing
    TKCDF TKLCAF TKXSBase TKXCAF TKCAF TKVCAF TKBin TKXml TKStd
    TKOffset TKRWMesh TKDESTEP TKDEIGES TKDESTL
    TKV3d TKService TKHLR TKFillet TKFeat
    TKBinL TKStdL TKXmlL TKBinXCAF TKXmlXCAF
    TKDE TKDECascade
    TKTObj TKBinTObj TKXmlTObj
    TKDEGLTF TKDEOBJ TKDEPLY TKDEVRML
)
target_link_libraries(occt INTERFACE ${OCCT_ESSENTIAL_LIBS})

# ============================================================================
# nlohmann/json — reused by both MKF and MVB++
# ============================================================================
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json/
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(json)

# ============================================================================
# Catch2 v3 for tests
# ============================================================================
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.8.0
)
FetchContent_MakeAvailable(Catch2)

# ============================================================================
# MKF — brings MAS.hpp (via quicktype) and physics/autocomplete logic
# Fetched from GitHub so the build is self-contained.
# ============================================================================
set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS" CACHE PATH "" FORCE)
set(INCLUDE_MKF_TESTS OFF CACHE BOOL "" FORCE)
set(INCLUDE_PYMKF     OFF CACHE BOOL "" FORCE)
set(ENABLE_MATPLOTPP  OFF CACHE BOOL "" FORCE)
set(ENABLE_NGSPICE    OFF CACHE BOOL "" FORCE)

# MKF main on GitHub depends on a MAS commit that is local-only (not yet pushed).
# Use file:// clones from the local repos so FetchContent pulls the exact working
# versions.  Switch both to https://github.com/OpenMagnetics/{MAS,MKF}.git once
# the pending MAS commit is pushed.
set(_MAS_LOCAL "https://github.com/OpenMagnetics/MAS.git")
# Pin to the last MKF commit compatible with the public MAS schema.
# Bump to main once MKF commits that reference BobbinFamily::T and
# get_inner_phase_shift3 are matched by a pushed MAS schema update.
set(_MKF_LOCAL "https://github.com/OpenMagnetics/MKF.git")

FetchContent_Declare(
    MAS
    GIT_REPOSITORY "${_MAS_LOCAL}"
    # Bumped to MAS origin/main HEAD: PFC topology variants, AHB/PSFB schema
    # updates, dutyCycle bounds, bulkCapacitance, plus MAS 1.0 enum migration
    # (twoPieceSet, wound, innerOrTop, etc. via migrate-to-1.0.py).
    GIT_TAG        a16a6007fa364e4af2fc81a1dc895da4fe393e7c
    GIT_SHALLOW    FALSE
)
FetchContent_MakeAvailable(MAS)
set(MAS_DIR "${mas_SOURCE_DIR}" CACHE PATH "" FORCE)

# MKF hard-codes MAS/data/ relative to its own source dir. We must plant the
# MAS symlink BEFORE add_subdirectory runs. FetchContent always places sources
# at _deps/<lowercase_name>-src, so we can use that convention directly.
set(_MAS_SRC "${CMAKE_BINARY_DIR}/_deps/mas-src")
set(_MKF_SRC "${CMAKE_BINARY_DIR}/_deps/mkf-src")
set(_MKF_BIN "${CMAKE_BINARY_DIR}/_deps/mkf-build")

FetchContent_Declare(
    MKF
    GIT_REPOSITORY "${_MKF_LOCAL}"
    # Bumped to MKF origin/main HEAD to match the MAS 2f37eea pin above.
    # PtP harness additions, topology expansions. MKF's magnetic_autocomplete
    # and Core/Coil constructors are what MVB++ actually uses.
    GIT_TAG        424bf300022a07b73b9bafd90a53bb817514fef4
    GIT_SHALLOW    FALSE
    GIT_SUBMODULES ""
)
FetchContent_GetProperties(MKF)
if(NOT mkf_POPULATED)
    FetchContent_Populate(MKF)
endif()
# git creates an empty MAS/ placeholder for the uninitialized submodule;
# remove it so we can plant a symlink to the separately-fetched MAS.
if(IS_DIRECTORY "${_MKF_SRC}/MAS" AND NOT IS_SYMLINK "${_MKF_SRC}/MAS")
    file(REMOVE_RECURSE "${_MKF_SRC}/MAS")
endif()
if(NOT EXISTS "${_MKF_SRC}/MAS")
    file(CREATE_LINK "${_MAS_SRC}" "${_MKF_SRC}/MAS" SYMBOLIC)
endif()
# Initialize submodules MKF needs (cci_coords for codegen, CAS for Topology.h)
if(IS_DIRECTORY "${_MKF_SRC}/cci_coords" AND NOT EXISTS "${_MKF_SRC}/cci_coords/coordinates")
    execute_process(
        COMMAND git submodule update --init cci_coords CAS EAS
        WORKING_DIRECTORY "${_MKF_SRC}"
        RESULT_VARIABLE _cci_result
    )
    if(NOT _cci_result EQUAL 0)
        message(WARNING "Failed to init MKF submodules (exit ${_cci_result})")
    endif()
endif()
if(NOT TARGET MKF)
    # MKF hardcodes add_library(MKF SHARED ...) — patch it to STATIC so the
    # Python wheel and WASM module are fully self-contained.
    file(READ "${_MKF_SRC}/CMakeLists.txt" _mkf_cmake)
    string(REPLACE "add_library(MKF SHARED" "add_library(MKF STATIC" _mkf_cmake "${_mkf_cmake}")
    # MKF generates CAS.hpp via the CAS custom target but forgets to declare
    # MKF depending on it — first parallel build fails with "CAS.hpp: No such
    # file or directory". Patch the dependency list to include CAS.
    string(REPLACE
        "add_dependencies(MKF MAS cci_data_gen)"
        "add_dependencies(MKF MAS cci_data_gen CAS)"
        _mkf_cmake "${_mkf_cmake}")
    file(WRITE "${_MKF_SRC}/CMakeLists.txt" "${_mkf_cmake}")
    include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann")
    add_subdirectory("${_MKF_SRC}" "${_MKF_BIN}" EXCLUDE_FROM_ALL)
endif()

# ============================================================================
# MVB++ library
# ============================================================================
# When building the Python wheel we want a self-contained mvbpp.so with no
# external libmvb++.so dependency.  Build mvb++ as STATIC + position-independent
# so it can be linked into the pybind11 module.
if(BUILD_PYTHON_BINDINGS)
    set(_mvbpp_lib_kind STATIC)
else()
    set(_mvbpp_lib_kind SHARED)
endif()
add_library(mvb++ ${_mvbpp_lib_kind}
    src/MagneticBuilder.cpp
    src/SectionBuilder.cpp
    src/SectionDrawing.cpp
    src/ProjectionDrawing.cpp
    src/SpacerBuilder.cpp
    src/FR4Builder.cpp
    src/StepExporter.cpp
    src/Symmetry.cpp
    src/Utils.cpp
    src/TurnBuilder.cpp
    src/BobbinBuilder.cpp
    src/shapes/ShapeBuilder.cpp
    src/shapes/ShapeE.cpp
    src/shapes/ShapeEp.cpp
    src/shapes/ShapeEpx.cpp
    src/shapes/ShapeLp.cpp
    src/shapes/ShapeEc.cpp
    src/shapes/ShapeEr.cpp
    src/shapes/ShapeP.cpp
    src/shapes/ShapePm.cpp
    src/shapes/ShapePQ.cpp
    src/shapes/ShapeRM.cpp
    src/shapes/ShapeT.cpp
    src/shapes/ShapeU.cpp
    src/shapes/ShapeUr.cpp
    src/shapes/ShapeUt.cpp
    src/shapes/Factory.cpp
)
set_target_properties(mvb++ PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(mvb++ PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${MAS_DIRECTORY}
    ${json_SOURCE_DIR}/include/nlohmann
    ${json_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/_deps/spline-src/src
    ${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include/magic_enum
    ${CMAKE_BINARY_DIR}/_deps/svg-src/src
)

target_link_libraries(mvb++ PUBLIC
    occt
    MKF
    nlohmann_json::nlohmann_json
)

# ============================================================================
# Tests
# ============================================================================
enable_testing()
add_executable(mvb_tests
    tests/test_main.cpp
    tests/test_core_shapes.cpp
    tests/test_all_shapes.cpp
    tests/test_all_gapped.cpp
    tests/test_step_comparison.cpp
    tests/test_integration.cpp
    tests/test_json_integration.cpp
    tests/test_symmetry.cpp
    tests/test_mas_battery.cpp
    tests/test_toroidal_topview.cpp
    tests/test_toroidal_single_turn_centered.cpp
)

target_compile_definitions(mvb_tests PRIVATE
    REFERENCE_STEPS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/reference_steps"
    MAS_EXAMPLES_DIR="${CMAKE_BINARY_DIR}/_deps/mas-src/examples"
    MAS_DATA_DIR="${CMAKE_BINARY_DIR}/_deps/mas-src/data"
    # Locally-patched copies of MAS/examples/complete with stale wire/shape
    # refs corrected to valid DB entries (see tools/fix-mas-complete-examples.py).
    MAS_COMPLETE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/mas_complete_fixtures"
)

target_link_libraries(mvb_tests PRIVATE mvb++ Catch2::Catch2WithMain)

add_executable(mvbpp_generate_steps tests/generate_mvbpp_steps.cpp)
target_link_libraries(mvbpp_generate_steps PRIVATE mvb++)

add_executable(mvbpp_step_generator tools/mvbpp_step_generator.cpp)
target_link_libraries(mvbpp_step_generator PRIVATE mvb++)

add_executable(mvbpp_dump_turns tools/dump_turns.cpp)
target_link_libraries(mvbpp_dump_turns PRIVATE mvb++)

add_executable(mvbpp_section tools/mvbpp_section.cpp)
target_link_libraries(mvbpp_section PRIVATE mvb++)

add_executable(mvbpp_generate_bobbins tools/generate_mvbpp_bobbins.cpp)
target_link_libraries(mvbpp_generate_bobbins PRIVATE occt)

add_test(NAME mvb_tests COMMAND mvb_tests)

# ============================================================================
# Python bindings (optional — cmake -DBUILD_PYTHON_BINDINGS=ON)
# ============================================================================
option(BUILD_PYTHON_BINDINGS "Build the pybind11 Python module (mvbpp)" OFF)
if(BUILD_PYTHON_BINDINGS)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)
    pybind11_add_module(mvbpp bindings/python/mvbpp_python.cpp)
    install(TARGETS mvbpp
            LIBRARY DESTINATION .
            COMPONENT PythonModule
            EXCLUDE_FROM_ALL)

    # Static OCCT requires --whole-archive when linked into a shared object
    # to pull in all vtable/typeinfo symbols the linker would otherwise drop.
    set(_occt_static_libs)
    foreach(_lib ${OCCT_ESSENTIAL_LIBS})
        list(APPEND _occt_static_libs "${OCCT_PREFIX}/lib/lib${_lib}.a")
    endforeach()
    target_link_options(mvbpp PRIVATE
        "LINKER:--whole-archive"
        ${_occt_static_libs}
        "LINKER:--no-whole-archive"
    )

    target_include_directories(mvbpp PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${MAS_DIRECTORY}
        ${json_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/_deps/spline-src/src
        ${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include/magic_enum
    )
    target_link_libraries(mvbpp PRIVATE
        "-Wl,--whole-archive" mvb++ "-Wl,--no-whole-archive")

    # Static OCCT's TKService references X11 symbols (XOpenDisplay,
    # XGetVisualInfo, XGetWindowAttributes) for offscreen visual queries.
    # Without an explicit link line, the .so loads but throws
    # "undefined symbol: XGetVisualInfo" on import, forcing users to
    # LD_PRELOAD libX11. Link X11 directly (Linux only).
    if(UNIX AND NOT APPLE)
        find_package(X11 REQUIRED)
        target_link_libraries(mvbpp PRIVATE ${X11_LIBRARIES})
    endif()
endif()

# Symlink MVB testData into our build tree for runtime access
add_custom_command(TARGET mvb_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../MVB/tests/testData
        ${CMAKE_CURRENT_BINARY_DIR}/testData
    COMMENT "Linking MVB testData"
)

