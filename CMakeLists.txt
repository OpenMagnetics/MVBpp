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
# OCCT install location. Deliberately OUTSIDE the build tree by default: OCCT takes ~20 minutes
# to compile, and putting it under build/ meant (a) any build-tree wipe destroyed it, and (b) every
# clean build produced a FRESH OCCT binary -- so the geometry kernel underneath the whole pipeline
# was not reproducible from one build to the next. A shared, built-once install makes OCCT a fixed
# input rather than a per-build variable. Override with -DOCCT_PREFIX=... to isolate a build.
set(OCCT_PREFIX "$ENV{HOME}/OpenMagnetics/occt-7.9.3-install"
    CACHE PATH "OCCT 7.9.3 install prefix (shared across build trees; built once)")

# The toolkits MVB++ links. This list is the SINGLE SOURCE OF TRUTH: the ExternalProject's
# BUILD_BYPRODUCTS below is derived from it, so the two can never drift apart. They HAD drifted --
# BUILD_BYPRODUCTS named 27 toolkits while this list links 46 -- which made a FROM-SCRATCH Ninja
# build impossible: ninja refuses a link input that neither exists nor has a rule ("libTKV3d.a,
# needed by libmvb++.so, missing and no known rule to make it"). It only ever worked because a
# previously-built occt-install was already on disk; wiping the build tree exposed it.
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
set(OCCT_BYPRODUCTS "")
foreach(_l ${OCCT_ESSENTIAL_LIBS})
    list(APPEND OCCT_BYPRODUCTS "${OCCT_PREFIX}/lib/lib${_l}.a")
endforeach()

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
    BUILD_BYPRODUCTS ${OCCT_BYPRODUCTS}
    USES_TERMINAL_BUILD TRUE
)

add_library(occt INTERFACE)
add_dependencies(occt occt_external)
target_include_directories(occt INTERFACE ${OCCT_PREFIX}/include/opencascade)
target_link_directories(occt INTERFACE ${OCCT_PREFIX}/lib)

# OCCT 7.9.3 toolkits have circular dependencies (e.g. TKGeomAlgo <-> TKMath),
# which no single left-to-right ordering of static archives can satisfy — GNU ld
# would drop symbols pulled later (undefined references to math_FunctionSet,
# GProp_PGProps, ...). Wrap the whole set in --start-group/--end-group so ld
# iterates the archives until all cross-references resolve. Use full .a paths so
# the group survives CMake's link-item processing.
set(_occt_group_libs)
foreach(_lib ${OCCT_ESSENTIAL_LIBS})
    list(APPEND _occt_group_libs "${OCCT_PREFIX}/lib/lib${_lib}.a")
endforeach()
target_link_libraries(occt INTERFACE
    -Wl,--start-group
    ${_occt_group_libs}
    -Wl,--end-group
)

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

# Both repos are fetched from GitHub: the commits MVB++ needs (MAS febef0f with
# the bobbin base/divider schema, MKF 6b23d51a with Bobbin::get_pin_rails and
# Coil::terminal_exit_slots) are pushed, so the build no longer has to read the
# developer's local clones -- which may be dirty or mid-rebase in a shared box.
set(_MAS_LOCAL "https://github.com/OpenMagnetics/MAS.git")
set(_MKF_LOCAL "https://github.com/OpenMagnetics/MKF.git")

FetchContent_Declare(
    MAS
    GIT_REPOSITORY "${_MAS_LOCAL}"
    # Track MAS main. Provides the magneticApplication / subApplication schema
    # that MKF main (below) requires. Moves in lockstep with the MKF main pin.
    GIT_TAG        main
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
    # Track MKF main (in lockstep with the MAS-main pin above). The header-only
    # get_wire_insulation_relative_permittivity that OMFEM's electrostatic mesh needs is
    # now upstreamed on main (StrayCapacitance.h), so the old 028b3b9c pin is no longer
    # needed; pinning it against MAS-main caused a quicktype schema-skew codegen failure.
    GIT_TAG        main
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
# Initialize submodules MKF needs (cci_coords for codegen, CAS for CAS.hpp /
# Topology.h, PEAS for the PEAS-integrated schemas). MKF origin/main no longer
# carries an EAS submodule, so initialising it would fail the whole command and
# leave CAS uninitialised — list only the submodules that actually exist.
if(IS_DIRECTORY "${_MKF_SRC}/cci_coords" AND NOT EXISTS "${_MKF_SRC}/cci_coords/coordinates")
    execute_process(
        COMMAND git submodule update --init cci_coords CAS PEAS
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
    src/PinBuilder.cpp
    src/BaseBuilder.cpp
    src/ShuntBuilder.cpp
    src/FR4Builder.cpp
    src/StepExporter.cpp
    src/Symmetry.cpp
    src/Utils.cpp
    src/TurnBuilder.cpp
    src/ConductorBuilder.cpp
    src/TerminalFillet.cpp
    # ABT #685: the ONE centreline -> copper assembler. Real-winding Coil geometry reaches
    # OCCT only through this TU (Alf's rule); the six chunk primitives emit centrelines.
    src/WireAssembler.cpp
    src/BobbinBuilder.cpp
    src/shapes/ShapeBuilder.cpp
    src/shapes/ShapeDrum.cpp
    src/shapes/ShapeDrumRing.cpp
    src/shapes/ShapeMolded.cpp
    src/shapes/ShapeDrumSemishielded.cpp
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
    tests/test_core_coating.cpp
    tests/test_symmetry.cpp
    tests/test_mas_battery.cpp
    tests/test_toroidal_topview.cpp
    tests/test_toroidal_single_turn_centered.cpp
    tests/test_oblong_planar_turn.cpp
    tests/test_real_winding.cpp
    tests/test_terminal_lead_lengths.cpp
    tests/test_toroid_mounting.cpp
    tests/test_section_projection.cpp
    tests/test_spacers.cpp
    tests/test_pins.cpp
    tests/test_pinrail.cpp
    tests/test_shunts.cpp
    tests/test_named_roles.cpp
    tests/test_conducting_default.cpp
    tests/test_enamel_gate_verdict.cpp
    tests/test_toroid_base.cpp
    tests/test_step_roundtrip.cpp
    tests/test_sleeves.cpp
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

add_executable(mvbpp_fitprobe tools/mvbpp_fitprobe.cpp)
target_link_libraries(mvbpp_fitprobe PRIVATE mvb++)

# ABT #970 step 1: the foil meshing spike (bands only, no primitive in the product path).
add_executable(mvbpp_foilspike tools/mvbpp_foilspike.cpp)
target_link_libraries(mvbpp_foilspike PRIVATE mvb++)

add_executable(mvbpp_pathprobe tools/mvbpp_pathprobe.cpp)
target_link_libraries(mvbpp_pathprobe PRIVATE mvb++)

add_executable(mvbpp_paint2d tools/mvbpp_paint2d.cpp)
target_link_libraries(mvbpp_paint2d PRIVATE mvb++)

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
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/tests
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/realwinding_fixtures
        ${CMAKE_CURRENT_BINARY_DIR}/tests/realwinding_fixtures
    COMMENT "Linking MVB testData and MVB++ realwinding fixtures"
)

