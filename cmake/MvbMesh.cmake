# mvbpp_mesh: MAS -> FEM mesh (ABT #1588). gmsh + MMG, built ONCE against MVB++'s own OCCT so the
# geometry kernel under the mesher is the one that drew the geometry (one OCCT per process: gmsh's
# OCC kernel and MVB++ share it, which is what lets the WASM module carry both).
#
# Like OCCT, both install OUTSIDE the build tree by default: a build-tree wipe must not rebuild
# the mesher underneath a reference corpus, and the identity gate compares meshes made by one
# fixed gmsh. Override the prefixes to isolate a build.
#
# gmsh comes from our fork (github.com/AlfVII/gmsh, branch om-4.15.2: upstream 4.15.2 plus the
# OMFEM patches -- Steiner re-homing, pyramid apices inside the tet region, the O(1) OCC _unbind,
# the RAII generator lock). Its options and compile flags are those of the shared gmsh OMFEM's
# reference meshes were made with (gmsh-om415-build), checked object by object: every ENABLE_*
# the same (METIS off, gmsh's default is on) and -fvisibility=hidden as gmsh gives its shared
# target. A module switched on or off changes which algorithms exist, and so the mesh.

set(MVBPP_MMG_PREFIX "$ENV{HOME}/OpenMagnetics/mmg-5.8.0-install"
    CACHE PATH "MMG install prefix (shared across build trees; built once)")
set(MVBPP_GMSH_PREFIX "$ENV{HOME}/OpenMagnetics/gmsh-om-4.15.2-install"
    CACHE PATH "gmsh (OpenMagnetics fork) install prefix (shared across build trees; built once)")

# An install that already exists is reused, never rebuilt (same reason as OCCT's).
if(EXISTS "${MVBPP_MMG_PREFIX}/lib/libmmg.a")
    message(STATUS "MMG: reusing the install at ${MVBPP_MMG_PREFIX}")
    add_custom_target(mmg_external)
else()
ExternalProject_Add(
    mmg_external
    GIT_REPOSITORY https://github.com/MmgTools/mmg.git
    GIT_TAG v5.8.0
    GIT_SHALLOW TRUE
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=OFF
        -DLIBMMG_STATIC=ON
        -DLIBMMG3D_STATIC=ON
        -DBUILD_TESTING=OFF
        -DMMG_PATTERN=OFF
        -DUSE_ELAS=OFF
        -DUSE_POINTMAP=OFF
        -DUSE_SCOTCH=OFF
        -DUSE_VTK=OFF
        -DCMAKE_INSTALL_PREFIX=${MVBPP_MMG_PREFIX}
    # The shared box builds at -j3 (memory: the machine was killed under pressure, 2026-09-12).
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> -j 3
    BUILD_BYPRODUCTS ${MVBPP_MMG_PREFIX}/lib/libmmg.a
    USES_TERMINAL_BUILD TRUE
)
endif()

if(EXISTS "${MVBPP_GMSH_PREFIX}/lib/libgmsh.a")
    message(STATUS "gmsh: reusing the install at ${MVBPP_GMSH_PREFIX}")
    add_custom_target(gmsh_external)
else()
ExternalProject_Add(
    gmsh_external
    DEPENDS occt_external mmg_external
    GIT_REPOSITORY https://github.com/AlfVII/gmsh.git
    GIT_TAG om-4.15.2
    GIT_SHALLOW TRUE
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DENABLE_BUILD_LIB=ON
        # gmsh compiles only its SHARED target with -fvisibility=hidden. OMFEM's reference meshes
        # were made with that shared library, and without the flag GCC inlines differently in
        # 460 of gmsh's 624 objects; with it, the static objects disassemble to the shared ones.
        "-DCMAKE_C_FLAGS=-fvisibility=hidden"
        "-DCMAKE_CXX_FLAGS=-fvisibility=hidden"
        -DENABLE_BUILD_SHARED=OFF
        -DENABLE_BUILD_DYNAMIC=OFF
        -DENABLE_FLTK=OFF
        -DENABLE_TESTS=OFF
        -DENABLE_BLAS_LAPACK=OFF
        -DENABLE_ALGLIB=ON -DENABLE_ANN=ON -DENABLE_BAMG=ON -DENABLE_BLOSSOM=ON
        -DENABLE_DINTEGRATION=ON -DENABLE_DOMHEX=ON -DENABLE_EIGEN=ON -DENABLE_GMM=ON
        -DENABLE_GMP=ON -DENABLE_HXT=ON -DENABLE_KBIPACK=ON -DENABLE_MATHEX=ON -DENABLE_MESH=ON -DENABLE_METIS=OFF
        -DENABLE_MMG=ON -DENABLE_NETGEN=ON -DENABLE_NII2MESH=ON -DENABLE_OCC=ON
        -DENABLE_OCC_CAF=ON -DENABLE_ONELAB=ON -DENABLE_ONELAB_METAMODEL=ON -DENABLE_OPENMP=ON
        -DENABLE_OPTHOM=ON -DENABLE_PARSER=ON -DENABLE_PLUGINS=ON -DENABLE_POST=ON
        -DENABLE_QUADMESHINGTOOLS=ON -DENABLE_QUADTRI=ON -DENABLE_SOLVER=ON -DENABLE_TINYXML2=ON
        -DENABLE_UNTANGLE=ON -DENABLE_VOROPP=ON -DENABLE_WINSLOWUNTANGLER=ON
        -DCMAKE_PREFIX_PATH=${OCCT_PREFIX}
        -DOCC_INC=${OCCT_PREFIX}/include/opencascade
        -DMMG_INC=${MVBPP_MMG_PREFIX}/include
        -DMMG_LIB=${MVBPP_MMG_PREFIX}/lib/libmmg.a
        -DCMAKE_INSTALL_PREFIX=${MVBPP_GMSH_PREFIX}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> -j 3
    BUILD_BYPRODUCTS ${MVBPP_GMSH_PREFIX}/lib/libgmsh.a
    USES_TERMINAL_BUILD TRUE
)
endif()

find_package(OpenMP REQUIRED)
find_library(MVBPP_GMP_LIBRARY gmp REQUIRED)

# What this library is compiled from, refreshed on every build (cmake/BuildRev.cmake), so a
# mesh's provenance names the code that was LINKED rather than whatever a checkout holds now.
set(_mvbpp_gen "${CMAKE_BINARY_DIR}/generated")
add_custom_target(mvbpp_buildrev
    COMMAND ${CMAKE_COMMAND} -DOUT=${_mvbpp_gen}/mvb/mesh/BuildRev.h
            -DMVBPP_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DMKF_DIR=${_MKF_SRC} -DMAS_DIR=${_MAS_SRC}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildRev.cmake
    BYPRODUCTS ${_mvbpp_gen}/mvb/mesh/BuildRev.h
    COMMENT "Recording the MVB++/MKF/MAS revisions mvbpp_mesh is built from")

add_library(mvbpp_mesh ${_mvbpp_lib_kind}
    src/mesh/Mesher.cpp
    src/mesh/SizeField.cpp
    src/mesh/MeshRecipe.cpp
    src/mesh/MeshIO.cpp
)
add_dependencies(mvbpp_mesh gmsh_external mvbpp_buildrev)
target_include_directories(mvbpp_mesh
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${MVBPP_GMSH_PREFIX}/include ${MVBPP_MMG_PREFIX}/include ${_mvbpp_gen})
target_link_libraries(mvbpp_mesh
    PUBLIC  mvb++ nlohmann_json::nlohmann_json
    PRIVATE ${MVBPP_GMSH_PREFIX}/lib/libgmsh.a ${MVBPP_MMG_PREFIX}/lib/libmmg.a
            occt ${MVBPP_GMP_LIBRARY} OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})

# The round-trip gate on real meshes (tools/mvbpp_mesh_roundtrip.cpp).
add_executable(mvbpp_mesh_roundtrip tools/mvbpp_mesh_roundtrip.cpp)
target_link_libraries(mvbpp_mesh_roundtrip PRIVATE mvbpp_mesh)
