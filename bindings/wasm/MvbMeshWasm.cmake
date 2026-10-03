# mvbpp_mesh for WASM (ABT #1588): the browser meshes with the same gmsh fork and MMG as the
# native build, compiled by Emscripten against this module's own OCCT (occt_wasm), so gmsh's OCC
# kernel and MVB++ share one OCCT inside the one .wasm.
#
# Differences from native, all forced by the platform rather than chosen:
#   - no OpenMP (single-threaded wasm): the mesh recipe pins threads = 1 anyway, which is what
#     makes meshes repeatable (OMFEM measured 6 threads non-repeatable on 3/3 designs);
#   - no GMP (gmsh's Kbipack homology only; not used by meshing).
# MMG's build runs its own `genheader` tool to write Fortran headers; under Emscripten that tool
# is a .js file the build cannot execute, so a host-compiled genheader is used instead (the same
# patch the existing mmg-wasm build carries).

set(MVBPP_MMG_WASM_PREFIX "${CMAKE_BINARY_DIR}/mmg-install")
set(MVBPP_GMSH_WASM_PREFIX "${CMAKE_BINARY_DIR}/gmsh-install")
set(_mmg_genheader "${CMAKE_BINARY_DIR}/mmg-genheader-host")

ExternalProject_Add(
    mmg_wasm
    GIT_REPOSITORY https://github.com/MmgTools/mmg.git
    GIT_TAG v5.8.0
    GIT_SHALLOW TRUE
    # Host genheader + the macro that calls it instead of the cross-compiled one.
    PATCH_COMMAND cc <SOURCE_DIR>/scripts/genheader.c -o ${_mmg_genheader}
        COMMAND ${CMAKE_COMMAND} -DSRC=<SOURCE_DIR> -P ${CMAKE_CURRENT_SOURCE_DIR}/mmg_wasm_patch.cmake
    CMAKE_ARGS
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_SHARED_LIBS=OFF
        -DLIBMMG_STATIC=ON
        -DLIBMMG3D_STATIC=ON
        -DBUILD_TESTING=OFF
        -DMMG_PATTERN=OFF
        -DUSE_ELAS=OFF
        -DUSE_POINTMAP=OFF
        -DUSE_SCOTCH=OFF
        -DUSE_VTK=OFF
        -DNATIVE_GENHEADER=${_mmg_genheader}
        # Emscripten's libc carries libm: there is no separate math library for MMG to find.
        -DM_LIB=
        -DCMAKE_INSTALL_PREFIX=${MVBPP_MMG_WASM_PREFIX}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> -j 3
    BUILD_BYPRODUCTS ${MVBPP_MMG_WASM_PREFIX}/lib/libmmg.a
    USES_TERMINAL_BUILD TRUE
)

ExternalProject_Add(
    gmsh_wasm
    DEPENDS occt_wasm mmg_wasm
    GIT_REPOSITORY https://github.com/AlfVII/gmsh.git
    GIT_TAG om-4.15.2
    GIT_SHALLOW TRUE
    CMAKE_ARGS
        -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
        -DCMAKE_BUILD_TYPE=Release
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
        -DENABLE_OPENMP=OFF
        -DENABLE_GMP=OFF
        -DENABLE_ALGLIB=ON -DENABLE_ANN=ON -DENABLE_BAMG=ON -DENABLE_BLOSSOM=ON
        -DENABLE_DINTEGRATION=ON -DENABLE_DOMHEX=ON -DENABLE_EIGEN=ON -DENABLE_GMM=ON
        -DENABLE_HXT=ON -DENABLE_KBIPACK=ON -DENABLE_MATHEX=ON -DENABLE_MESH=ON -DENABLE_METIS=OFF
        -DENABLE_MMG=ON -DENABLE_NETGEN=ON -DENABLE_NII2MESH=ON -DENABLE_OCC=ON
        -DENABLE_OCC_CAF=ON -DENABLE_ONELAB=ON -DENABLE_ONELAB_METAMODEL=ON
        -DENABLE_OPTHOM=ON -DENABLE_PARSER=ON -DENABLE_PLUGINS=ON -DENABLE_POST=ON
        -DENABLE_QUADMESHINGTOOLS=ON -DENABLE_QUADTRI=ON -DENABLE_SOLVER=ON -DENABLE_TINYXML2=ON
        -DENABLE_UNTANGLE=ON -DENABLE_VOROPP=ON -DENABLE_WINSLOWUNTANGLER=ON
        -DCMAKE_PREFIX_PATH=${OCCT_PREFIX}
        -DCMAKE_FIND_ROOT_PATH=${OCCT_PREFIX}
        -DOCC_INC=${OCCT_PREFIX}/include/opencascade
        -DMMG_INC=${MVBPP_MMG_WASM_PREFIX}/include
        -DMMG_LIB=${MVBPP_MMG_WASM_PREFIX}/lib/libmmg.a
        -DCMAKE_INSTALL_PREFIX=${MVBPP_GMSH_WASM_PREFIX}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> -j 3 --target lib
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>
    BUILD_BYPRODUCTS ${MVBPP_GMSH_WASM_PREFIX}/lib/libgmsh.a
    USES_TERMINAL_BUILD TRUE
)

set(_mvbpp_root "${CMAKE_CURRENT_SOURCE_DIR}/../..")
# What this module is compiled from, refreshed on every build (see ../../cmake/BuildRev.cmake).
set(_mvbpp_gen "${CMAKE_BINARY_DIR}/generated")
add_custom_target(mvbpp_buildrev
    COMMAND ${CMAKE_COMMAND} -DOUT=${_mvbpp_gen}/mvb/mesh/BuildRev.h
            -DMVBPP_DIR=${_mvbpp_root} -DMKF_DIR=${_MKF_SRC} -DMAS_DIR=${_MAS_SRC}
            -P ${_mvbpp_root}/cmake/BuildRev.cmake
    BYPRODUCTS ${_mvbpp_gen}/mvb/mesh/BuildRev.h
    COMMENT "Recording the MVB++/MKF/MAS revisions mvbpp_mesh is built from")
add_library(mvbpp_mesh STATIC
    ${_mvbpp_root}/src/mesh/Mesher.cpp
    ${_mvbpp_root}/src/mesh/SizeField.cpp
    ${_mvbpp_root}/src/mesh/MeshRecipe.cpp
    ${_mvbpp_root}/src/mesh/MeshIO.cpp
    ${_mvbpp_root}/src/mesh/MeshParts.cpp
    ${_mvbpp_root}/src/mesh/MeshSupport.cpp
    ${_mvbpp_root}/src/mesh/Mesh3d.cpp
)
add_dependencies(mvbpp_mesh gmsh_wasm mvbpp_buildrev)
target_compile_options(mvbpp_mesh PRIVATE -fexceptions)
target_include_directories(mvbpp_mesh
    PUBLIC  ${_mvbpp_root}/include
    PRIVATE ${MVBPP_GMSH_WASM_PREFIX}/include ${MVBPP_MMG_WASM_PREFIX}/include ${_mvbpp_gen})
target_link_libraries(mvbpp_mesh
    PUBLIC  mvb++ nlohmann_json::nlohmann_json
    PRIVATE ${MVBPP_GMSH_WASM_PREFIX}/lib/libgmsh.a ${MVBPP_MMG_WASM_PREFIX}/lib/libmmg.a occt)
