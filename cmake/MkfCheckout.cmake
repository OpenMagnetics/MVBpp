# The MKF checkout under _deps/mkf-src carries two edits that configure makes to it: the
# CMakeLists SHARED->STATIC (and CAS dependency) patch, and a symlink to the separately fetched MAS
# in place of MKF's uninitialised MAS submodule. FetchContent's git update stashes local edits
# before moving to the new main and pops them after, and the pop fails on every MKF main advance
# ("CONFLICT (distinct types): MAS" -- a symlink against a gitlink), which stopped native and WASM
# configures alike.
#
# mvbpp_mkf_unplant() is called BEFORE FetchContent_Populate(MKF): it takes those two edits back
# out, so the update step sees a clean checkout, and configure plants them again afterwards. Any
# OTHER change in the checkout stops the configure with the list, rather than being stashed,
# reverted or built: it is not ours to discard, and a build of an edited MKF records the wrong
# revision.
function(mvbpp_mkf_unplant src)
    if(NOT EXISTS "${src}/.git")
        return()  # first configure: nothing fetched yet
    endif()
    if(IS_SYMLINK "${src}/MAS")
        file(REMOVE "${src}/MAS")          # removes the link, never the MAS checkout it points to
        file(MAKE_DIRECTORY "${src}/MAS")  # what git leaves for an uninitialised submodule
    endif()
    execute_process(COMMAND git -C "${src}" checkout -- CMakeLists.txt RESULT_VARIABLE _rc)
    if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "mvbpp_mkf_unplant: cannot restore ${src}/CMakeLists.txt from git")
    endif()
    execute_process(COMMAND git -C "${src}" status --porcelain --untracked-files=no
                    OUTPUT_VARIABLE _st OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE _rc)
    if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "mvbpp_mkf_unplant: git status failed in ${src}")
    endif()
    if(NOT _st STREQUAL "")
        message(FATAL_ERROR "mvbpp_mkf_unplant: ${src} has changes this build did not make; "
                            "the MKF update would stash them. Resolve them by hand:\n${_st}")
    endif()
endfunction()
