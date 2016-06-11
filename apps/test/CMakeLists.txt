# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/alf/CMakeLists.txt
#
# CMakeLists.txt file for r_test.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_r_test)
message (STATUS "Configuring apps/r_test")

set (SEQAN_APP_VERSION "1.1.1")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (r_test r_test.cpp)
add_executable (r_test1 r_test1.cpp)
add_executable (r_test2 r_test2.cpp)
add_executable (r_test3 r_test3.cpp)


# Add dependencies found by find_package (SeqAn).
target_link_libraries (r_test ${SEQAN_LIBRARIES})
target_link_libraries (r_test1 ${SEQAN_LIBRARIES})
target_link_libraries (r_test2 ${SEQAN_LIBRARIES})
target_link_libraries (r_test3 ${SEQAN_LIBRARIES})


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (r_test)
  seqan_setup_install_vars (r_test1)
  seqan_setup_install_vars (r_test2)
  seqan_setup_install_vars (r_test3)

endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install alf in ${PREFIX}/bin directory
install (TARGETS r_test 
         DESTINATION bin)
install (TARGETS r_test1 
         DESTINATION bin)
install (TARGETS r_test2
        DESTINATION bin)
install (TARGETS r_test3
        DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/alf for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/small.fasta
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_r_test (r_test)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------
	
# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} alf CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:r_test")
  set (CPACK_PACKAGE_NAME "r_test")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "r_test - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (r_test "r_test")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:r_test")

