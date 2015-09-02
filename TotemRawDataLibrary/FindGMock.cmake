# Locate the Google C++ Mocking Framework.
#
# Defines the following variables:
#
#   GMOCK_FOUND - Found the Google Mocking framework
#   GMOCK_INCLUDE_DIRS - Include directories
#
# Also defines the library variables below as normal
# variables.  These contain debug/optimized keywords when
# a debugging library is found.
#
#   GMOCK_BOTH_LIBRARIES - Both libgmock & libgmock-main
#   GMOCK_LIBRARIES - libgmock
#   GMOCK_MAIN_LIBRARIES - libgmock-main
#
# Accepts the following variables as input:
#
#   GMOCK_ROOT - (as CMake or env. variable)
#                The root directory of the gmock install prefix
#

find_path(GMOCK_INCLUDE_DIRS NAMES gmock/gmock.h HINTS
  ${GMOCK_ROOT}/include
  $ENV{GMOCK_ROOT}/include
)

find_library (GMOCK_LIBRARIES NAMES gmock HINTS ${GMOCK_ROOT} $ENV{GMOCK_ROOT})
find_library (GMOCK_MAIN_LIBRARIES NAMES gmock_main HINTS ${GMOCK_ROOT} $ENV{GMOCK_ROOT})

if(GMOCK_INCLUDE_DIRS AND GMOCK_LIBRARIES AND GMOCK_MAIN_LIBRARIES)
   set(GMOCK_FOUND TRUE)
   set(GMOCK_BOTH_LIBRARIES "${GMOCK_LIBRARIES};${GMOCK_MAIN_LIBRARIES}")
   message(STATUS "Found GMock: ${GMOCK_LIBRARIES}")
else()
   set(GMOCK_FOUND FALSE)
   message(STATUS "GMock not found!")
endif()

mark_as_advanced(
    GMOCK_INCLUDE_DIRS
    GMOCK_BOTH_LIBRARIES
    GMOCK_LIBRARIES
    GMOCK_MAIN_LIBRARIES
    GMOCK_FOUND
)
