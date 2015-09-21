# - Try to find CASTOR
#  define CASTOR_DIR to specify own path for search
#  (See http://savannah.cern.ch/files/?group=castor)
#  Check for rfio_api.h, stager_api.h for CASTOR 2 and libshift
#
#  CASTOR_INCLUDE_DIR - where to find rfio_api.h, etc.
#  CASTOR_LIBRARIES   - List of libraries when using ....
#  CASTOR_FOUND       - True if CASTOR 2  libraries found.

set(CASTOR_FOUND FALSE)
set(CASTOR_LIBRARIES)

find_path(CASTOR_INCLUDE_DIR NAMES patchlevel.h PATHS
  ${CASTOR_DIR}/include/shift
  $ENV{CASTOR_DIR}/include/shift
  /cern/pro/include
  /cern/new/include
  /cern/old/include
  /opt/shift/include
  /usr/local/shift/include
  /usr/include/shift
  /usr/local/include/shift 
  /usr/include 
  /usr/local/include
)

if(CASTOR_INCLUDE_DIR)
  file(READ ${CASTOR_INCLUDE_DIR}/patchlevel.h contents)
  string(REGEX MATCH   "BASEVERSION[ ]*[\"][ ]*([^ \"]+)" cont ${contents})
  string(REGEX REPLACE "BASEVERSION[ ]*[\"][ ]*([^ \"]+)" "\\1" CASTOR_VERSION ${cont})
endif()

set(LOCATIONS 
  ${CASTOR_DIR}/lib
  ${CASTOR_DIR}/lib64
  $ENV{CASTOR_DIR}/lib
  $ENV{CASTOR_DIR}/lib64
  /cern/pro/lib 
  /cern/new/lib
  /cern/old/lib 
  /opt/shift/lib 
  /usr/local/shift/lib
  /usr/lib/shift
  /usr/local/lib/shift
  /usr/lib64
  /usr/lib
  /usr/local/lib
)
find_library(CASTOR_shift_LIBRARY NAMES shift shiftmd HINTS ${LOCATIONS})
find_library(CASTOR_rfio_LIBRARY NAMES castorrfio PATHS ${LOCATIONS})
find_library(CASTOR_common_LIBRARY NAMES castorcommon PATHS ${LOCATIONS})
find_library(CASTOR_client_LIBRARY NAMES castorclient castorClient PATHS ${LOCATIONS})
find_library(CASTOR_ns_LIBRARY NAMES castorns PATHS ${LOCATIONS})

if(CASTOR_shift_LIBRARY)
  set(CASTOR_LIBRARIES ${CASTOR_LIBRARIES} ${CASTOR_shift_LIBRARY})
endif()

if(CASTOR_rfio_LIBRARY)
  set(CASTOR_LIBRARIES ${CASTOR_LIBRARIES} ${CASTOR_shift_LIBRARY})
endif()

if(CASTOR_common_LIBRARY)
  set(CASTOR_LIBRARIES ${CASTOR_LIBRARIES} ${CASTOR_shift_LIBRARY})
endif()

if(CASTOR_client_LIBRARY)
  set(CASTOR_LIBRARIES ${CASTOR_LIBRARIES} ${CASTOR_shift_LIBRARY})
endif()

if(CASTOR_ns_LIBRARY)
  set(CASTOR_LIBRARIES ${CASTOR_LIBRARIES} ${CASTOR_shift_LIBRARY})
endif()

if(CASTOR_shift_LIBRARY)
	get_filename_component(CASTOR_LIBRARY_DIR ${CASTOR_shift_LIBRARY} DIRECTORY)
endif()

if(CASTOR_INCLUDE_DIR AND CASTOR_LIBRARIES)
  set(CASTOR_FOUND TRUE)
  message(STATUS "Found Castor version ${CASTOR_VERSION}")
  message(STATUS "Found Castor include dir ${CASTOR_INCLUDE_DIR}")
  message(STATUS "Found Castor library dir ${CASTOR_LIBRARY_DIR}")
else()
  message(STATUS "Castor not found!")
endif()

set(CASTOR_INCLUDE_DIR ${CASTOR_INCLUDE_DIR} ${CASTOR_INCLUDE_DIR}/../)

mark_as_advanced(
  CASTOR_LIBRARIES
  CASTOR_FOUND
  CASTOR_INCLUDE_DIR
)
