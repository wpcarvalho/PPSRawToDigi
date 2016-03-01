#use TOTEM_RAW_DATA_EXTRA_INCLUDES to pass include paths from outside project (for example if -DMONITOR is set add ../monitor/interface)
#use TOTEM_RAW_DATA_EXTRA_SOURCES to pass src paths from outside project (for example if -DMONITOR is set add ../monitor/src)

cmake_minimum_required (VERSION 2.8)
project (TotemRawDataLibrary)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)
set(GENERATED_SOURCES_DIR ${CMAKE_BINARY_DIR}/generated)

# defining rpath to avoid clashes with other version on LD_LIBRARY_PATH
# and avoid adjusting LD_LIBRARY_PATH manually
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

## -- -- -- -- --   SETTING UP LIBRARIES    -- -- -- -- -- ##

#CASTOR
include(FindCastor.cmake)
IF(CASTOR_FOUND)
    ADD_DEFINITIONS(-DUSE_CASTOR)
    INCLUDE_DIRECTORIES(${CASTOR_INCLUDE_DIR})
	LINK_DIRECTORIES(${CASTOR_LIBRARY_DIR})
ENDIF(CASTOR_FOUND)

#XRootD
include(FindXRootD.cmake)
IF(XROOTD_FOUND)
    ADD_DEFINITIONS(-DUSE_XROOTD)
    INCLUDE_DIRECTORIES(${XROOTD_INCLUDES})
    LINK_DIRECTORIES(${XROOTD_LIB_DIR})
ENDIF(XROOTD_FOUND)

## -- -- -- -- --   SETTING UP ENVIRONMENT    -- -- -- -- -- ##

# setting up build env properties
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -pthread -g")
set(CMAKE_CXX_COMPILER gcc)

#switch using -DCMAKE_BUILD_TYPE=[Debug|Release]
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXbtX_FLAGS_RELEASE} -O3 -g0")


## -- -- -- -- --   CREATING LIBRARIES        -- -- -- -- -- ##
## -- -- -- -- --   ONE PER SOURCE FOLDER     -- -- -- -- -- ##
# structuring sources into something alike packages
# and building it one-library-per-package is done only for convenience (listing all sources for each binary is tedious)
# and lets see which binary depends on which packages

# as those libraries structure shouldn't be exposed we don't want link them at runtime
set(BUILD_SHARED_LIBS FALSE)

INCLUDE_DIRECTORIES(Readers/interface)
FILE(GLOB TotemRawDataLibrary_Readers_src "Readers/src/*.cc")

if(NOT XROOTD_FOUND)
    list(REMOVE_ITEM TotemRawDataLibrary_Readers_src "${CMAKE_CURRENT_LIST_DIR}/Readers/src/XRootStorageFile.cc")
endif(NOT XROOTD_FOUND)

if(NOT CASTOR_FOUND)
    list(REMOVE_ITEM TotemRawDataLibrary_Readers_src "${CMAKE_CURRENT_LIST_DIR}/Readers/src/RFIOStorageFile.cc")
endif(NOT CASTOR_FOUND)

add_library (TotemRawDataLibrary_Readers_lib ${TotemRawDataLibrary_Readers_src})


INCLUDE_DIRECTORIES(DataFormats/interface)
FILE(GLOB TotemRawDataLibrary_DataFormats_src "DataFormats/src/*.cc")
add_library (TotemRawDataLibrary_DataFormats_lib ${TotemRawDataLibrary_DataFormats_src})

INCLUDE_DIRECTORIES(Utilities/interface)
FILE(GLOB TotemRawDataLibrary_Utilities_src "Utilities/src/*.cc")
add_library (TotemRawDataLibrary_Utilities_lib ${TotemRawDataLibrary_Utilities_src})


## -- -- -- -- --   CREATING EXECUTABLES      -- -- -- -- -- ##
INCLUDE_DIRECTORIES(..)
INCLUDE_DIRECTORIES(${TOTEM_RAW_DATA_EXTRA_INCLUDES})

#analyzeFile (Utilities) exec
ADD_EXECUTABLE(analyzeFile Utilities/bin/analyzeFile.cc ${TOTEM_RAW_DATA_EXTRA_SOURCES})
TARGET_LINK_LIBRARIES(analyzeFile
    stdc++
    TotemRawDataLibrary_Readers_lib TotemRawDataLibrary_DataFormats_lib
)

#checkFile (Utilities) exec
ADD_EXECUTABLE(checkFile Utilities/bin/checkFile.cc ${TOTEM_RAW_DATA_EXTRA_SOURCES})
TARGET_LINK_LIBRARIES(checkFile
    stdc++
    TotemRawDataLibrary_Utilities_lib TotemRawDataLibrary_Readers_lib
    TotemRawDataLibrary_DataFormats_lib
)

#readFile (Utilities) exec
ADD_EXECUTABLE(readFile Utilities/bin/readFile.cc ${TOTEM_RAW_DATA_EXTRA_SOURCES})
TARGET_LINK_LIBRARIES(readFile
    stdc++
    TotemRawDataLibrary_Utilities_lib TotemRawDataLibrary_Readers_lib
    TotemRawDataLibrary_DataFormats_lib
	m
)

#checkCounterMajority (Utilities) exec
ADD_EXECUTABLE(checkCounterMajority Utilities/bin/checkCounterMajority.cc ${TOTEM_RAW_DATA_EXTRA_SOURCES})
TARGET_LINK_LIBRARIES(checkCounterMajority
    stdc++
    TotemRawDataLibrary_Utilities_lib TotemRawDataLibrary_Readers_lib
    TotemRawDataLibrary_DataFormats_lib
	m
)

IF(CASTOR_FOUND)
    TARGET_LINK_LIBRARIES(analyzeFile shift castorrfio)
    TARGET_LINK_LIBRARIES(checkFile shift castorrfio)
    TARGET_LINK_LIBRARIES(readFile shift castorrfio)
    TARGET_LINK_LIBRARIES(checkCounterMajority shift castorrfio)
ENDIF(CASTOR_FOUND)

if(XROOTD_FOUND)
    TARGET_LINK_LIBRARIES(analyzeFile XrdPosix XrdPosixPreload) #XrdClient
    TARGET_LINK_LIBRARIES(checkFile XrdPosix XrdPosixPreload)
    TARGET_LINK_LIBRARIES(readFile XrdPosix XrdPosixPreload)
    TARGET_LINK_LIBRARIES(checkCounterMajority XrdPosix XrdPosixPreload)
endif(XROOTD_FOUND)




## -- -- -- -- --   Testing libraries      -- -- -- -- -- ##

# Google C++ Mock
include(FindGMock.cmake)

# Google C++ Test Framework
# set GTEST_ROOT - (as a CMake or environment variable)
# or rely on GMock is shipped with GTest bundled
find_package(GTest)
if(GMOCK_FOUND AND NOT GTEST_FOUND)
    message(STATUS "Looking for GTest bundled with GMock...")
    set(GTEST_ROOT ${GMOCK_INCLUDE_DIRS}/../gtest)
    find_package(GTest)
endif()



## -- -- -- -- --   Test cases addition     -- -- -- -- -- ##
if(XROOTD_FOUND AND CASTOR_FOUND)
    add_executable(StorageFile_test Readers/test/StorageFileTest.cc)
    target_link_libraries(StorageFile_test TotemRawDataLibrary_Readers_lib
        XrdPosix
        XrdPosixPreload
        XrdClient
        castorrfio
        shift
        stdc++)
endif(XROOTD_FOUND AND CASTOR_FOUND)

IF(GTEST_FOUND AND GMOCK_FOUND AND NOT CASTOR_FOUND)

    enable_testing()
    include_directories(${GMOCK_INCLUDE_DIRS})
    include_directories(${GTEST_INCLUDE_DIRS})

    add_executable(OptoRxSupplementalData_test DataFormats/test/OptoRxSupplementalData.cc)
    target_link_libraries(OptoRxSupplementalData_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES} )
    add_test(DataFormats_OptoRxSupplementalData ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/OptoRxSupplementalData_test)

    add_executable(RawVFATFrame_test DataFormats/test/RawVFATFrame.cc)
    target_link_libraries(RawVFATFrame_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES}  TotemRawDataLibrary_DataFormats_lib)
    add_test(DataFormats_RawVFATFrame ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/RawVFATFrame_test)
    
    add_executable(ClusterizationVFATFrame_test DataFormats/test/ClusterizationVFATFrame.cc)
    target_link_libraries(ClusterizationVFATFrame_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES}  TotemRawDataLibrary_DataFormats_lib)
    add_test(DataFormats_ClusterizationVFATFrame ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ClusterizationVFATFrame_test)
    
    add_executable(PositionedVFATFrameCollection_test DataFormats/test/PositionedVFATFrameCollection.cc)
    target_link_libraries(PositionedVFATFrameCollection_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES}  TotemRawDataLibrary_DataFormats_lib)
    add_test(DataFormats_PositionedVFATFrameCollection ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/PositionedVFATFrameCollection_test)
    
    add_executable(DynamicOptoRxVMEBFile_test Readers/test/DynamicOptoRxVMEBFile.cc)
    target_link_libraries(DynamicOptoRxVMEBFile_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES} TotemRawDataLibrary_Readers_lib TotemRawDataLibrary_DataFormats_lib m)
    add_test(Readers_DynamicOptoRxVMEBFile ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/DynamicOptoRxVMEBFile_test)

    add_executable(Integration_DynamicOptoRxVMEBFile_test Readers/test/Integration_DynamicOptoRxVMEBFile.cc)
    target_link_libraries(Integration_DynamicOptoRxVMEBFile_test stdc++ ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES} TotemRawDataLibrary_Readers_lib TotemRawDataLibrary_DataFormats_lib m)
    add_test(Readers_Integration_DynamicOptoRxVMEBFile ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Integration_DynamicOptoRxVMEBFile_test)
ELSE()
    IF(CASTOR_FOUND)
        MESSAGE(STATUS "Using castor - tests not included (castor clashes with gtest)")
    ELSE()
        MESSAGE(STATUS "No GTest/GMock - tests not included")
    ENDIF()
ENDIF()
