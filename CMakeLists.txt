cmake_minimum_required(VERSION 3.16)
project(ARCEDB VERSION 1.0 LANGUAGES CXX)
# make details
# set (CMAKE_VERBOSE_MAKEFILE on) 
set(CMAKE_CXX_STANDARD 20)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DARCEDB_DEBUG)
    message(STATUS "STRAT ARCEDB DEBUG MODE")
    set(CMAKE_CXX_FLAGS "-O0 -g")
else() 
    set(CMAKE_CXX_FLAGS "-march=native -O3 -funroll-loops")
ENDIF()

#--------------------------Source tree-----------------------------
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")


#------------------------------Install-----------------------------
set(ARCEDB_INCLUDES_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/ARCEDB)

# ------Add source files to library and header files to install----
set(ARCEDB_SOURCE_FILES "")
add_subdirectory(src/ARCEDB/utils)
add_subdirectory(src/ARCEDB/conversion)
add_subdirectory(src/ARCEDB/comparison)
# -------------------------Add thirdparty--------------------------
add_subdirectory(thirdparty/SEAL)
add_subdirectory(thirdparty/TFHEpp)

find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#-------------------Build only a static library--------------------
add_library(ARCEDB STATIC ${ARCEDB_SOURCE_FILES})
target_include_directories(ARCEDB PUBLIC src)
target_include_directories(ARCEDB PUBLIC thirdparty/SEAL/natice/src)

# find_package(NTL)
# IF(${NTL_FOUND})
#     target_include_directories(ARCEDB PUIBLIC ${CMAKE_INSTALL_PREFIX}/include)
#     target_link_libraries(ARCEDB ${NTL_LIBRARIES})
# endif(NTL_FOUND)

target_link_libraries(ARCEDB tfhe++ seal ntl)
add_subdirectory(test)