cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0077 NEW)

include(FetchContent)
include(ExternalProject)

# Set global build options
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries" FORCE)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


FetchContent_Declare(
    radiance
    GIT_REPOSITORY https://github.com/lbnl-eta/radiance
    GIT_TAG master
)

# Populate the content without adding the entire project
FetchContent_GetProperties(radiance)
if(NOT radiance_POPULATED)
    FetchContent_Populate(radiance)
endif()

# Get the path to the common directory
set(RADIANCE_COMMON_DIR ${radiance_SOURCE_DIR}/src/common)

# OpenEXR
ExternalProject_Add(openexr_external
    GIT_REPOSITORY https://github.com/mitsuba-renderer/openexr
    GIT_TAG snapshot_2
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
	-DBUILD_TESTING=OFF
)

ExternalProject_Add(libjpeg
    GIT_REPOSITORY https://github.com/libjpeg-turbo/libjpeg-turbo
    GIT_TAG 3.1.0
    CMAKE_ARGS
	-DENABLE_SHARED=FALSE
	-DENABLE_STATIC=TRUE
	-DWITH_TURBOJPEG=FALSE
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
	-DBUILD_TESTING=OFF
)



# Configure libtiff options
set(BUILD_SHARED_LIBS OFF)
set(tiff-tools OFF)
set(tiff-tests OFF)
set(tiff-docs OFF)
set(tiff-contrib OFF)
set(tiff-install OFF)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} ${MATH_LIBRARY}")
endif()

set(zlib FALSE CACHE BOOL "" FORCE)
set(pixarlog FALSE CACHE BOOL "" FORCE)
# set(ZLIB_SUPPORT FALSE CACHE BOOL "Enable Zlib usage" FORCE)
set(JPEG_SUPPORT OFF CACHE BOOL "Enable JPEG compression" FORCE)
set(OJPEG_SUPPORT OFF CACHE BOOL "Enable old-style JPEG compression" FORCE)
set(JBIG_SUPPORT OFF CACHE BOOL "Enable JBIG compression" FORCE)
set(WEBP_SUPPORT OFF CACHE BOOL "Enable WebP compression" FORCE)
set(ZSTD_SUPPORT OFF CACHE BOOL "Enable ZSTD compression" FORCE)
set(DEFLATE_SUPPORT OFF CACHE BOOL "Enable DEFLATE compression" FORCE)

# Fetch and configure libtiff
FetchContent_Declare(
    libtiff
    GIT_REPOSITORY https://gitlab.com/libtiff/libtiff.git
    GIT_TAG v4.7.0
)
FetchContent_MakeAvailable(libtiff)

# Set up include directories
ExternalProject_Get_Property(libjpeg INSTALL_DIR)
set(LIBJPEG_INCLUDE_DIR ${INSTALL_DIR}/include)
set(LIBJPEG_LIBRARY_DIR ${INSTALL_DIR}/lib)
set(LIBJPEG_INCLUDE_DIR2 ${INSTALL_DIR}/src/libjpeg-build)
ExternalProject_Get_Property(libjpeg SOURCE_DIR)
set(LIBJPEG_INCLUDE_DIR3 ${SOURCE_DIR}/src)

ExternalProject_Get_Property(openexr_external INSTALL_DIR)
set(OPENEXR_INCLUDE_DIR ${INSTALL_DIR}/include/OpenEXR)
set(OPENEXR_LIBRARY_DIR ${INSTALL_DIR}/lib)

# ExternalProject_Get_Property(libtiff_external INSTALL_DIR)
set(TIFF_INCLUDE_DIR ${INSTALL_DIR}/include)
# set(TIFF_LIBRARY_DIR ${INSTALL_DIR}/lib)
# Verify the directory exists

set(RADIANCE_COMMON_DIR ${radiance_SOURCE_DIR}/src/common)
if(NOT EXISTS ${RADIANCE_COMMON_DIR})
    message(FATAL_ERROR "Radiance common directory not found at ${RADIANCE_COMMON_DIR}")
endif()


# Set source files
set(SOURCES
	src/cache.cpp
	src/dbaccess.cpp
	src/dbheader.cpp
	src/dbhlink.cpp
	src/dbrecord.cpp
	src/exif.cpp
	src/exifthumb.cpp
	src/jsfidctflt.cpp
	src/jstreamsrc.cpp
	src/panimage.cpp
	src/panwriter.cpp
	src/pdbase.cpp
	src/pdequant.cpp
	src/pdispimg.cpp
	src/pdraw.cpp
	src/pfeatures.cpp
	src/phdalign.cpp
	src/phdflare.cpp
	src/phdrimg.cpp
	src/photophile.cpp
	src/picache.cpp
	src/ppano.cpp
	src/pstrings.cpp
	src/pthumb.cpp
	src/readexr.cpp
	src/readjpeg.cpp
	src/readtiff.cpp
	src/textform.cpp
	src/textmap.cpp
	src/thumbicon.cpp
	src/tiffin.cpp
	src/win_platform.cpp
	src/writeexr.cpp
	src/writejpeg.cpp
)


# Create library
add_library(pan ${SOURCES})

# Include directories
target_include_directories(pan
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(pan PRIVATE jpeg-turbo tiff openexr m)
target_link_directories(pan PRIVATE ${JPEG_TURBO_LIBRARY_DIR} ${OPENEXR_LIBRARY_DIR} ${TIFF_LIBRARY_DIR})
target_include_directories(pan PRIVATE ${RADIANCE_COMMON_DIR} ${libtiff_SOURCE_DIR}/libtiff ${LIBJPEG_INCLUDE_DIR} ${LIBJPEG_INCLUDE_DIR2} ${LIBJPEG_INCLUDE_DIR3} ${OPENEXR_INCLUDE_DIR})

# Installation
install(TARGETS pan
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
