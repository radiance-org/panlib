cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

include(ExternalProject)
# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0077 NEW)


ExternalProject_Add(libjpeg-turbo 
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libjpeg-turbo
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    CMAKE_ARGS
        -DENABLE_SHARED=FALSE
        -DENABLE_STATIC=TRUE
        -DBUILD_TESTING=OFF
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Add(openexr
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openexr
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DBUILD_TESTING=OFF
        -DOPENEXR_BUILD_BOTH_STATIC_SHARED=ON
        -DPYILMBASE_ENABLE=OFF
)

ExternalProject_Get_Property(libjpeg-turbo INSTALL_DIR)
ExternalProject_Get_Property(libjpeg-turbo SOURCE_DIR)
set(LIBJPEG_INCLUDE_DIR ${INSTALL_DIR}/src/libjpeg-turbo-build ${SOURCE_DIR}/src)

ExternalProject_Get_Property(openexr INSTALL_DIR)
set(OPENEXR_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openexr)
set(OPENEXR_INCLUDE_DIR 
${OPENEXR_SOURCE_DIR}/OpenEXR/IlmImf 
${OPENEXR_SOURCE_DIR}/IlmBase/Imath 
${OPENEXR_SOURCE_DIR}/IlmBase/Iex 
${OPENEXR_SOURCE_DIR}/IlmBase/Half
${INSTALL_DIR}/src/openexr-build/OpenEXR/config 
${INSTALL_DIR}/src/openexr-build/IlmBase/config
)

set(RAD_COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/radiance/src/common)


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
target_include_directories(pan PRIVATE ${CMAKE_SOURCE_DIR}/src ${RAD_COMMON_DIR} ${LIBJPEG_INCLUDE_DIR} ${OPENEXR_INCLUDE_DIR})
