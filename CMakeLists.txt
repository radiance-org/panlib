cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0077 NEW)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(ExternalProject)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(rtrad)

set(IMATH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Imath)
set(IMATH_IS_SUBPROJECT TRUE)
add_subdirectory(${IMATH_DIR})
set(IMATH_INCLUDE_DIR ${IMATH_DIR}/src/Imath ${CMAKE_CURRENT_BINARY_DIR}/external/Imath/config)

set(OPENEXR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openexr)
set(OPENEXR_IS_SUBPROJECT TRUE)
set(OPENEXR_INSTALL ON)
set(OPENEXR_INSTALL_TOOLS OFF)
set(OPENEXR_BUILD_LIBS TRUE)
set(OPENEXR_BUILD_EXAMPLES FALSE)
set(OPENEXR_BUILD_TOOLS FALSE)
add_subdirectory(${OPENEXR_DIR})
set(OPENEXR_SOURCE_DIR ${OPENEXR_DIR}/src/lib)
set(OPENEXR_INCLUDE_DIR ${OPENEXR_SOURCE_DIR}/OpenEXR ${OPENEXR_SOURCE_DIR}/Iex ${CMAKE_CURRENT_BINARY_DIR}/external/openexr/cmake)


set(tiff-tools OFF CACHE BOOL "" FORCE)
set(tiff-tests OFF CACHE BOOL "" FORCE)
set(tiff-docs OFF CACHE BOOL "" FORCE)
set(tiff-contrib OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(zlib OFF CACHE BOOL "" FORCE)
set(jpeg OFF CACHE BOOL "" FORCE)
set(pixarlog OFF CACHE BOOL "" FORCE)
set(libdeflate OFF CACHE BOOL "" FORCE) 
set(lerc OFF CACHE BOOL "" FORCE)
set(lzma OFF CACHE BOOL "" FORCE)
set(webp OFF CACHE BOOL "" FORCE)
set(zstd OFF CACHE BOOL "" FORCE)
set(CMAKE_CROSSCOMPILING OFF CACHE BOOL "" FORCE)
add_subdirectory(external/libtiff)


set(LIBJPEG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libjpeg-turbo)
ExternalProject_Add(libjpeg-turbo
SOURCE_DIR ${LIBJPEG_SOURCE_DIR}
CMAKE_ARGS
	-DENABLE_STATIC=TRUE
	-DENABLE_SHARED=FALSE
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

ExternalProject_Get_Property(libjpeg-turbo INSTALL_DIR)
set(LIBJPEG_INCLUDE_DIR ${INSTALL_DIR}/src/libjpeg-turbo-build ${LIBJPEG_SOURCE_DIR}/src)
# Create an internal imported target for jpeg
add_library(jpeg STATIC IMPORTED)
# add_dependencies(jpeg libjpeg-turbo)

set(LIBJPEG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libjpeg-turbo)
set_target_properties(jpeg PROPERTIES
    IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}jpeg${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/src/libjpeg-turbo-build
)

set(RAD_COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/radiance/src/common)

# Set source files
set(SOURCES
src/dmessage.c
src/jhcomp.c
src/jhdecomp.c
src/jhresamp.c
src/jhtonemap2.c
src/jhtonemap.c
src/memobject.c
src/pblendpano.c
src/picolor.c
src/pimage.c
src/pireader.c
src/piresamp.c
src/readbmp.c
src/readrad.c
src/rmdirectory.c
src/system.c
src/tiffmsg.c
src/writerad.c
src/writetiff.c
src/cache.cpp
src/dbaccess.cpp
src/dbheader.cpp
src/dbhlink.cpp
src/dbrecord.cpp
src/exif.cpp
src/jsfidctflt.cpp
src/jstreamsrc.cpp
src/panimage.cpp
src/panwriter.cpp
src/pdbase.cpp
src/pdispimg.cpp
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
src/thumbicon.cpp
src/tiffin.cpp
src/writeexr.cpp
src/writejpeg.cpp
src/piwarp.c
src/pisum.c
src/pidequant.c
src/pbilat.c
src/phisto.c
src/phistomatch.c
src/pconvolve.c
src/exifthumb.cpp
src/pdilate.c
src/pfeatures.cpp
src/radheader.c
src/readnrm.c
src/readdpt.c
src/writenrm.c
src/writedpt.c
src/pdequant.cpp
src/pdraw.cpp
src/pmipmap.c
src/readmtx.c
src/writemtx.c
src/textmap.cpp
src/phistadj.c
)



# Create library
add_library(pan STATIC ${SOURCES})
target_link_libraries(pan PRIVATE cpprad rtrad OpenEXR tiff jpeg)
target_include_directories(pan PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${RAD_COMMON_DIR} ${LIBJPEG_INCLUDE_DIR} ${OPENEXR_INCLUDE_DIR} ${IMATH_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
