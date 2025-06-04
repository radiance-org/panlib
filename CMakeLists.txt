cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


cmake_policy(SET CMP0077 NEW)

# Build type
set(CMAKE_BUILD_TYPE Release)

include(ExternalProject)

if(MSVC)
    # add_definitions(/MP)
    add_definitions(-D_HAS_STD_BYTE=0)
    add_definitions(-DNOMINMAX)
    add_compile_options(/wd"4244" /wd"4305")
endif()
if(WIN32 AND "${CMAKE_C_COMPILER_ID}" MATCHES "^(Intel)$")
  set(_INTEL_WINDOWS 1)
endif()
if(MSVC OR _INTEL_WINDOWS)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
endif()
if(WIN32)
  add_definitions(-DEZXML_NOMMAP -Dfseeko=fseek -Dstrcasecmp=_stricmp)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(rtrad)

# Configure OpenEXR (which includes Imath)
set(OPENEXR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openexr)
set(OPENEXR_IS_SUBPROJECT TRUE)
set(OPENEXR_INSTALL ON)
set(OPENEXR_INSTALL_TOOLS OFF)
set(OPENEXR_BUILD_LIBS TRUE)
set(OPENEXR_BUILD_EXAMPLES FALSE)
set(OPENEXR_BUILD_TOOLS FALSE)
# Let OpenEXR handle Imath configuration
set(OPENEXR_FORCE_INTERNAL_IMATH ON)
add_subdirectory(${OPENEXR_DIR})

# Set include directories for OpenEXR and Imath
set(OPENEXR_SOURCE_DIR ${OPENEXR_DIR}/src/lib)
set(OPENEXR_INCLUDE_DIR ${OPENEXR_SOURCE_DIR}/OpenEXR ${OPENEXR_SOURCE_DIR}/Iex ${CMAKE_CURRENT_BINARY_DIR}/external/openexr/cmake)
set(IMATH_INCLUDE_DIR ${OPENEXR_DIR}/src/lib/Imath ${CMAKE_CURRENT_BINARY_DIR}/external/openexr/src/lib/Imath)

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
set(LIBJPEG_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libjpeg-turbo-build)
set(LIBJPEG_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/libjpeg-turbo-install)

# Set the correct library name based on platform
if(WIN32)
    set(JPEG_LIBRARY_NAME "jpeg-static")
    set(LIB_SUFFIX "")
else()
    set(JPEG_LIBRARY_NAME "jpeg")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS "/etc/redhat-release")
        set(LIB_SUFFIX "64")
    else()
        set(LIB_SUFFIX "")
    endif()
endif()

ExternalProject_Add(libjpeg-turbo
    SOURCE_DIR ${LIBJPEG_SOURCE_DIR}
    BINARY_DIR ${LIBJPEG_BINARY_DIR}
    INSTALL_DIR ${LIBJPEG_INSTALL_DIR}
    CMAKE_ARGS
        -DENABLE_STATIC=ON
        -DENABLE_SHARED=OFF
        -DCMAKE_INSTALL_PREFIX=${LIBJPEG_INSTALL_DIR}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
)

# Create imported target for libjpeg
add_library(jpeg STATIC IMPORTED GLOBAL)
# We need to delay setting the properties until after the build
file(MAKE_DIRECTORY ${LIBJPEG_INSTALL_DIR}/include)  # Create include directory

ExternalProject_Get_Property(libjpeg-turbo install_dir)
set_target_properties(jpeg PROPERTIES
    IMPORTED_LOCATION_RELEASE ${install_dir}/lib${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}${JPEG_LIBRARY_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}
    IMPORTED_LOCATION ${install_dir}/lib${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}${JPEG_LIBRARY_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_INCLUDE_DIRECTORIES ${install_dir}/include
)
add_dependencies(jpeg libjpeg-turbo)


set(LIBJPEG_INCLUDE_DIR ${BINARY_DIR})

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
add_dependencies(pan libjpeg-turbo)
target_link_libraries(pan PRIVATE cpprad rtrad OpenEXR tiff jpeg)
target_include_directories(pan PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${RAD_COMMON_DIR} ${LIBJPEG_INCLUDE_DIR} ${LIBJPEG_SOURCE_DIR}/src ${LIBJPEG_BINARY_DIR} ${OPENEXR_INCLUDE_DIR} ${IMATH_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
