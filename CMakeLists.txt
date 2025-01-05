cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

include(ExternalProject)
# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0077 NEW)


set(IMATH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Imath)
add_subdirectory(${IMATH_DIR})
set(IMATH_INCLUDE_DIR ${IMATH_DIR}/src/Imath ${CMAKE_CURRENT_BINARY_DIR}/external/Imath/config)

set(OPENEXR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openexr)
set(OPENEXR_IS_SUBPROJECT TRUE)
set(OPENEXR_INSTALL ON)
set(OPENEXR_INSTALL_TOOLS OFF)
set(OPENEXR_BUILD_LIBS FALSE)
set(OPENEXR_BUILD_EXAMPLES FALSE)
set(OPENEXR_BUILD_TOOLS FALSE)
add_subdirectory(${OPENEXR_DIR})
set(OPENEXR_SOURCE_DIR ${OPENEXR_DIR}/src/lib)
set(OPENEXR_INCLUDE_DIR ${OPENEXR_SOURCE_DIR}/OpenEXR ${OPENEXR_SOURCE_DIR}/Iex ${CMAKE_CURRENT_BINARY_DIR}/external/openexr/cmake)

if(WITH_JPEG8)
  set(JPEG_LIB_VERSION 80)
elseif(WITH_JPEG7)
  set(JPEG_LIB_VERSION 70)
else()
  set(JPEG_LIB_VERSION 62)
endif()
math(EXPR JPEG_LIB_VERSION_DIV10 "${JPEG_LIB_VERSION} / 10")
math(EXPR JPEG_LIB_VERSION_MOD10 "${JPEG_LIB_VERSION} % 10")
if(JPEG_LIB_VERSION STREQUAL "62")
  set(DEFAULT_SO_MAJOR_VERSION ${JPEG_LIB_VERSION})
else()
  set(DEFAULT_SO_MAJOR_VERSION ${JPEG_LIB_VERSION_DIV10})
endif()
if(JPEG_LIB_VERSION STREQUAL "80")
  set(DEFAULT_SO_MINOR_VERSION 2)
else()
  set(DEFAULT_SO_MINOR_VERSION 0)
endif()

set(JPEG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libjpeg-turbo/src)
configure_file(${JPEG_SOURCE_DIR}/jversion.h.in jversion.h)
configure_file(${JPEG_SOURCE_DIR}/jconfig.h.in jconfig.h)
configure_file(${JPEG_SOURCE_DIR}/jconfigint.h.in jconfigint.h)


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
target_include_directories(pan PRIVATE ${CMAKE_SOURCE_DIR}/src ${RAD_COMMON_DIR} ${JPEG_SOURCE_DIR} ${OPENEXR_INCLUDE_DIR} ${IMATH_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
