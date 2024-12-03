cmake_minimum_required(VERSION 3.10)
project(panimage VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options for format support
option(WITH_JPEG "Enable JPEG support" ON)
option(WITH_TIFF "Enable TIFF support" ON)
option(WITH_OPENEXR "Enable OpenEXR support" ON)

# Find required packages based on options
if(WITH_JPEG)
    find_package(JPEG REQUIRED)
endif()

if(WITH_TIFF)
    find_package(TIFF REQUIRED)
endif()

if(WITH_OPENEXR)
    find_package(OpenEXR REQUIRED)
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
	src/pfeatures_dbg.cpp
	src/pfeatures_v1.cpp
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

# Conditionally exclude format-specific sources
if(NOT WITH_JPEG)
    list(REMOVE_ITEM SOURCES 
        src/readjpeg.cpp
        src/writejpeg.cpp
        src/jstreamsrc.cpp
        src/jsfidctflt.cpp
    )
endif()

if(NOT WITH_TIFF)
    list(REMOVE_ITEM SOURCES 
        src/readtiff.cpp
        src/writetiff.cpp
        src/tiffin.cpp
    )
endif()

if(NOT WITH_OPENEXR)
    list(REMOVE_ITEM SOURCES 
        src/readexr.cpp
        src/writeexr.cpp
    )
endif()

# Create library
add_library(pan ${SOURCES})

# Include directories
target_include_directories(pan
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# Set compile definitions and link libraries based on options
if(WITH_JPEG)
    target_compile_definitions(pan PRIVATE HAVE_JPEG_SUPPORT)
    target_link_libraries(pan PRIVATE JPEG::JPEG)
endif()

if(WITH_TIFF)
    target_compile_definitions(pan PRIVATE HAVE_TIFF_SUPPORT)
    target_link_libraries(pan PRIVATE TIFF::TIFF)
endif()

if(WITH_OPENEXR)
    target_compile_definitions(pan PRIVATE HAVE_OPENEXR_SUPPORT)
    target_link_libraries(pan PRIVATE OpenEXR::OpenEXR)
endif()

# Installation
install(TARGETS pan
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
