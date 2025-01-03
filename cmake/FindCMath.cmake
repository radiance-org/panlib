
include(CheckLibraryExists)
include(CheckSymbolExists)

# First try to find math functions without explicit library
check_symbol_exists(pow "math.h" HAVE_POW_NO_LIB)

if(HAVE_POW_NO_LIB)
    # Math functions available without explicit library
    set(CMath_LIBRARY "")
    set(HAVE_POW TRUE)
else()
    # Try to find math functions in libm
    if(UNIX)
        check_library_exists(m pow "" HAVE_POW_IN_LIBM)
        if(HAVE_POW_IN_LIBM)
            set(CMath_LIBRARY "m")
            set(HAVE_POW TRUE)
        endif()
    endif()
endif()

if(HAVE_POW)
    if(NOT TARGET CMath::CMath)
        add_library(CMath::CMath INTERFACE IMPORTED)
        if(CMath_LIBRARY)
            set_target_properties(CMath::CMath PROPERTIES
                INTERFACE_LINK_LIBRARIES "${CMath_LIBRARY}"
            )
        endif()
    endif()
    
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(CMath
        REQUIRED_VARS HAVE_POW
        FAIL_MESSAGE "Could not find math library functions"
    )
    
    mark_as_advanced(CMath_LIBRARY)
else()
    set(CMath_FOUND FALSE)
    if(CMath_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find required math library functions")
    endif()
endif()
