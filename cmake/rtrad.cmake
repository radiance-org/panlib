set(RADIANCE_COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/radiance/src/common)

set(rtrad_SOURCES
  ${RADIANCE_COMMON_DIR}/addobjnotify.c
  ${RADIANCE_COMMON_DIR}/badarg.c
  ${RADIANCE_COMMON_DIR}/biggerlib.c
  ${RADIANCE_COMMON_DIR}/bmalloc.c
  ${RADIANCE_COMMON_DIR}/bmpfile.c
  ${RADIANCE_COMMON_DIR}/bsdf.c
  ${RADIANCE_COMMON_DIR}/bsdf_m.c
  ${RADIANCE_COMMON_DIR}/bsdf_t.c
  ${RADIANCE_COMMON_DIR}/byteswap.c
  ${RADIANCE_COMMON_DIR}/caldefn.c
  ${RADIANCE_COMMON_DIR}/calexpr.c
  ${RADIANCE_COMMON_DIR}/calfunc.c
  ${RADIANCE_COMMON_DIR}/calprnt.c
  ${RADIANCE_COMMON_DIR}/ccolor.c
  ${RADIANCE_COMMON_DIR}/ccyrgb.c
  ${RADIANCE_COMMON_DIR}/chanvalue.c
  ${RADIANCE_COMMON_DIR}/clip.c
  ${RADIANCE_COMMON_DIR}/color.c
  ${RADIANCE_COMMON_DIR}/colrops.c
  ${RADIANCE_COMMON_DIR}/cone.c
  ${RADIANCE_COMMON_DIR}/cvtcmd.c
  ${RADIANCE_COMMON_DIR}/data.c
  ${RADIANCE_COMMON_DIR}/depthcodec.c
  ${RADIANCE_COMMON_DIR}/dircode.c
  ${RADIANCE_COMMON_DIR}/disk2square.c
  ${RADIANCE_COMMON_DIR}/dmessage.c
  ${RADIANCE_COMMON_DIR}/ealloc.c
  ${RADIANCE_COMMON_DIR}/eputs.c
  ${RADIANCE_COMMON_DIR}/erf.c
  ${RADIANCE_COMMON_DIR}/error.c
  ${RADIANCE_COMMON_DIR}/expandarg.c
  ${RADIANCE_COMMON_DIR}/ezxml.c
  ${RADIANCE_COMMON_DIR}/face.c
  ${RADIANCE_COMMON_DIR}/falsecolor.c
  ${RADIANCE_COMMON_DIR}/fdate.c
  ${RADIANCE_COMMON_DIR}/fgetline.c
  ${RADIANCE_COMMON_DIR}/fgetval.c
  ${RADIANCE_COMMON_DIR}/fgetword.c
  ${RADIANCE_COMMON_DIR}/fixargv0.c
  ${RADIANCE_COMMON_DIR}/fltdepth.c
  ${RADIANCE_COMMON_DIR}/font.c
  ${RADIANCE_COMMON_DIR}/fputword.c
  ${RADIANCE_COMMON_DIR}/free_os.c
  ${RADIANCE_COMMON_DIR}/fropen.c
  ${RADIANCE_COMMON_DIR}/fvect.c
  ${RADIANCE_COMMON_DIR}/gethomedir.c
  ${RADIANCE_COMMON_DIR}/getlibpath.c
  ${RADIANCE_COMMON_DIR}/getpath.c
  ${RADIANCE_COMMON_DIR}/header.c
  ${RADIANCE_COMMON_DIR}/hilbert.c
  ${RADIANCE_COMMON_DIR}/idmap.c
  ${RADIANCE_COMMON_DIR}/image.c
  ${RADIANCE_COMMON_DIR}/instance.c
  ${RADIANCE_COMMON_DIR}/interp2d.c
  ${RADIANCE_COMMON_DIR}/invmat4.c
  ${RADIANCE_COMMON_DIR}/jitteraperture.c
  ${RADIANCE_COMMON_DIR}/lamps.c
  ${RADIANCE_COMMON_DIR}/linregr.c
  ${RADIANCE_COMMON_DIR}/loadbsdf.c
  ${RADIANCE_COMMON_DIR}/loadvars.c
  ${RADIANCE_COMMON_DIR}/lookup.c
  ${RADIANCE_COMMON_DIR}/mat4.c
  ${RADIANCE_COMMON_DIR}/mesh.c
  ${RADIANCE_COMMON_DIR}/modobject.c
  ${RADIANCE_COMMON_DIR}/multisamp.c
  ${RADIANCE_COMMON_DIR}/myhostname.c
  ${RADIANCE_COMMON_DIR}/normcodec.c
  ${RADIANCE_COMMON_DIR}/objset.c
  ${RADIANCE_COMMON_DIR}/octree.c
  ${RADIANCE_COMMON_DIR}/otypes.c
  ${RADIANCE_COMMON_DIR}/paths.c
  ${RADIANCE_COMMON_DIR}/plocate.c
  ${RADIANCE_COMMON_DIR}/portio.c
  ${RADIANCE_COMMON_DIR}/process.c
  ${RADIANCE_COMMON_DIR}/quit.c
  ${RADIANCE_COMMON_DIR}/readfargs.c
  ${RADIANCE_COMMON_DIR}/readmesh.c
  ${RADIANCE_COMMON_DIR}/readobj.c
  ${RADIANCE_COMMON_DIR}/readoct.c
  ${RADIANCE_COMMON_DIR}/resolu.c
  ${RADIANCE_COMMON_DIR}/rexpr.c
  ${RADIANCE_COMMON_DIR}/savestr.c
  ${RADIANCE_COMMON_DIR}/savqstr.c
  ${RADIANCE_COMMON_DIR}/sceneio.c
  ${RADIANCE_COMMON_DIR}/spec_rgb.c
  ${RADIANCE_COMMON_DIR}/tcos.c
  ${RADIANCE_COMMON_DIR}/timegm.c
  ${RADIANCE_COMMON_DIR}/tmap16bit.c
  ${RADIANCE_COMMON_DIR}/tmapcolrs.c
  ${RADIANCE_COMMON_DIR}/tmapluv.c
  ${RADIANCE_COMMON_DIR}/tmaptiff.c
  ${RADIANCE_COMMON_DIR}/tmesh.c
  ${RADIANCE_COMMON_DIR}/tonemap.c
  ${RADIANCE_COMMON_DIR}/triangulate.c
  ${RADIANCE_COMMON_DIR}/urand.c
  ${RADIANCE_COMMON_DIR}/urind.c
  ${RADIANCE_COMMON_DIR}/wordfile.c
  ${RADIANCE_COMMON_DIR}/words.c
  ${RADIANCE_COMMON_DIR}/wputs.c
  ${RADIANCE_COMMON_DIR}/xf.c
  ${RADIANCE_COMMON_DIR}/zeroes.c
)

if(UNIX)
	list(APPEND rtrad_SOURCES ${RADIANCE_COMMON_DIR}/unix_process.c)
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	list(APPEND rtrad_SOURCES ${RADIANCE_COMMON_DIR}/strnstr.c ${RADIANCE_COMMON_DIR}/strlcpy.c)
  endif()
else()
	list(APPEND rtrad_SOURCES ${RADIANCE_COMMON_DIR}/win_process.c ${RADIANCE_COMMON_DIR}/win_popen.c ${RADIANCE_COMMON_DIR}/win_usleep.c ${RADIANCE_COMMON_DIR}/strnstr.c ${RADIANCE_COMMON_DIR}/strlcpy.c)
endif()

add_library(rtrad STATIC ${rtrad_SOURCES})

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
	target_link_libraries(rtrad m)
endif()

if(WIN32)
  target_link_libraries(rtrad ws2_32)
endif()

add_library(cpprad STATIC
  ${RADIANCE_COMMON_DIR}/abitmap.cpp
  ${RADIANCE_COMMON_DIR}/abitmapio.cpp
)
