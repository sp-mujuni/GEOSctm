ecbuild_add_executable (
  TARGET GEOSctm.x
  SOURCES GEOSctm.F90
  LIBS GEOSctm_GridComp
  )

set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

target_link_libraries (GEOSctm.x OpenMP::OpenMP_Fortran)
target_include_directories (GEOSctm.x PUBLIC ${INC_ESMF})

file (GLOB templates CONFIGURE_DEPENDS *.tmpl)

set (programs
   ctm_run.j
   ctm_regress.j
   ctm_post.j
   ctm_archive.j
   ctm_convert.j
   ctm_quickplot.csh
   ctm_moveplot.j
   )

install (
   FILES ${templates}
   DESTINATION etc
   )

install (
   PROGRAMS ${programs}
   DESTINATION bin
   )

install (
   FILES fvcore_layout.rc logging.yaml
   DESTINATION etc
   )

if(HYDROSTATIC)
   set(CFG_HYDROSTATIC TRUE)
else()
   set(CFG_HYDROSTATIC FALSE)
endif()

# Did we build for AMD Rome hardware (aka EPYC)?
cmake_host_system_information(RESULT proc_decription QUERY PROCESSOR_DESCRIPTION)
if (${proc_decription} MATCHES "EPYC")
  set(CFG_BUILT_ON_ROME TRUE)
else ()
  set(CFG_BUILT_ON_ROME FALSE)
endif ()

set (setup_scripts
   ctm_setup
   )

foreach (file ${setup_scripts})
   configure_file(${file} ${file} @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${file} DESTINATION bin)
endforeach ()


configure_file(.GEOSCTM_VERSION .GEOSCTM_VERSION @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/.GEOSCTM_VERSION DESTINATION etc)
