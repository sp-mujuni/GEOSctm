ecbuild_add_executable (
  TARGET GEOSctm.x
  SOURCES GEOSctm.F90
  LIBS GEOSctm_GridComp
  )  

set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

#link_directories (${MKL_LIBRARIES})
target_link_libraries (GEOSctm.x ${OpenMP_Fortran_LIBRARIES})
set_target_properties(GEOSctm.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
target_include_directories (GEOSctm.x PUBLIC ${INC_ESMF})

file (GLOB templates *.tmpl)

set(otheretc
   fvcore_layout.rc
   .GEOSCTM_VERSION
   )

set (programs
   ctm_setup
   ctm_run.j
   ctm_regress.j
   ctm_post.j
   ctm_archive.j
   ctm_convert.j
   ctm_quickplot.csh
   ctm_moveplot.j
   )


install (
   FILES ${templates} ${otheretc}
   DESTINATION etc
   )

install (
   PROGRAMS ${programs}
   DESTINATION bin
   )

install (
   FILES fvcore_layout.rc
   DESTINATION etc
   )

