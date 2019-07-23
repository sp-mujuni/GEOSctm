esma_set_this ()

set (SRCS
   GEOSctm.F90
   )

esma_add_library (${this}
   SRCS ${SRCS}
   DEPENDENCIES GEOSctm_GridComp MAPL_Base
   )

ecbuild_add_executable (
   TARGET GEOSctm.x
   SOURCES GEOSctm.F90
   LIBS ${this}
   )  
