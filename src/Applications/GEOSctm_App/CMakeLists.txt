esma_set_this ()

set (SRCS
  GEOSctm.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSctm_GridComp MAPL_Base
  )
  
