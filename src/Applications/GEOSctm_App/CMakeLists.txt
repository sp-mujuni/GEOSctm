esma_set_this ()

set (SRCS
  GEOSctm.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base
  )
  
