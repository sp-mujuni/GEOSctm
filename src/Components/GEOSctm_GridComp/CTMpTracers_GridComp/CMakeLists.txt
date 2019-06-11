esma_set_this ()

set (SRCS
  CTM_pTracersGridCompMod.F90
  CTM_pTracersStubCompMod.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base
  )
  
