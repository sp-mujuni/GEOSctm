esma_set_this ()

set (SRCS
  CTM_ConvectionGridCompMod.F90
  CTM_ConvectionStubCompMod.F90
  GenericConvectionMethod_mod.F90
  GmiConvectionMethod_mod.F90
  convectiveTransport_mod.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSchem_GridComp Chem_Shared MAPL_Base
  )
  
