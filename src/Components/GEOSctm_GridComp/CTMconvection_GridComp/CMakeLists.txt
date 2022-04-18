 esma_set_this ()

if (Convection_MODE MATCHES stub)
  set (SRCS
    convectiveTransport_mod.F90
    CTM_rasCalculationsMod.F90
    )
else ()
  set (SRCS
    CTM_ConvectionGridCompMod.F90
    GmiConvectionMethod_mod.F90
    GenericConvectionMethod_mod.F90
    convectiveTransport_mod.F90
    CTM_rasCalculationsMod.F90
    )
endif ()  

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSchem_GridComp Chem_Shared MAPL
  )
  
