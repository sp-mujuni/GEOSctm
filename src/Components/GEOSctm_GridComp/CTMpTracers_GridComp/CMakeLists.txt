esma_set_this ()

if (pTracers_MODE MATCHES "stub")
  set (SRCS
    CTM_pTracersStubCompMod.F90
    )
else ()
  set (SRCS
    CTM_pTracersGridCompMod.F90
    )
endif ()  

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES FVdycoreCubed_GridComp GMIchem_GridComp Chem_Shared MAPL GMAO_hermes
  )
  
