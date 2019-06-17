esma_set_this ()

set (SRCS
  CTM_pTracersGridCompMod.F90
  CTM_pTracersStubCompMod.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES FVdycoreCubed_GridComp GMIchem_GridComp Chem_Shared MAPL_Base GMAO_hermes
  )
  
