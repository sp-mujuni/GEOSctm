esma_set_this ()

set (SRCS
  CTM_DiffusionGridCompMod.F90
  CTM_DiffusionStubCompMod.F90
  GmiDiffusionMethod_mod.F90
  updateDiffusion_mod.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSchem_GridComp Chem_Shared MAPL_Base
  )
  
