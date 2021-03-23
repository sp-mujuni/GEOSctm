esma_set_this ()

if (Diffusion_MODE MATCHES "stub")
  set (SRCS
    CTM_DiffusionGridStubCompMod.F90
    GmiDiffusionMethod_mod.F90
    updateDiffusion_mod.F90
    )
else ()
  set (SRCS
    CTM_DiffusionGridCompMod.F90
    GmiDiffusionMethod_mod.F90
    updateDiffusion_mod.F90
    )
endif ()  



esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSchem_GridComp Chem_Shared MAPL
  )
  
