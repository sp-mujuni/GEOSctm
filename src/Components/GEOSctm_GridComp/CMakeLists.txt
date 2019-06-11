esma_set_this ()

set (alldirs
  CTMconvection_GridComp
  CTMdiffusion_GridComp
  CTMpTracers_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_ctmEnvGridComp.F90 GEOS_ctmGridCompMod.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL_Base)
