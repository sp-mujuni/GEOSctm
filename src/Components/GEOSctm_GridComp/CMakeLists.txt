esma_set_this ()

set (alldirs
  @GEOSgcm_GridComp
  @GEOSchem_GridComp
  CTMconvection_GridComp
  CTMdiffusion_GridComp
  CTMpTracers_GridComp
  )

set (FV_PRECISION R8)
esma_add_library (${this}
  SRCS GEOS_ctmEnvGridComp.F90 GEOS_ctmGridCompMod.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES FVdycoreCubed_Gridcomp MAPL_Base)
