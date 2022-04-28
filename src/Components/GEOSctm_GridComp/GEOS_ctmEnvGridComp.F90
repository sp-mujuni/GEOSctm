#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOS_ctmEnvGridComp -- Prepares derived variables for GEOSctm
!
! !INTERFACE:
!
      module GEOS_ctmEnvGridComp
!
! !USES:
      use ESMF
      use MAPL
      use FV_StateMod, only : calcCourantNumberMassFlux => fv_computeMassFluxes
      use fv_arrays_mod , only: FVPRC
      use m_set_eta,  only : set_eta
      use GEOS_FV3_UtilitiesMod, only: A2D2c
      use CTM_rasCalculationsMod, only: INIT_RASPARAMS, DO_RAS, RASPARAM_Type

      implicit none
      private

! !PUBLIC MEMBER FUNCTIONS:

      public SetServices
      public compAreaWeightedAverage

      interface compAreaWeightedAverage
         module procedure compAreaWeightedAverage_2d
         module procedure compAreaWeightedAverage_3d
      end interface
!
! !DESCRIPTION:
! This GC is used to derive variables needed by the CTM GC children.

      !Derived types for the internal state
      type T_CTMenv_STATE
         private
         logical                    :: enable_pTracers     = .FALSE. ! do idealized Passive Tracers?
         logical                    :: read_advCoreFields  = .FALSE. ! read courant numbers and mass fluxes?
         logical                    :: read_PLE            = .FALSE. ! read PLE instead of PS?
         logical                    :: enable_rasCalculations = .FALSE. ! do RAS calculations?
         logical                    :: do_ctmAdvection        = .TRUE.  ! do Advection?
         TYPE(RASPARAM_Type)        :: RASPARAMS
         character(len=ESMF_MAXSTR) :: metType                       ! MERRA2 or MERRA1 or FPIT or FP
      end type T_CTMenv_STATE

      type CTMenv_WRAP
         type (T_CTMenv_STATE), pointer :: PTR
      end type CTMenv_WRAP
!
! !AUTHORS:
! Jules.Kouatchou-1@nasa.gov
!
!EOP
!-------------------------------------------------------------------------
      integer,  parameter :: r8     = SELECTED_REAL_KIND(14,300)
      integer,  parameter :: r4     = SELECTED_REAL_KIND(6,30)

      real(r8), parameter :: RADIUS = MAPL_RADIUS
      real(r8), parameter :: PI     = MAPL_PI_R8
      real(r8), parameter :: GPKG   = 1000.0d0
      real(r8), parameter :: MWTAIR =   28.96d0
      real(r8), parameter :: SecondsPerMinute = 60.0d0

!-------------------------------------------------------------------------
      CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices -- Sets ESMF services for this component
!
! !INTERFACE:
!
      subroutine SetServices ( GC, RC )
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
!
! !OUTPUT PARAMETERS:
      integer, intent(OUT)               :: RC  ! return code
!
! !DESCRIPTION:  
!   The SetServices for the CTM Der GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs. 
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer                    :: STATUS
      type (ESMF_Config)         :: CF
      type (ESMF_Config)         :: configFile
      character(len=ESMF_MAXSTR) :: COMP_NAME
      CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen = 'CTM_GridComp.rc'
      character(len=ESMF_MAXSTR) :: IAm = 'SetServices'
      type (T_CTMenv_STATE), pointer :: state
      type (CTMenv_WRAP)             :: wrap

     ! Get my name and set-up traceback handle
     ! ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
      Iam = trim(COMP_NAME) // TRIM(Iam)

      ! Wrap internal state for storing in GC; rename legacyState
      ! -------------------------------------
      allocate ( state, stat=STATUS )
      VERIFY_(STATUS)
      wrap%ptr => state

     ! Register services for this component
     ! ------------------------------------
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__ )
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        __RC__ )

      ! Save pointer to the wrapped internal state in the GC
      !-----------------------------------------------------

      call ESMF_UserCompSetInternalState ( GC, 'CTMenv', wrap, STATUS )
      VERIFY_(STATUS)

      configFile = ESMF_ConfigCreate(__RC__)

      call ESMF_ConfigLoadFile(configFile, TRIM(rcfilen), __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%enable_pTracers,        &
                             Default  = .FALSE.,                       &
                             Label    = "ENABLE_pTracers:",     __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%read_PLE,     &
                             Default  = .FALSE.,                       &
                             Label    = "read_PLE:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%read_advCoreFields,     &
                             Default  = .FALSE.,                       &
                             Label    = "read_advCoreFields:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%enable_rasCalculations,     &
                             Default  = .FALSE.,                       &
                             Label    = "enable_rasCalculations:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%do_ctmAdvection,       &
                                     Default  = .TRUE.,                     &
                                     Label    = "do_ctmAdvection:",  __RC__ )

      ! Type of meteological fields (MERRA2 or MERRA1 or FPIT or FP)
      call ESMF_ConfigGetAttribute(configFile, state%metType,             &
                                     Default  = 'MERRA2',           &
                                     Label    = "metType:",  __RC__ )

      IF ((TRIM(state%metType) == "F515_516") .OR. &
          (TRIM(state%metType) == "F5131"))        state%metType = "FP"

! !IMPORT STATE:
      !------------------------------------------
      ! If we want vertical remapping must add PS
      !------------------------------------------
      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME        = 'PS',                           &
           LONG_NAME         = 'surface_pressure',             &
           UNITS             = 'Pa',                           &
           DIMS              = MAPL_DimsHorzOnly,              &
           VLOCATION         = MAPL_VLocationNone,    __RC__)

      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME        = 'AREA',                         &
           LONG_NAME         = 'agrid_cell_area',              &
           UNITS             = 'm+2',                          &
           DIMS              = MAPL_DimsHorzOnly,              &
           VLOCATION         = MAPL_VLocationNone,    __RC__)

      IF (state%read_advCoreFields) THEN
         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'MFX',                                       &
              LONG_NAME  = 'pressure_weighted_eastward_mass_flux',      &
              UNITS      = 'Pa m+2 s-1',                                &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'MFY',                                       &
              LONG_NAME  = 'pressure_weighted_northward_mass_flux',     &
              UNITS      = 'Pa m+2 s-1',                                &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'CX',                                        &
              LONG_NAME  = 'eastward_accumulated_courant_number',       &
              UNITS      = '',                                          &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'CY',                                        &
              LONG_NAME  = 'northward_accumulated_courant_number',      &
              UNITS      = '',                                          &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )
      ENDIF

      IF (state%read_PLE) THEN
         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'PLE0',                                      &
              LONG_NAME  = 'pressure_at_layer_edges_before_advection',  &
              UNITS      = 'Pa',                                        &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'PLE1',                                      &
              LONG_NAME  = 'pressure_at_layer_edges_after_advection',   &
              UNITS      = 'Pa',                                        &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,             __RC__  )
      ELSE
         call MAPL_AddImportSpec(GC,                                    &
              SHORT_NAME        = 'PS0',                                &
              LONG_NAME         = 'surface_pressure_before_advection',  &
              UNITS             = 'Pa',                                 &
              DIMS              = MAPL_DimsHorzOnly,                    &
              VLOCATION         = MAPL_VLocationNone,    __RC__)

         call MAPL_AddImportSpec(GC,                                    &
              SHORT_NAME        = 'PS1',                                &
              LONG_NAME         = 'surface_pressure_after_advection',   &
              UNITS             = 'Pa',                                 &
              DIMS              = MAPL_DimsHorzOnly,                    &
              VLOCATION         = MAPL_VLocationNone,    __RC__)
      ENDIF

      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'UC0',                                       &
           LONG_NAME  = 'eastward_wind_on_C-Grid_before_advection',  &
           UNITS      = 'm s-1',                                     &
           STAGGERING = MAPL_CGrid,                                  &
           ROTATION   = MAPL_RotateCube,                             & 
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'UC1',                                       &
           LONG_NAME  = 'eastward_wind_on_C-Grid_after_advection',   &
           UNITS      = 'm s-1',                                     &
           STAGGERING = MAPL_CGrid,                                  &
           ROTATION   = MAPL_RotateCube,                             &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'VC0',                                       &
           LONG_NAME  = 'northward_wind_on_C-Grid_before_advection', &
           UNITS      = 'm s-1',                                     &
           STAGGERING = MAPL_CGrid,                                  &
           ROTATION   = MAPL_RotateCube,                             &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'VC1',                                       &
           LONG_NAME  = 'northward_wind_on_C-Grid_after_advection',  &
           UNITS      = 'm s-1',                                     &
           STAGGERING = MAPL_CGrid,                                  &
           ROTATION   = MAPL_RotateCube,                             &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'T',                                         &
           LONG_NAME  = 'air_temperature',                           &
           UNITS      = 'K',                                         &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

     call MAPL_AddImportSpec(GC,                             &
          SHORT_NAME = 'TS',                                        &
          LONG_NAME  = 'surface temperature',                       &
          UNITS      = 'K',                                         &
          DIMS       = MAPL_DimsHorzOnly,                           &
          VLOCATION  = MAPL_VLocationNone,                 __RC__  )

     call MAPL_AddImportSpec(GC,                             &
          SHORT_NAME = 'FROCEAN',                                   &
          LONG_NAME  = 'areal_ocean_fraction',                      &
          UNITS      = '1',                                         &
          DIMS       = MAPL_DimsHorzOnly,                           &
          VLOCATION  = MAPL_VLocationNone,                  __RC__  )

      call MAPL_AddImportSpec(GC,                                    &
           SHORT_NAME = 'Q',                                         &
           LONG_NAME  = 'specific_humidity',                         &
           UNITS      = 'kg kg-1',                                   &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'FRLAKE',  &
           LONG_NAME          = 'fraction_of_lake',  &
           UNITS              = '1', &
           DIMS               = MAPL_DimsHorzOnly,    &
           VLOCATION          = MAPL_VLocationNone,     __RC__  )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'FRACI',  &
           LONG_NAME          = 'ice_covered_fraction_of_tile',  &
           UNITS              = '1', &
           DIMS               = MAPL_DimsHorzOnly,    &
           VLOCATION          = MAPL_VLocationNone,    __RC__  )

      IF ( (TRIM(state%metType) == 'MERRA2') .OR.  &
           (TRIM(state%metType) == 'FPIT')   .OR.  &
           (TRIM(state%metType) == 'FP')    ) THEN
         call MAPL_AddImportSpec(GC,                                 &
              SHORT_NAME         = 'ZLE',                            &
              LONG_NAME          = 'geopotential_height',            &
              UNITS              = 'm',                              &
              DIMS               = MAPL_DimsHorzVert,                &
              VLOCATION          = MAPL_VLocationEdge,       __RC__  )

         call MAPL_AddImportSpec(GC, &
              SHORT_NAME         = 'RH2',  &
              LONG_NAME          = 'relative_humidity_after_moist',  &
              UNITS              = '1', &
              DIMS               = MAPL_DimsHorzVert,    &
              VLOCATION          = MAPL_VLocationCenter,   __RC__  )
      END IF

      ! Only doing the Imports if we are not doing Idealized Passive Tracer
      !----------------------------------------------------------
      IF (.NOT. state%enable_pTracers) THEN
         IF (state%enable_rasCalculations) THEN
            call MAPL_AddImportSpec(GC, &
                 SHORT_NAME = 'ZPBL',  &
                 LONG_NAME  = 'planetary_boundary_layer_height', &
                 UNITS      = 'm',                               &
                 DIMS       = MAPL_DimsHorzOnly,                 &
                 VLOCATION  = MAPL_VLocationNone,         __RC__ )

            call MAPL_AddImportSpec(GC,                            &
                 SHORT_NAME = 'FRLAND',                            &
                 LONG_NAME  = 'fraction_of_land',                  &
                 UNITS      = '1',                                 &
                 DIMS       = MAPL_DimsHorzOnly,                   &
                 VLOCATION  = MAPL_VLocationNone,          __RC__  )

            call MAPL_AddImportSpec ( gc,                          &
                 SHORT_NAME = 'KH',                                &
                 LONG_NAME  = 'scalar_diffusivity',                &
                 UNITS      = 'm+2 s-1',                           &
                 DIMS       = MAPL_DimsHorzVert,                   &
                 VLOCATION  = MAPL_VLocationEdge,           __RC__ )
         ELSE
            call MAPL_AddImportSpec(GC, &
                 SHORT_NAME         = 'CNV_MFC',  &
                 LONG_NAME          = 'cumulative_mass_flux',  &
                 UNITS              = 'kg m-2 s-1', &
                 DIMS               = MAPL_DimsHorzVert,    &
                 VLOCATION          = MAPL_VLocationEdge,    __RC__ )

            call MAPL_AddImportSpec(GC, &
                 SHORT_NAME         = 'CNV_MFD',  &
                 LONG_NAME          = 'detraining_mass_flux',  &
                 UNITS              = 'kg m-2 s-1', &
                 DIMS               = MAPL_DimsHorzVert,    &
                 VLOCATION          = MAPL_VLocationCenter,    __RC__ )
         END IF

         call MAPL_AddImportSpec(GC,                             &
              SHORT_NAME  = 'QL',                                  &
              LONG_NAME   = 'mass_fraction_of_cloud_liquid_water', &
              UNITS       = 'kg kg-1',                             &
              DIMS        = MAPL_DimsHorzVert,                     &
              VLOCATION   = MAPL_VLocationCenter, __RC__)

         call MAPL_AddImportSpec(GC,                             &
              SHORT_NAME  = 'QI',                                  &
              LONG_NAME   = 'mass_fraction_of_cloud_ice_water',    &
              UNITS       = 'kg kg-1',                             &
              DIMS        = MAPL_DimsHorzVert,                     &
              VLOCATION   = MAPL_VLocationCenter, __RC__)

         call MAPL_AddImportSpec(GC,                             &
              SHORT_NAME  = 'PRECCON',                             &
              LONG_NAME   = 'convective_precipitation',            &
              UNITS       = 'km m-2 s-1',                          &
              DIMS        = MAPL_DimsHorzOnly,                     &
              VLOCATION   = MAPL_VLocationNone, __RC__)

         call MAPL_AddImportSpec(GC,                             &
              SHORT_NAME  = 'PRECANV',                             &
              LONG_NAME   = 'anvil_precipitation',                 &
              UNITS       = 'km m-2 s-1',                          &
              DIMS        = MAPL_DimsHorzOnly,                     &
              VLOCATION   = MAPL_VLocationNone, __RC__)

         call MAPL_AddImportSpec(GC,                             &
              SHORT_NAME  = 'PRECLSC',                             &
              LONG_NAME   = 'nonanvil_large_scale_precipitation',  &
              UNITS       = 'km m-2 s-1',                          &
              DIMS        = MAPL_DimsHorzOnly,                     &
              VLOCATION   = MAPL_VLocationNone, __RC__)

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'QLTOT',                                     &
              LONG_NAME  = 'mass_fraction_of_cloud_liquid_water',       &
              UNITS      = 'kg kg-1',                                   &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'QITOT',                                     &
              LONG_NAME  = 'mass_fraction_of_cloud_ice_water',          &
              UNITS      = 'kg kg-1',                                   &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'QLTOT1',                                    &
              LONG_NAME  = 'mass_fraction_of_cloud_liquid_water_before',&
              UNITS      = 'kg kg-1',                                   &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec ( gc,                                  &
              SHORT_NAME = 'QITOT1',                                    &
              LONG_NAME  = 'mass_fraction_of_cloud_ice_water_before',   &
              UNITS      = 'kg kg-1',                                   &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddImportSpec(GC,                                    &
              SHORT_NAME ='CN_PRCP',                                    &
              LONG_NAME ='convective_precipitation',                    &
              UNITS     ='kg m-2 s-1',                                  &
              DIMS      = MAPL_DimsHorzOnly,                            &
              VLOCATION = MAPL_VLocationNone,                __RC__  )
      END IF

! !EXPORT STATE:

      ! Default exports
      !----------------
      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'CXr8',                                      &
           LONG_NAME  = 'eastward_accumulated_courant_number',       &
           UNITS      = '',                                          &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'CYr8',                                      &
           LONG_NAME  = 'northward_accumulated_courant_number',      &
           UNITS      = '',                                          &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'MFXr8',                                     &
           LONG_NAME  = 'pressure_weighted_accumulated_eastward_mass_flux', &
           UNITS      = 'Pa m+2 s-1',                                &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'MFYr8',                                     &
           LONG_NAME  = 'pressure_weighted_accumulated_northward_mass_flux', &
           UNITS      = 'Pa m+2 s-1',                                &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'PLE1r8',                                    &
           LONG_NAME  = 'pressure_at_layer_edges_after_advection',   &
           UNITS      = 'Pa',                                        &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationEdge,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'PLE0r8',                                    &
           LONG_NAME  = 'pressure_at_layer_edges_before_advection',  &
           UNITS      = 'Pa',                                        &
           PRECISION  = ESMF_KIND_R8,                                &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationEdge,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'PLE',                                       &
           LONG_NAME  = 'pressure_at_layer_edges',                   &
           UNITS      = 'Pa',                                        &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationEdge,             __RC__  )

      call MAPL_AddExportSpec ( gc,                                  &
           SHORT_NAME = 'TH',                                        &
           LONG_NAME  = 'potential_temperature',                     &
           UNITS      = 'K',                                         &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationCenter,             __RC__  )

      call MAPL_AddExportSpec(GC,                               &
           SHORT_NAME         = 'AIRDENS',                      &
           LONG_NAME          = 'air_density',                  &
           UNITS              = 'kg m-3',                       &
           DIMS               = MAPL_DimsHorzVert,              &
           VLOCATION          = MAPL_VLocationCenter,  __RC__)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME         = 'LWI',  &
           LONG_NAME          = 'land-ocean-ice_mask',  &
           UNITS              = '1', &
           DIMS               = MAPL_DimsHorzOnly,    &
           VLOCATION          = MAPL_VLocationNone,  __RC__  )

      call MAPL_AddExportSpec(GC,                               &
           SHORT_NAME         = 'MASS',                         &
           LONG_NAME          = 'total_mass',                   &
           UNITS              = 'kg',                           &
           DIMS               = MAPL_DimsHorzVert,              &
           VLOCATION          = MAPL_VLocationCenter,  __RC__)

      call MAPL_AddExportSpec(GC,                                    &
           SHORT_NAME         = 'ZLE',                               &
           LONG_NAME          = 'geopotential_height',               &
           UNITS              = 'm',                                 &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationEdge,       __RC__  )

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME         = 'RH2',  &
           LONG_NAME          = 'relative_humidity_after_moist',  &
           UNITS              = '1', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationCenter,   __RC__  )

      ! Exports if not doing Passive Tracer experiment
      !-----------------------------------------------
      IF (.NOT. state%enable_pTracers) THEN
         call MAPL_AddExportSpec(GC,                                 &
              SHORT_NAME         = 'CNV_MFC',                        &
              LONG_NAME          = 'cumulative_mass_flux',           &
              UNITS              = 'kg m-2 s-1',                     &
              DIMS               = MAPL_DimsHorzVert,                &
              VLOCATION          = MAPL_VLocationEdge,        __RC__ )

         call MAPL_AddExportSpec(GC,                                 &
              SHORT_NAME         = 'CNV_MFD',                        &
              LONG_NAME          = 'detraining_mass_flux',           &
              UNITS              = 'kg m-2 s-1',                     &
              DIMS               = MAPL_DimsHorzVert,                &
              VLOCATION          = MAPL_VLocationCenter,      __RC__ )

         call MAPL_AddExportSpec ( gc,                            &
              SHORT_NAME = 'U',                                   &
              LONG_NAME  = 'eastward_wind',                       &
              UNITS      = 'm s-1',                               &
              DIMS       = MAPL_DimsHorzVert,                     &
              VLOCATION  = MAPL_VLocationCenter,           __RC__ )

         call MAPL_AddExportSpec ( gc,                            &
              SHORT_NAME = 'V',                                   &
              LONG_NAME  = 'northward_wind',                      &
              UNITS      = 'm s-1',                               &
              DIMS       = MAPL_DimsHorzVert,                     &
              VLOCATION  = MAPL_VLocationCenter,           __RC__ )
         call MAPL_AddExportSpec(GC,                              &
              SHORT_NAME= 'ITY',                                  &
              LONG_NAME = 'vegetation_type',                      &
              UNITS     = '1',                                    &
              DIMS      = MAPL_DimsHorzOnly,                      &
              VLOCATION = MAPL_VLocationNone,             __RC__  )

         call MAPL_AddExportSpec(GC,                              &
              SHORT_NAME='BYNCY',                                 &
              LONG_NAME ='buoyancy_of surface_parcel',            &
              UNITS     ='m s-2',                                 &
              DIMS      = MAPL_DimsHorzVert,                      &
              VLOCATION = MAPL_VLocationCenter,           __RC__  )

         call MAPL_AddExportSpec ( gc,                            &
              SHORT_NAME = 'CNV_QC',                              &
              LONG_NAME  = 'grid_mean_convective_condensate',     &
              UNITS      = 'kg kg-1',                             &
              DIMS       = MAPL_DimsHorzVert,                     &
              VLOCATION  = MAPL_VLocationCenter,          __RC__  )

         call MAPL_AddExportSpec ( gc,                            &
              SHORT_NAME = 'QCTOT',                               &
              LONG_NAME  = 'mass_fraction_of_total_cloud_water',  &
              UNITS      = 'kg kg-1',                             &
              DIMS       = MAPL_DimsHorzVert,                     &
              VLOCATION  = MAPL_VLocationCenter,          __RC__  )

         call MAPL_AddExportSpec(GC,                              &
              SHORT_NAME         = 'LFR',                         &
              LONG_NAME          = 'lightning_flash_rate',        &
              UNITS              = 'km-2 s-1',                    &
              DIMS               = MAPL_DimsHorzOnly,             &
              VLOCATION          = MAPL_VLocationNone,    __RC__  )

         call MAPL_AddExportSpec(GC,                           &
              SHORT_NAME = 'QLCN',                              &
              LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
              UNITS      = 'kg kg-1',                           &
              DIMS       = MAPL_DimsHorzVert,                   &
              VLOCATION  = MAPL_VLocationCenter,  __RC__)

         call MAPL_AddExportSpec(GC,                           &
              SHORT_NAME = 'QICN',                              &
              LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
              UNITS      = 'kg kg-1',                           &
              DIMS       = MAPL_DimsHorzVert,                   &
              VLOCATION  = MAPL_VLocationCenter, __RC__)

      END IF

      IF (.NOT. state%read_advCoreFields) THEN
         call MAPL_AddExportSpec ( gc,                                  &
              SHORT_NAME = 'MFX',                                       &
              LONG_NAME  = 'pressure_weighted_eastward_mass_flux',      &
              UNITS      = 'Pa m+2 s-1',                                &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddExportSpec ( gc,                                  &
              SHORT_NAME = 'MFY',                                       &
              LONG_NAME  = 'pressure_weighted_northward_mass_flux',     &
              UNITS      = 'Pa m+2 s-1',                                &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddExportSpec ( gc,                                  &
              SHORT_NAME = 'CX',                                        &
              LONG_NAME  = 'eastward_accumulated_courant_number',       &
              UNITS      = '',                                          &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )

         call MAPL_AddExportSpec ( gc,                                  &
              SHORT_NAME = 'CY',                                        &
              LONG_NAME  = 'northward_accumulated_courant_number',      &
              UNITS      = '',                                          &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationCenter,             __RC__  )
      END IF


      ! Set the Profiling timers
      !-------------------------
      call MAPL_TimerAdd(GC,    name="INITIALIZE"  , __RC__ )
      call MAPL_TimerAdd(GC,    name="RUN"         , __RC__ )

      ! Create children's gridded components and invoke their SetServices
      ! -----------------------------------------------------------------
      call MAPL_GenericSetServices    ( GC,  __RC__  )

      RETURN_(ESMF_SUCCESS)
  
      end subroutine SetServices
!
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize -- Initialized method for composite the CTMder
!
! !INTERFACE:
!
      subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT VARIABLES:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION: 
!  The Initialize method of the CTM Cinderella Component.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      __Iam__('Initialize')
      character(len=ESMF_MAXSTR)    :: COMP_NAME
      type(ESMF_Grid)               :: esmfGrid
      type (ESMF_VM)                :: VM
      type(MAPL_MetaComp), pointer  :: ggState      ! GEOS Generic State
      type (ESMF_Config)            :: CF
      type (T_CTMenv_STATE), pointer :: CTMenv_STATE
      type (CTMenv_WRAP)             :: WRAP
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !  Get my name and set-up traceback handle
      !  ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, VM=VM, __RC__ )
      Iam = TRIM(COMP_NAME)//"::Initialize"

      !  Initialize GEOS Generic
      !  ------------------------
      call MAPL_GenericInitialize ( gc, IMPORT, EXPORT, clock,  __RC__ )

      !  Get my internal MAPL_Generic state
      !  -----------------------------------
      call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

      call MAPL_TimerOn(ggSTATE,"TOTAL")
      call MAPL_TimerOn(ggSTATE,"INITIALIZE")

      ! Get my private state from the component
      !----------------------------------------
      call ESMF_UserCompGetInternalState(gc, 'CTMenv', WRAP, STATUS)
      VERIFY_(STATUS)

      CTMenv_STATE => WRAP%PTR

      IF (MAPL_AM_I_ROOT()) THEN
         PRINT*
         PRINT*, '<---------------------------------------------->'
         PRINT*, '<----            CTM Env Settings          ---->'
         PRINT*, '<---------------------------------------------->'
         PRINT*, TRIM(Iam)//' enable_pTracers:        ', CTMenv_STATE%enable_pTracers
         PRINT*, TRIM(Iam)//' do_ctmAdvection:        ', CTMenv_STATE%do_ctmAdvection
         PRINT*, TRIM(Iam)//' read_advCoreFields:     ', CTMenv_STATE%read_advCoreFields
         PRINT*, TRIM(Iam)//' enable_rasCalculations: ', CTMenv_STATE%enable_rasCalculations
         PRINT*, '<---------------------------------------------->'
         PRINT*
      END IF

      IF (CTMenv_STATE%enable_rasCalculations) THEN
         IF (MAPL_AM_I_ROOT()) &
                   PRINT*, TRIM(Iam)//': Doing RAS Calculations'
         CALL INIT_RASPARAMS(CTMenv_STATE%RASPARAMS)
      ENDIF

      call MAPL_TimerOff(ggSTATE,"INITIALIZE")
      call MAPL_TimerOff(ggSTATE,"TOTAL")

      RETURN_(ESMF_SUCCESS)

      end subroutine Initialize
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run -- Run method
!
! !INTERFACE:
!
      subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION: 
! The Run method of the CTM Cinderalla Component.
!
!EOP
!-------------------------------------------------------------------------
!BOC 
!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)      :: IAm = "Run"
      integer                         :: STATUS
      character(len=ESMF_MAXSTR)      :: COMP_NAME
      type (MAPL_MetaComp), pointer   :: ggState
      type (ESMF_Grid)                :: esmfGrid
      type (T_CTMenv_STATE), pointer  :: CTMenv_STATE
      type (CTMenv_wrap)              :: WRAP

      ! Imports
      !--------

      real, pointer, dimension(:,:,:) ::      PLE1 => null()
      real, pointer, dimension(:,:,:) ::      PLE0 => null()
      real, pointer, dimension(:,:,:) ::       UC0 => null()
      real, pointer, dimension(:,:,:) ::       UC1 => null()
      real, pointer, dimension(:,:,:) ::       VC0 => null()
      real, pointer, dimension(:,:,:) ::       VC1 => null()
      real, pointer, dimension(:,:,:) ::         T => null()
      real, pointer, dimension(:,:,:) ::         Q => null()
      real, pointer, dimension(:,:,:) ::       ZLE => null()
      real, pointer, dimension(:,:,:) ::     QITOT => null()
      real, pointer, dimension(:,:,:) ::     QLTOT => null()
      real, pointer, dimension(:,:,:) ::    QITOT1 => null()
      real, pointer, dimension(:,:,:) ::    QLTOT1 => null()
      real, pointer, dimension(:,:)   ::  cellArea => null()
      real, pointer, dimension(:,:)   ::       PS0 => null()
      real, pointer, dimension(:,:)   ::       PS1 => null()

      real, pointer, dimension(:,:)   ::        TS => null()
      real, pointer, dimension(:,:)   ::   FROCEAN => null()
      real, pointer, dimension(:,:)   ::   CN_PRCP => null()
      real, pointer, dimension(:,:,:) ::   CNV_MFC => null()
      real, pointer, dimension(:,:)   ::    FRLAKE => null()
      !real, pointer, dimension(:,:,:) ::     PFLCU => null()

      ! Exports
      !--------
      real,     pointer, dimension(:,:,:) ::        TH => null()
      real,     pointer, dimension(:,:,:) ::       PLE => null()
      real,     pointer, dimension(:,:,:) ::   AIRDENS => null()
      real,     pointer, dimension(:,:,:) ::      MASS => null()
      real(r8), pointer, dimension(:,:,:) ::      CXr8 => null()
      real(r8), pointer, dimension(:,:,:) ::      CYr8 => null()
      real(r8), pointer, dimension(:,:,:) ::    PLE1r8 => null()
      real(r8), pointer, dimension(:,:,:) ::    PLE0r8 => null()
      real(r8), pointer, dimension(:,:,:) ::     MFXr8 => null()
      real(r8), pointer, dimension(:,:,:) ::     MFYr8 => null()

      real(r4), pointer, dimension(:,:,:) ::      CXr4 => null()
      real(r4), pointer, dimension(:,:,:) ::      CYr4 => null()
      real(r4), pointer, dimension(:,:,:) ::     MFXr4 => null()
      real(r4), pointer, dimension(:,:,:) ::     MFYr4 => null()

      real(r8), pointer, dimension(:,:,:) ::      UCr8 => null()
      real(r8), pointer, dimension(:,:,:) ::      VCr8 => null()
      real(r8), pointer, dimension(:,:,:) ::     PLEr8 => null()

      !real,     pointer, dimension(:,:,:) ::         U => null()
      !real,     pointer, dimension(:,:,:) ::         V => null()
      real,     pointer, dimension(:,:,:) ::     QCTOT => null()
      real,     pointer, dimension(:,:,:) ::    CNV_QC => null()

      real,     pointer, dimension(:,:)   ::     FRACI => null()
      real,     pointer, dimension(:,:)   ::       LFR => null()

      real,     pointer, dimension(:,:,:) ::    RH2imp => null()
      real,     pointer, dimension(:,:,:) ::       RH2 => null()
      real,     pointer, dimension(:,:)   ::       LWI => null()
      real,     pointer, dimension(:,:)   ::       ITY => null()
      !real,     pointer, dimension(:,:,:) ::   totflux => null()
      real,     pointer, dimension(:,:,:) ::    expZLE => null()

      real,     pointer, dimension(:,:,:) ::        QI => null()
      real,     pointer, dimension(:,:,:) ::        QL => null()
      real,     pointer, dimension(:,:,:) ::      QICN => null()
      real,     pointer, dimension(:,:,:) ::      QLCN => null()
      real,     pointer, dimension(:,:)   ::   PRECCON => null()
      real,     pointer, dimension(:,:)   ::   PRECANV => null()
      real,     pointer, dimension(:,:)   ::   PRECLSC => null()
      !real,     pointer, dimension(:,:)   ::     TPREC => null()

      real,     pointer, dimension(:,:,:) ::      Uexp => null()
      real,     pointer, dimension(:,:,:) ::      Vexp => null()
      real,     pointer, dimension(:,:,:) ::  QITOTexp => null()
      real,     pointer, dimension(:,:,:) ::  QLTOTexp => null()

      real, pointer, dimension(:,:,:)     ::CNV_MFCexp => null()
      real, pointer, dimension(:,:,:)     ::CNV_MFDexp => null()
      real, pointer, dimension(:,:,:)     ::   CNV_MFD => null()
      real, pointer, dimension(:,:,:)     ::        KH => null()
      real, pointer, dimension(:,:)       ::    FRLAND => null()
      real, POINTER, dimension(:,:)       ::      PBLH => null()
      real, pointer, dimension(:)         ::      PREF => null()
      real, pointer, dimension(:,:,:)     ::CNV_MFCras => null()
      real, pointer, dimension(:,:,:)     ::CNV_MFDras => null()
      real(r8), pointer, dimension(:,:)   ::      LATS => null()
      real(r8), pointer, dimension(:,:)   ::      LONS => null()

      real(r8), allocatable             :: AK(:)
      real(r8), allocatable             :: BK(:)

      integer :: IM, JM, LM
      integer :: k, is, ie, js, je, nc
      integer :: ndt, isd, ied, jsd, jed, i, j, l
      real(r8) :: DT

      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------
      call ESMF_GridCompGet ( GC, name=COMP_NAME, Grid=esmfGrid, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)

      ! Get my internal MAPL_Generic state
      !-----------------------------------
      call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )

      call MAPL_TimerOn(ggState,"TOTAL")
      call MAPL_TimerOn(ggState,"RUN")

      ! Get my private state from the component
      !----------------------------------------
      call ESMF_UserCompGetInternalState(gc, 'CTMenv', WRAP, STATUS)
      VERIFY_(STATUS)

      CTMenv_STATE => WRAP%PTR

      ! Get the time-step
      ! -----------------------
      call MAPL_GetResource( ggState, ndt, 'RUN_DT:', default=0, __RC__ )
      DT = ndt

      !-----------------------------
      ! Required Imports and Exports
      !-----------------------------
      call MAPL_GetPointer ( IMPORT,      UC0,   'UC0', __RC__ )
      call MAPL_GetPointer ( IMPORT,      UC1,   'UC1', __RC__ )
      call MAPL_GetPointer ( IMPORT,      VC0,   'VC0', __RC__ )
      call MAPL_GetPointer ( IMPORT,      VC1,   'VC1', __RC__ )
      call MAPL_GetPointer ( IMPORT, cellArea,  'AREA', __RC__ )

      is = lbound(UC0,1); ie = ubound(UC0,1)
      js = lbound(UC0,2); je = ubound(UC0,2)
      IM = ie-is+1
      JM = je-js+1
      LM = size  (UC0,3)
      nc = (ie-is+1)*(je-js+1)

      call MAPL_GetPointer ( EXPORT,     PLE,    'PLE', __RC__ )
      call MAPL_GetPointer ( EXPORT,  PLE0r8, 'PLE0r8', __RC__ )
      call MAPL_GetPointer ( EXPORT,  PLE1r8, 'PLE1r8', __RC__ )

      IF (CTMenv_STATE%read_PLE) THEN
         call MAPL_GetPointer ( IMPORT,     PLE0,  'PLE0', __RC__ )
         call MAPL_GetPointer ( IMPORT,     PLE1,  'PLE1', __RC__ )
         PLE    = PLE0
         PLE0r8 = PLE0
         PLE1r8 = PLE1
      ELSE
         call MAPL_GetPointer ( IMPORT,      PS0,   'PS0', __RC__ )
         call MAPL_GetPointer ( IMPORT,      PS1,   'PS1', __RC__ )
         ! Compute the 3D pressure using surface pressure
         IF ( ASSOCIATED(PS0) .AND. ASSOCIATED(PS1) ) THEN

            allocate(AK(LM+1),stat=rc); VERIFY_(rc)
            allocate(BK(LM+1),stat=rc); VERIFY_(rc)

            call ESMF_AttributeGet(esmfGrid, name="GridAK", valuelist=AK, __RC__)
            call ESMF_AttributeGet(esmfGrid, name="GridBK", valuelist=BK, __RC__)

            call computeEdgePressure(PLE0r8, PS0, AK, BK, LM)
            call computeEdgePressure(PLE1r8, PS1, AK, BK, LM)

            DEALLOCATE(AK, BK)

            IF (ASSOCIATED(PLE)) PLE = PLE0r8
         ENDIF
      ENDIF

      !--------------------------------------------
      ! courant numbers and mass fluxes for AdvCore
      !--------------------------------------------

      IF (CTMenv_STATE%do_ctmAdvection) THEN
         call MAPL_GetPointer ( EXPORT,   MFXr8,  'MFXr8', __RC__ )
         call MAPL_GetPointer ( EXPORT,   MFYr8,  'MFYr8', __RC__ )
         call MAPL_GetPointer ( EXPORT,    CXr8,   'CXr8', __RC__ )
         call MAPL_GetPointer ( EXPORT,    CYr8,   'CYr8', __RC__ )

         IF (CTMenv_STATE%read_advCoreFields) THEN
            ! Get the courant numbers and mass fluxes from files
            call MAPL_GetPointer ( IMPORT,   MFXr4,  'MFX', __RC__ )
            call MAPL_GetPointer ( IMPORT,   MFYr4,  'MFY', __RC__ )
            call MAPL_GetPointer ( IMPORT,    CXr4,   'CX', __RC__ )
            call MAPL_GetPointer ( IMPORT,    CYr4,   'CY', __RC__ )

            MFXr8 = MFXr4
            MFYr8 = MFYr4
             CXr8 =  CXr4
             CYr8 =  CYr4
         ELSE
            ! Compute the courant numbers and mass fluxes
            ALLOCATE( UCr8(is:ie,js:je,lm),   STAT=STATUS); VERIFY_(STATUS)
            ALLOCATE( VCr8(is:ie,js:je,lm),   STAT=STATUS); VERIFY_(STATUS)
            ALLOCATE(PLEr8(is:ie,js:je,lm+1), STAT=STATUS); VERIFY_(STATUS)

            call A2D2C(UC1, VC1, LM, .TRUE.)
            call A2D2c(UC0, VC0, LM, .TRUE.)
            UCr8  = 0.50d0*(UC1 + UC0)
            VCr8  = 0.50d0*(VC1 + VC0)
            IF (CTMenv_STATE%read_PLE) THEN
               PLEr8 = 0.50d0*(PLE1 + PLE0)
            ELSE
               PLEr8 = 0.50d0*(PLE1r8 + PLE0r8)
            ENDIF

            call calcCourantNumberMassFlux(UCr8, VCr8, PLEr8, &
                                   MFXr8, MFYr8, CXr8, CYr8, DT)

            DEALLOCATE(UCr8, VCr8, PLEr8)
         ENDIF
      ENDIF

      !-----------
      ! Derive LWI
      !-----------
      call derive_LWI()

      !--------------------------------------
      ! Derive the Potential Temperature (TH)
      !--------------------------------------
      call derive_TH()

      !-----------------------------------
      ! Derive ZLE and RH2 if using MERRA1
      !-----------------------------------
      IF ( (TRIM(CTMenv_STATE%metType) == 'MERRA2') .OR.  &
           (TRIM(CTMenv_STATE%metType) == 'FPIT')   .OR.  &
           (TRIM(CTMenv_STATE%metType) == 'FP')    ) THEN
         call MAPL_GetPointer ( IMPORT,    ZLE,      'ZLE', __RC__ )
         call MAPL_GetPointer ( EXPORT, expZLE,      'ZLE', ALLOC=.TRUE., __RC__ )
         IF (ASSOCIATED(expZLE) .AND. ASSOCIATED(ZLE)) expZLE = ZLE

         call MAPL_GetPointer ( IMPORT, RH2imp,      'RH2', __RC__ )
         call MAPL_GetPointer ( EXPORT,    RH2,      'RH2', ALLOC=.TRUE., __RC__ )
         IF (ASSOCIATED(RH2imp) .AND. ASSOCIATED(RH2)) RH2 = RH2imp
      ELSEIF ( TRIM(CTMenv_STATE%metType) == 'MERRA1') THEN
         !---------------------------------
         ! RH2 and ZLE if using MERRA1 data
         !---------------------------------
         call MAPL_GetPointer ( EXPORT,    RH2,    'RH2', __RC__ )
         call MAPL_GetPointer ( EXPORT,    ZLE,    'ZLE', __RC__ )

         IF (ASSOCIATED(RH2) .AND. ASSOCIATED(ZLE)) THEN
            call compute_ZLE_RH2 (ZLE, RH2, TH, Q, PLE, ie-is+1, je-js+1, LM)
         END IF
      END IF

      ! ---------------------------------------
      ! Derive Air Density and Atmospheric Mass
      ! ---------------------------------------
      call derive_AIRDENS_MASS()

      IF (.NOT. CTMenv_STATE%enable_pTracers) THEN

         IF (ASSOCIATED(UC0)) THEN
            call MAPL_GetPointer ( EXPORT,     Uexp,      'U', ALLOC=.TRUE., __RC__ )
            Uexp = UC0
         ENDIF

         IF (ASSOCIATED(VC0)) THEN
            call MAPL_GetPointer ( EXPORT,     Vexp,      'V', ALLOC=.TRUE., __RC__ )
            Vexp = VC0
         ENDIF

         !---------------------
         ! Derive QICN and QLCN
         !---------------------
         call derive_QICN_QLCN()

         !-------------------------------------------
         ! Mass Fraction of Total Cloud Water (QCTOT)
         !-------------------------------------------
         call derive_QCTOT()

         !-------------------------------------------
         ! Grid Mean Convective Condensate (CNV_QC)
         !-------------------------------------------
         call derive_CNV_QC()

         !----------------
         ! Vegetation Type
         !----------------
         call MAPL_GetPointer ( EXPORT, ITY, 'ITY', ALLOC=.TRUE., __RC__ )
         ITY = 1.0

         !-----------------------
         ! Convective Mass Fluxes
         !-----------------------
         call derive_convective_mass_fluxes()

         !------------------------------------------------
         ! Flash Rate (LFR) for Lighting Parameterization
         ! Buoyancy
         !------------------------------------------------
         call derive_LFR_BYNCY()

      END IF ! .NOT. enable_pTracers

   
      IF (.NOT. CTMenv_STATE%read_advCoreFields) THEN
         call export_advCoreFields()
      END IF

      call MAPL_TimerOff(ggState,"RUN")
      call MAPL_TimerOff(ggState,"TOTAL")

      ! All Done
      ! --------
      RETURN_(ESMF_SUCCESS)

      !----------------------
      ! - Intenal Subroutines
      !----------------------
         CONTAINS
         !
         !-------------------------------------------------------
         !
         subroutine derive_LWI()
         call MAPL_GetPointer ( IMPORT, FROCEAN, 'FROCEAN', __RC__ )
         call MAPL_GetPointer ( IMPORT,   FRACI,   'FRACI', __RC__ )
         call MAPL_GetPointer ( IMPORT,  FRLAKE,  'FRLAKE', __RC__ )
         call MAPL_GetPointer ( IMPORT,       Q,       'Q', __RC__ )
         call MAPL_GetPointer ( IMPORT,      TS,      'TS', __RC__ )

         IF ( ASSOCIATED(FROCEAN) .AND. ASSOCIATED(FRACI) .AND. &
              ASSOCIATED(FRLAKE)  .AND. ASSOCIATED(TS)    .AND. &
              ASSOCIATED(Q) ) THEN
            call MAPL_GetPointer ( EXPORT,     LWI,     'LWI', ALLOC=.TRUE., __RC__ )
            call computeLWI      (LWI, TS, FRLAKE, FROCEAN, FRACI)
         ENDIF
         end subroutine derive_LWI
         !
         !-------------------------------------------------------
         !
         subroutine derive_TH()
         real,     pointer, dimension(:,:,:) ::    PKAPPA => null()
         call MAPL_GetPointer ( IMPORT,      T,     'T',  __RC__  )
         IF ( ASSOCIATED(T) .AND. ASSOCIATED(PLE) ) THEN
            call MAPL_GetPointer ( EXPORT,     TH,    'TH', ALLOC=.TRUE.,  __RC__  )

            ALLOCATE( PKAPPA(is:ie,js:je,LM),   STAT=STATUS); VERIFY_(STATUS)

            PKAPPA(:,:,:) = ((0.5*(PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) ))/100000.)**(MAPL_RGAS/MAPL_CP)
            TH(:,:,:) = T(:,:,:)/PKAPPA(:,:,:)

            DEALLOCATE(PKAPPA)
         ENDIF
         end subroutine derive_TH
         !
         !-------------------------------------------------------
         !
         subroutine derive_AIRDENS_MASS()

         IF ( ASSOCIATED(PLE) .AND. ASSOCIATED(TH) .AND. ASSOCIATED(Q) .AND. &
              ASSOCIATED(ZLE) ) THEN
            call MAPL_GetPointer ( EXPORT, AIRDENS, 'AIRDENS', ALLOC=.TRUE., __RC__ )
            call MAPL_GetPointer ( EXPORT,    MASS,    'MASS', ALLOC=.TRUE., __RC__ )

            ! Compute air density
            call airdens_ ( AIRDENS, PLE, TH, Q, ie-is+1, je-js+1, LM)

            ! Compute the total mass
            DO k = 1, LM
               MASS(:,:,k) = AIRDENS(:,:,k)*cellArea(:,:)*(ZLE(:,:,k-1)-ZLE(:,:,k))
            END DO
         ENDIF
         end subroutine derive_AIRDENS_MASS
         !
         !-------------------------------------------------------
         !
         subroutine derive_LFR_BYNCY()
            ! Imports: CNV_MFC, CN_PRCP
            ! Exports: LFR, BYNCY
         real,     pointer, dimension(:,:,:) ::     BYNCY => null()
         real,     pointer, dimension(:,:)   ::  CNV_TOPP => null()
         real,     pointer, dimension(:,:)   ::      CAPE => null()
         real,     pointer, dimension(:,:)   :: flashRate => null()

         call MAPL_GetPointer ( IMPORT,  CNV_MFC, 'CNV_MFC', __RC__ )
         call MAPL_GetPointer ( IMPORT,  CN_PRCP, 'CN_PRCP', __RC__ )

         IF (ASSOCIATED(CNV_MFC) .AND. ASSOCIATED(CN_PRCP)) THEN
            call MAPL_GetPointer ( EXPORT,      LFR,     'LFR', ALLOC=.TRUE., __RC__ )
            call MAPL_GetPointer ( EXPORT,    BYNCY,   'BYNCY', ALLOC=.TRUE., __RC__ )

            ! Determine the pressure at convective cloud top
            ALLOCATE( CNV_TOPP(is:ie,js:je),   STAT=STATUS); VERIFY_(STATUS)
            CNV_TOPP(:,:) = MAPL_UNDEF
            do j=js, je
               do i=is, ie
                  do l=1,lm
                     if (CNV_MFC(i,j,l)/=0.0) then
                        CNV_TOPP(i,j) = PLE(i,j,l)
                        exit
                     endif
                  enddo
               enddo
            enddo

            ALLOCATE( CAPE(is:ie,js:je),   STAT=STATUS); VERIFY_(STATUS)
            call computeCAPE (TH, Q, PLE, CAPE, BYNCY, ie-is+1, je-js+1, LM)

            ALLOCATE( flashRate(is:ie,js:je),   STAT=STATUS); VERIFY_(STATUS)

            call computeFlashRate (ggState, nc, LM, TS, CNV_TOPP, FROCEAN, &
                          CN_PRCP, CAPE, CNV_MFC, TH, PLE, ZLE, flashRate, RC=STATUS)
            VERIFY_(STATUS)

            LFR(:,:) = flashRate(:,:)

            DEALLOCATE(CNV_TOPP, CAPE, flashRate)
         ENDIF
         end subroutine derive_LFR_BYNCY
         !
         !-------------------------------------------------------
         !
         subroutine derive_convective_mass_fluxes()
            ! Imports:
            ! Exports: CNV_MFC, CNV_MFD
         IF (CTMenv_STATE%enable_rasCalculations) THEN

            ! Get the LATS and LONS
            call ESMF_GridGetCoord(esmfGrid, coordDim=2, localDE=0, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   farrayPtr=LATS, rc=status)
            VERIFY_(status)

            call ESMF_GridGetCoord(esmfGrid, coordDim=1, localDE=0, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   farrayPtr=LONS, rc=status)

            !-----------------------------------------
            ! Compute the reference air pressure in Pa
            ! AK in Pa
            ! MAPL_P00 = 100000 Pa
            !-----------------------------------------
            ALLOCATE(PREF(LM+1))
            allocate(AK(LM+1),stat=rc); VERIFY_(rc)
            allocate(BK(LM+1),stat=rc); VERIFY_(rc)

            call ESMF_AttributeGet(esmfGrid, name="GridAK", valuelist=AK, __RC__)
            call ESMF_AttributeGet(esmfGrid, name="GridBK", valuelist=BK, __RC__)

            PREF = AK + BK * MAPL_P00

            ALLOCATE(CNV_MFDras(IM,JM,1:LM))
            ALLOCATE(CNV_MFCras(IM,JM,0:LM))

            call MAPL_GetPointer ( IMPORT,      KH,        'KH',  __RC__ )
            call MAPL_GetPointer ( IMPORT,    PBLH,      'ZPBL',  __RC__ )
            call MAPL_GetPointer ( IMPORT,  FRLAND,    'FRLAND',  __RC__ )

            IF ( ASSOCIATED(KH) .AND. ASSOCIATED(PBLH) .AND. &
                 ASSOCIATED(FRLAND) ) THEN
               CALL DO_RAS(CTMenv_STATE%RASPARAMS, PREF, TH, PLE, KH, PBLH, Q, &
                           FRLAND, TS, CNV_MFDras, CNV_MFCras,      &
                           REAL(LATS), REAL(LONS), REAL(DT), IM, JM, LM)

               call MAPL_GetPointer ( EXPORT, CNV_MFC, 'CNV_MFC',  __RC__ )
               call MAPL_GetPointer ( EXPORT, CNV_MFD, 'CNV_MFD',  __RC__ )

               CNV_MFC(:,:,:) = CNV_MFCras(:,:,:)
               CNV_MFD(:,:,:) = CNV_MFDras(:,:,:)
            ENDIF
            DEALLOCATE(PREF)
            DEALLOCATE(CNV_MFDras, CNV_MFCras)
            DEALLOCATE(AK, BK)
         ELSE
            call MAPL_GetPointer ( IMPORT,  CNV_MFC,  'CNV_MFC',  __RC__ )
            IF (ASSOCIATED(CNV_MFC)) THEN
               call MAPL_GetPointer ( EXPORT,  CNV_MFCexp,  'CNV_MFC',  __RC__ )
               CNV_MFCexp(:,:,:) = CNV_MFC(:,:,:)
            ENDIF

            call MAPL_GetPointer ( IMPORT,  CNV_MFD,  'CNV_MFD',  __RC__ )
            IF (ASSOCIATED(CNV_MFD)) THEN
               call MAPL_GetPointer ( EXPORT,  CNV_MFDexp,  'CNV_MFD',  __RC__ )
               CNV_MFDexp(:,:,:) = CNV_MFD(:,:,:)
            ENDIF
         ENDIF
         end subroutine derive_convective_mass_fluxes
         !
         !-------------------------------------------------------
         !
         subroutine derive_QCTOT()
            ! Imports: QITOT, QLTOT
            ! Exports: QCTOT
         call MAPL_GetPointer ( IMPORT,    QITOT,  'QITOT', __RC__  )
         call MAPL_GetPointer ( IMPORT,    QLTOT,  'QLTOT',  __RC__  )

         IF (ASSOCIATED(QITOT) .AND. ASSOCIATED(QLTOT) ) THEN
            call MAPL_GetPointer ( EXPORT,    QCTOT,  'QCTOT', ALLOC=.TRUE.,  __RC__  )
            QCTOT(:,:,:) = QLTOT(:,:,:) + QITOT(:,:,:)
         ENDIF
         end subroutine derive_QCTOT
         !
         !-------------------------------------------------------
         !
         subroutine derive_CNV_QC()
            ! Expected: QCTOT derived in derive_QCTOT
            ! Imports: QITOT1, QLTOT1 
            ! Exports: CNV_QC
         call MAPL_GetPointer ( IMPORT,   QITOT1, 'QITOT1', __RC__  )
         call MAPL_GetPointer ( IMPORT,   QLTOT1, 'QLTOT1', __RC__  )
         IF (ASSOCIATED(QITOT) .AND. ASSOCIATED(QLTOT1) .AND. ASSOCIATED(QITOT1) ) THEN
            call MAPL_GetPointer ( EXPORT,   CNV_QC, 'CNV_QC',  ALLOC=.TRUE., __RC__  )
            CNV_QC(:,:,:) = QCTOT(:,:,:) - ( QLTOT1(:,:,:) + QITOT1(:,:,:) )

            WHERE ( CNV_QC(:,:,:) < 0.0 ) CNV_QC(:,:,:) = 0.0
            CNV_QC(:,:,:) = CNV_QC(:,:,:) *(DT/(30.0*SecondsPerMinute))
         ENDIF
         end subroutine derive_CNV_QC
         !
         !-------------------------------------------------------
         !
         subroutine  derive_QICN_QLCN()

         real,     pointer, dimension(:,:)   ::     TPREC => null()
         real,     pointer, dimension(:,:)   ::  tempFrac => null()

         call MAPL_GetPointer ( IMPORT,  PRECCON, 'PRECCON', __RC__ )
         call MAPL_GetPointer ( IMPORT,  PRECANV, 'PRECANV', __RC__ )
         call MAPL_GetPointer ( IMPORT,  PRECLSC, 'PRECLSC', __RC__ )
         call MAPL_GetPointer ( IMPORT,       QI,      'QI', __RC__ )
         call MAPL_GetPointer ( IMPORT,       QL,      'QL', __RC__ )

         IF (ASSOCIATED(QI) .AND. ASSOCIATED(QL)) THEN
            call MAPL_GetPointer ( EXPORT,     QICN,    'QICN',  __RC__ )
            call MAPL_GetPointer ( EXPORT,     QLCN,    'QLCN',  __RC__ )

            isd = lbound(QLCN,1); ied = ubound(QLCN,1)
            jsd = lbound(QLCN,2); jed = ubound(QLCN,2)

            ALLOCATE( tempFrac(isd:ied,jsd:jed), STAT=STATUS); VERIFY_(STATUS)
            ALLOCATE(    TPREC(isd:ied,jsd:jed), STAT=STATUS); VERIFY_(STATUS)

            call computeTotalPrecip(TPREC, PRECANV, PRECCON, PRECLSC)

            tempFrac(:,:) = 0.0
            WHERE (TPREC(:,:) .NE. 0.0)
                tempFrac(:,:) = PRECANV(:,:) /TPREC(:,:)
            END WHERE

            DO k=1,LM
               QICN(:,:,k) = QI(:,:,k)*tempFrac(:,:)
               QLCN(:,:,k) = QL(:,:,k)*tempFrac(:,:)
            ENDDO
            DEALLOCATE(tempFrac, TPREC)
         ENDIF

         end subroutine  derive_QICN_QLCN
         !
         !-------------------------------------------------------
         !
         subroutine export_advCoreFields()
         call MAPL_GetPointer ( EXPORT,   MFXr4,  'MFX', __RC__ )
         IF (ASSOCIATED(MFXr8).AND.ASSOCIATED(MFXr4)) MFXr4 = MFXr8

         call MAPL_GetPointer ( EXPORT,   MFYr4,  'MFY', __RC__ )
         IF (ASSOCIATED(MFYr8).AND.ASSOCIATED(MFYr4)) MFYr4 = MFYr8

         call MAPL_GetPointer ( EXPORT,   CXr4,  'CX', __RC__ )
         IF (ASSOCIATED(CXr8).AND.ASSOCIATED(CXr4)) CXr4 = CXr8

         call MAPL_GetPointer ( EXPORT,   CYr4,  'CY', __RC__ )
         IF (ASSOCIATED(CYr8).AND.ASSOCIATED(CYr4)) CYr4 = CYr8
         end subroutine export_advCoreFields

      end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine computeEdgePressure(PLE, PS, AK, BK, km)
!
! !INPUT PARAMETERS:
      INTEGER,  intent(in) :: km      ! number of vertical levels
      REAL(r4), intent(in) :: PS(:,:) ! Surface pressure (Pa)
      REAL(r8), intent(in) :: ak(km+1), bk(km+1)
!
! !OUTPUT PARAMETERS:
      REAL(r8), intent(out) :: PLE(:,:,:)  ! Edge pressure (Pa)
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER  :: L

      DO L = 1, km+1
         PLE(:,:,L) = ak(L) + bk(L)*PS(:,:)
      END DO

      RETURN

      end subroutine computeEdgePressure
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine computeTotalPrecip(TPREC, PRECANV, PRECCON, PRECLSC)
!
! !INPUT PARAMETERS:
     REAL(r4), intent(in) :: PRECANV(:,:) ! Surface precipitation flux from anvils (kg/m2/s)
     REAL(r4), intent(in) :: PRECCON(:,:) ! Surface precipitation flux from convection (kg/m2/s)
     REAL(r4), intent(in) :: PRECLSC(:,:) ! Surface precipitation flux from large-scale (kg/m2/s)
!
! !OUTPUT PARAMETERS:
     REAL(r4), intent(out) :: TPREC(:,:)  ! Total precipitation (kg/m2/s)
!EOP
!------------------------------------------------------------------------------
!BOC
      TPREC = PRECANV + PRECCON + PRECLSC

      RETURN

      end subroutine computeTotalPrecip
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine computeLWI(LWI, TSKIN, FRLAKE, FROCEAN, FRACI)
!
! !INPUT PARAMETERS:
     REAL(r4), intent(in) :: TSKIN(:,:)    ! Surface skin temperature (K)
     REAL(r4), intent(in) :: FRLAKE(:,:)   ! Fraction of lake type in grid box (1)
     REAL(r4), intent(in) :: FROCEAN(:,:)  ! Fraction of ocean in grid box (1)
     REAL(r4), intent(in) :: FRACI(:,:)    ! Ice covered fraction of tile (1)
!
! !OUTPUT PARAMETERS:
     REAL(r4), intent(out) :: LWI(:,:) ! Land water ice flag (1)
!
!EOP
!------------------------------------------------------------------------------
!BOC

                                          LWI = 1.0  ! Land
      where ( FROCEAN+FRLAKE >= 0.6     ) LWI = 0.0  ! Water
      where ( LWI==0 .and. FRACI>0.5    ) LWI = 2.0  ! Ice
      where ( LWI==0 .and. TSKIN<271.40 ) LWI = 2.0  ! Ice

      RETURN

      end subroutine computeLWI
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine computeRelativeHumidity(RH2, PRESS3D, T, QV)

!
! !INPUT PARAMETERS:
      REAL, intent(in) :: PRESS3D(:,:,:)  ! Pressure (Pa)
      REAL, intent(in) :: T      (:,:,:)  ! Air temperature (K)
      REAL, intent(in) :: QV     (:,:,:)  ! Specific humidity (kg/kg)
!
! !OUTPUT PARAMETERS:
      REAL, intent(out) :: RH2(:,:,:) ! Relative humidity (1)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! -----------------------------------------------------------------
      ! First calculate relative humidity from Seinfeld (1986) p. 181.
      ! The first  RH2 is the temperature dependent parameter a.
      ! The second RH2 is the saturation vapor pressure of water.
      ! The third  RH2 is the actual relative humidity as a fraction.
      ! Then make sure RH2 is between 0 and 0.95.
      !-----------------------------------------------------------------

      RH2(:,:,:) = 1.0d0 - (373.15d0 / T(:,:,:))

      RH2(:,:,:) =  &
             1013.25d0 * Exp (13.3185d0 * RH2(:,:,:)    -  &
                               1.9760d0 * RH2(:,:,:)**2 -  &
                               0.6445d0 * RH2(:,:,:)**3 -  &
                               0.1299d0 * RH2(:,:,:)**4)

      RH2(:,:,:) = QV(:,:,:) * MWTAIR / 18.0d0 /  &
                      GPKG * PRESS3D(:,:,:) / RH2(:,:,:)

      RH2(:,:,:) = Max (Min (RH2(:,:,:), 0.95d0), 0.0d0)

      RETURN 

      end subroutine computeRelativeHumidity
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINES: airdens
!
! !INTERFACE:

      subroutine airdens_ ( AIRDENS, PLE, TH, Q, im, jm, lm )
!
! !INPUT PARAMETERS:
      integer, intent(in) :: im, jm , lm
      real,    intent(in) :: PLE(im,jm,lm+1)   ! pressure edges
      real,    intent(in) :: TH(im,jm,lm)      ! (dry) potential temperature
      real,    intent(in) :: Q(im,jm,lm)       ! apecific humidity
!
! !OUTPUT PARAMETERS:
      real,    intent(out) :: AIRDENS(:,:,:)    ! air density [kg/m3]
!
! !DESCRIPTION:
! Computes the air density that might be needed when GEOSchem is not
! exercised.
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: k
      real :: eps
      integer :: STATUS, RC
      character(len=ESMF_MAXSTR)      :: IAm = "airdens_"
      real, allocatable :: npk(:,:,:) ! normalized pk = (PLE/p0)^kappa

      allocate(npk(im,jm,lm+1),stat=STATUS) ! work space
      VERIFY_(STATUS)

      eps = MAPL_RVAP / MAPL_RGAS - 1.0

      ! Compute normalized PLE**Kappa
      ! ----------------------------
      npk = (PLE/MAPL_P00)**MAPL_KAPPA

      ! Compute AIRDENS from hydrostatic equation
      ! -------------------------------------
      do k = 1, lm
         AIRDENS(:,:,k) =       ( PLE(:,:,k+1) - PLE(:,:,k) ) /      &
                      ( MAPL_CP * ( TH(:,:,k)*(1. + eps*Q(:,:,k) ) ) &
                              * ( npk(:,:,k+1) - npk(:,:,k) ) )
      end do

      deallocate(npk)

      end subroutine airdens_
!EOC
!-----------------------------------------------------------------------
!BOP
      function compAreaWeightedAverage_2d (var2D, vm, cellArea) result(wAverage)
!
! !INPUT PARAMETER:
      real            :: var2D(:,:)
      real            :: cellArea(:,:)
      type (ESMF_VM)  :: VM
!
! RETURNED VALUE:
      real  :: wAverage
!
! DESCRIPTION:
! Computes the area weighted average of a 2d variable.
!
!EOP
!-----------------------------------------------------------------------
!BOC
      logical, save :: first = .true.
      real(r8) , save :: sumArea
      real(r8) :: sumWeight
      integer :: im, jm, STATUS, RC
      real(r8), pointer :: weightVals(:,:)
      real(r8) :: sumWeight_loc, sumArea_loc
      character(len=ESMF_MAXSTR) :: IAm = 'compAreaWeightedAverage_2d'

      ! Determine the earth surface area
      if (first) then
         sumArea_loc   = SUM( cellArea  (:,:)  )
         call MAPL_CommsAllReduceSum(vm, sendbuf= sumArea_loc, &
                                         recvbuf= sumArea, &
                                         cnt=1, RC=status)
         VERIFY_(STATUS)

         first = .false.
      end if

      im = size(cellArea,1)
      jm = size(cellArea,2)

      allocate(weightVals(im,jm))
      weightVals(:,:) = cellArea(:,:)*var2D(:,:)

      sumWeight_loc = SUM( weightVals(:,:) )

      call MAPL_CommsAllReduceSum(vm, sendbuf= sumWeight_loc, recvbuf= sumWeight, &
         cnt=1, RC=status)
      VERIFY_(STATUS)

      wAverage = sumWeight/sumArea

      deallocate(weightVals)

      return

      end function compAreaWeightedAverage_2d
!EOC
!-----------------------------------------------------------------------
!BOP
      function compAreaWeightedAverage_3d (var3D, vm, cellArea) result(wAverage)
!
! !INPUT PARAMETER:
      real            :: var3D(:,:,:)
      real            :: cellArea(:,:)
      type (ESMF_VM)  :: VM
!
! RETURNED VALUE:
      real  :: wAverage
!
! DESCRIPTION:
! Computes the area weighted average of a 3d variable.
!
!EOP
!-----------------------------------------------------------------------
!BOC
      logical, save :: first = .true.
      real(r8) , save :: sumArea
      real(r8) :: sumWeight
      integer :: ik, im, jm, STATUS, RC
      real(r8), pointer :: weightVals(:,:)
      real(r8) :: sumWeight_loc, sumArea_loc
      character(len=ESMF_MAXSTR) :: IAm = 'compAreaWeightedAverage_3d'

      ! Determine the earth surface area
      if (first) then
         sumArea_loc   = SUM( cellArea  (:,:)  )
         call MAPL_CommsAllReduceSum(vm, sendbuf= sumArea_loc, &
                                         recvbuf= sumArea, &
                                         cnt=1, RC=status)
         VERIFY_(STATUS)

         first = .false.
      end if

      im = size(cellArea,1)
      jm = size(cellArea,2)

      allocate(weightVals(im,jm))
      weightVals(:,:) = 0.0d0
      DO ik = lbound(var3D,3), ubound(var3D,3)
         weightVals(:,:) = weightVals(:,:) + cellArea(:,:)*var3D(:,:,ik)
      END DO

      sumWeight_loc = SUM( weightVals(:,:) )

      call MAPL_CommsAllReduceSum(vm, sendbuf= sumWeight_loc, recvbuf= sumWeight, &
         cnt=1, RC=status)
      VERIFY_(STATUS)

      wAverage = sumWeight/sumArea

      deallocate(weightVals)

      return

      end function compAreaWeightedAverage_3d
!EOC
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: computeFlashRate
!
! !INTERFACE:
!
      subroutine computeFlashRate (STATE, nc, lm, TS, CCTP, FROCEAN, CN_PRCP, &
                        CAPE, CNV_MFC, TH, PLE, ZLE, strokeRate, RC)
!
! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: nc     ! Number of cells
      INTEGER, INTENT(IN) :: lm     ! Number of layers
    
      REAL, INTENT(IN), DIMENSION(nc) :: TS       ! Surface temperature [K]
      REAL, INTENT(IN), DIMENSION(nc) :: CCTP     ! Convective cloud top pressure [Pa] with MAPL_UNDEFs
      REAL, INTENT(IN), DIMENSION(nc) :: FROCEAN  ! Areal ocean fraction
      REAL, INTENT(IN), DIMENSION(nc) :: CN_PRCP  ! Convective precipitation [kg m^{-2} s^{-1}]
      REAL, INTENT(IN), DIMENSION(nc) :: CAPE     ! Convective available potential energy [J m^{-2}]

      REAL, INTENT(IN), DIMENSION(nc,lm) :: TH        ! Potential temperature [K]
      REAL, INTENT(IN), DIMENSION(nc,0:lm) :: CNV_MFC ! Convective mass flux [kg m^{-2} s^{-1}]
      REAL, INTENT(IN), DIMENSION(nc,0:lm) :: PLE     ! Layer interface pressures  [Pa]
      REAL, INTENT(IN), DIMENSION(nc,0:lm) :: ZLE     ! Layer depths [m]

!
! !OUTPUT PARAMETERS:
      REAL, INTENT(OUT), DIMENSION(nc) :: strokeRate  ! Flashes per second
      INTEGER, OPTIONAL, INTENT(OUT) :: RC
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(MAPL_MetaComp), POINTER :: STATE ! Internal MAPL_Generic state
!
! !DESCRIPTION:
!  Generate lightning flash rates [km$^{-2}$ s$^{-1}$] using a six-variable polynomial fit.\\
!
!
!  ORIGIN AND CONTACT\\
!  Dr. Dale Allen, Associate Research Scientist\\
!  Dept. of Atmospheric and Oceanic Science\\
!  University of Maryland\\
!  College Park, MD 20742\\
!  301-405-7629 (ph); 301-314-9482 (fax)\\
!  http://www.meto.umd.edu/~allen\\
!  
!
!  FORMULATION NOTES\\
!  Predictor variables are set to zero where CN\_PRCP is zero or where the 
!   optical depth cloud top height is less than 5.5 km.
!  The fit returns flash rates in units km$^{-2}$ day$^{-1}$.  Convert to 
!   km$^{-2}$ s$^{-1}$ for the export state.\\
!
!
!  OTHER NOTES OF INTEREST\\
!  MOIST sets CNV\_TOPP to zero if there is an absence of convection.
!
! !REVISION HISTORY:
! 30 Nov 2011 Nielsen     First crack
! 29 Feb 2012 Nielsen     Accomodate CNV\_TOPP MAPL\_UNDEF for and after Fortuna-2\_5\_p4
! 04 Nov 2014 Kouatchou   Adapted the subroutine for GEOSctm
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      ! Error log variables
      ! -------------------
      INTEGER :: STATUS

      REAL, DIMENSION(nc) :: A1X1
      REAL, DIMENSION(nc) :: A2X2
      REAL, DIMENSION(nc) :: A3X3
      REAL, DIMENSION(nc) :: A4X4
      REAL, DIMENSION(nc) :: A5X5

      ! Local variables
      ! ---------------
      INTEGER :: i            ! General-purpose integers
      INTEGER :: k
      INTEGER :: n
    
      REAL :: a0c,a0m         ! Coefficients at continental and marine locations
      REAL :: a1c,a1m
      REAL :: a2c,a2m
      REAL :: a3c,a3m
      REAL :: a4c,a4m
      REAL :: a5c,a5m
    
      REAL :: x1Divisor       ! Divisors for x1-x5.
      REAL :: x2Divisor
      REAL :: x3Divisor
      REAL :: x4Divisor
      REAL :: x5Divisor
    
      REAL :: x5Power         ! Exponent for the surface temperature deviation predictor
    
      REAL :: sfcTLimit       ! Temperature thresholds
      REAL :: airTLimit
    
      REAL :: hPaCldTop       ! Cloud top limiter for weak/no convection
    
      REAL, ALLOCATABLE, DIMENSION(:) :: x1         ! Five independent variables
      REAL, ALLOCATABLE, DIMENSION(:) :: x2
      REAL, ALLOCATABLE, DIMENSION(:) :: x3
      REAL, ALLOCATABLE, DIMENSION(:) :: x4
      REAL, ALLOCATABLE, DIMENSION(:) :: x5
    
      REAL, ALLOCATABLE, DIMENSION(:) :: cloudTopAG ! Cloud top height above ground
      REAL, ALLOCATABLE, DIMENSION(:) :: cnv_topp   ! Convective cloud top pressure with MAPL_UNDEFs
                                                    ! changed to zero
    
      REAL, ALLOCATABLE, DIMENSION(:,:) :: dZ       ! Layer depths [m]
      REAL, ALLOCATABLE, DIMENSION(:,:) :: p        ! Pressure at middle of layer [Pa]
      REAL, ALLOCATABLE, DIMENSION(:,:) :: T        ! Air temperature at middle of layer [K]
    
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: weakCnvMask   ! Weak or no convection mask
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask          ! Working mask
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cloudTopMask  ! Mask is 1 below cloud top

      character(len=ESMF_MAXSTR)    :: IAm = 'computeFlashRate'


      ! Preliminaries
      ! -------------
      RC = 0
      strokeRate(:) = 0.0

      ! Coefficients of the predictors, marine locations
      ! ------------------------------------------------
      CALL MAPL_GetResource(STATE,a0m,'MARINE_A0:',DEFAULT= 0.0139868, __RC__)
      CALL MAPL_GetResource(STATE,a1m,'MARINE_A1:',DEFAULT= 0.0358764, __RC__)
      CALL MAPL_GetResource(STATE,a2m,'MARINE_A2:',DEFAULT=-0.0610214, __RC__)
      CALL MAPL_GetResource(STATE,a3m,'MARINE_A3:',DEFAULT=-0.0102320, __RC__)
      CALL MAPL_GetResource(STATE,a4m,'MARINE_A4:',DEFAULT= 0.0031352, __RC__)
      CALL MAPL_GetResource(STATE,a5m,'MARINE_A5:',DEFAULT= 0.0346241, __RC__)
    
      ! Coefficients of the predictors, continental locations
      ! -----------------------------------------------------
      CALL MAPL_GetResource(STATE,a0c,'CONTINENT_A0:',DEFAULT=-0.0183172, __RC__)
      CALL MAPL_GetResource(STATE,a1c,'CONTINENT_A1:',DEFAULT=-0.0562338, __RC__)
      CALL MAPL_GetResource(STATE,a2c,'CONTINENT_A2:',DEFAULT= 0.1862740, __RC__)
      CALL MAPL_GetResource(STATE,a3c,'CONTINENT_A3:',DEFAULT=-0.0023363, __RC__)
      CALL MAPL_GetResource(STATE,a4c,'CONTINENT_A4:',DEFAULT=-0.0013838, __RC__)
      CALL MAPL_GetResource(STATE,a5c,'CONTINENT_A5:',DEFAULT= 0.0114759, __RC__)
    
      ! Divisors for nondimensionalization of the predictors
      ! ----------------------------------------------------
      CALL MAPL_GetResource(STATE,x1Divisor,'X1_DIVISOR:',DEFAULT=4.36, __RC__)
      CALL MAPL_GetResource(STATE,x2Divisor,'X2_DIVISOR:',DEFAULT=9.27, __RC__)
      CALL MAPL_GetResource(STATE,x3Divisor,'X3_DIVISOR:',DEFAULT=34.4, __RC__)
      CALL MAPL_GetResource(STATE,x4Divisor,'X4_DIVISOR:',DEFAULT=21.4, __RC__)
      CALL MAPL_GetResource(STATE,x5Divisor,'X5_DIVISOR:',DEFAULT=14600., __RC__)
    
      ! Exponent for the surface temperature deviation predictor
      ! --------------------------------------------------------
      CALL MAPL_GetResource(STATE,x5Power,'X5_EXPONENT:',DEFAULT=3.00, __RC__)
    
      ! Threshold temperatures
      ! ----------------------
      CALL MAPL_GetResource(STATE,sfcTLimit,'SFC_T_LIMIT:',DEFAULT=273.0, __RC__)
      CALL MAPL_GetResource(STATE,airTLimit,'AIR_T_LIMIT:',DEFAULT=263.0, __RC__)
    
      ! Cloud-top pressure limiter
      ! --------------------------
      CALL MAPL_GetResource(STATE,hPaCldTop,'CLOUD_TOP_LIMIT:',DEFAULT=500., __RC__)
    
      ! Layer depths [m]
      ! ----------------
      ALLOCATE(dZ(nc,lm),STAT=STATUS)
      VERIFY_(STATUS)
      dZ = zle(:,0:lm-1)-zle(:,1:lm)
    
      ! Pressure at mid-layer [Pa]
      ! --------------------------
      ALLOCATE(p(nc,lm),STAT=STATUS)
      VERIFY_(STATUS)
      p = (ple(:,1:lm)+ple(:,0:lm-1))*0.50
    
      ! Temperature at mid-layer [K]
      ! ----------------------------
      ALLOCATE(T(nc,lm),STAT=STATUS)
      VERIFY_(STATUS)
      T = TH*((p*1.00E-05)**(MAPL_RGAS/MAPL_CP))
    
      ! Reset CNV_TOPP's MAPL_UNDEFs to zeroes
      ! --------------------------------------
      ALLOCATE(cnv_topp(nc),STAT=STATUS)
      WHERE(CCTP == MAPL_UNDEF)
         cnv_topp = 0.00
      ELSEWHERE
         cnv_topp = CCTP
      END WHERE
    
      ! Set weak/no convection mask
      ! ---------------------------
      ALLOCATE(weakCnvMask(nc),STAT=STATUS)
      VERIFY_(STATUS)
      weakCnvMask = 0
      WHERE(cn_prcp == 0.00 .OR. cnv_topp >= hPaCldTop*100.00 .OR. CAPE >= MAPL_UNDEF) weakCnvMask = 1
    
      ! Convective cloud top mask
      ! -------------------------
      ALLOCATE(cloudTopMask(nc,lm),STAT=STATUS)
      VERIFY_(STATUS)
      cloudTopMask = 0
      DO k = 1,lm
         WHERE(ple(1:nc,k) > cnv_topp(1:nc) .AND. cnv_topp(1:nc) > 0.00) cloudTopMask(1:nc,k) = 1
      END DO
    
      ! Cloud top distance above ground [m]
      ! -----------------------------------
      ALLOCATE(cloudTopAG(nc),STAT=STATUS)
      VERIFY_(STATUS)
      cloudTopAG = 0.00
      DO i = 1,nc
         n = SUM(cloudTopMask(i,1:lm))
         IF(n > 0) cloudTopAG(i) = SUM(dZ(i,lm-n+1:lm))
      END DO
    
      ! X1: Cold cloud depth: Vertical extent [km] where T < airTLimit and p > cnv_topp
      ! -------------------------------------------------------------------------------
      ALLOCATE(x1(nc),STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(mask(nc,lm),STAT=STATUS)
      VERIFY_(STATUS)
    
      mask = 0
      WHERE(T < airTLimit .AND. cloudTopMask == 1) mask = 1
    
      x1 = 0.00
      DO i = 1,nc
         DO k = 1,lm
            IF(mask(i,k) == 1) x1(i) = x1(i)+dZ(i,k)*0.001
         END DO
      END DO
      WHERE(weakCnvMask == 1) x1 = 0.00
      x1 = x1/x1Divisor
    
      ! X4: Integrated convective mass flux
      ! -----------------------------------
      ALLOCATE(x4(nc),STAT=STATUS)
      VERIFY_(STATUS)
      x4 = 0.00
      DO i = 1,nc
         DO k = 1,lm
            IF(mask(i,k) == 1) x4(i) = x4(i)+cnv_mfc(i,k)*dZ(i,k)
         END DO
      END DO
      WHERE(weakCnvMask == 1) x4 = 0.00
      x4 = x4/x4Divisor
    
      ! X5: Surface temperature deviation from sfcTLimit, positive only.
      ! Note: UNDEF TS test retains the ability to boot-strap moist_import_rst.
      ! -----------------------------------------------------------------------
      ALLOCATE(x5(nc),STAT=STATUS)
      VERIFY_(STATUS)
      WHERE(TS == MAPL_UNDEF)
         x5 = 0.00
      ELSEWHERE
         x5 = TS-sfcTLimit
      END WHERE
      WHERE(weakCnvMask == 1) x5 = 0.00
      WHERE(x5 < 0.00) x5 = 0.00
      x5 = x5**x5Power/x5Divisor

      ! X2: Total cloud depth [km]
      ! --------------------------
      ALLOCATE(x2(nc),STAT=STATUS)
      VERIFY_(STATUS)
      x2 = cloudTopAG*0.001
      WHERE(weakCnvMask == 1) x2 = 0.00
      x2 = x2/x2Divisor
    
      ! X3: CAPE
      ! --------
      ALLOCATE(x3(nc),STAT=STATUS)
      VERIFY_(STATUS)
      x3 = CAPE
      WHERE(weakCnvMask == 1) x3 = 0.00
      x3 = x3/x3Divisor

      ! Polynomial fit [units: km^{-2} s^{-1}] and individual
      ! terms including marine and continental discrimination
      ! -----------------------------------------------------
      WHERE(frOcean >= 0.01)
         strokeRate = (a0m + a1m*x1 + a2m*x2 + a3m*x3 + a4m*x4 + a5m*x5)/86400.00
         A1X1 = a1m*x1/86400.00
         A2X2 = a2m*x2/86400.00
         A3X3 = a3m*x3/86400.00
         A4X4 = a4m*x4/86400.00
         A5X5 = a5m*x5/86400.00
      ELSEWHERE
         strokeRate = (a0c + a1c*x1 + a2c*x2 + a3c*x3 + a4c*x4 + a5c*x5)/86400.00
         A1X1 = a1c*x1/86400.00
         A2X2 = a2c*x2/86400.00
         A3X3 = a3c*x3/86400.00
         A4X4 = a4c*x4/86400.00
         A5X5 = a5c*x5/86400.00
      END WHERE

      ! Eliminate negatives
      ! -------------------
      WHERE(strokeRate < 0.00) strokeRate = 0.00

      ! Set rate to zero where any of x1 through x5 are zero
      ! ----------------------------------------------------
      WHERE(x1 == 0.00) strokeRate = 0.00
      WHERE(x2 == 0.00) strokeRate = 0.00
      WHERE(x3 == 0.00) strokeRate = 0.00
      WHERE(x4 == 0.00) strokeRate = 0.00
      WHERE(x5 == 0.00) strokeRate = 0.00

      ! Clean up
      ! --------
      DEALLOCATE(x1,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(x2,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(x3,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(x4,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(x5,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(cnv_topp,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(dZ,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(p,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(T,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(cloudTopAG,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(mask,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(cloudTopMask,STAT=STATUS)
      VERIFY_(STATUS)
      DEALLOCATE(weakCnvMask,STAT=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine computeFlashRate
!EOC
!-----------------------------------------------------------------------
!BOP
!
      subroutine computeCAPE (TH, Q, PLE, CAPE, BUOY, IM, JM, LM)

      use GEOS_UtilsMod

      integer,                     intent(in)  :: IM,JM,LM
      real, dimension(IM,JM,LM),   intent(in)  :: TH  ! potential temperature
      real, dimension(IM,JM,LM),   intent(in)  :: Q   ! specific humidity
      real, dimension(IM,JM,0:LM), intent(in)  :: PLE   ! pressure
      real, dimension(IM,JM),      intent(out) :: CAPE
      real, dimension(IM,JM,LM),   intent(out) :: BUOY
!EOP
!-----------------------------------------------------------------------
!BOC
      integer                         :: L
      real,    dimension(IM,JM,  LM)  :: DQS, QSS, PLO, TEMP, PK, DM, DP, DZET
      real,    dimension(IM,JM,  LM)  :: ZLO
      real,    dimension(IM,JM,0:LM)  :: ZLE
      real,    dimension(IM,JM,0:LM)    :: CNV_PLE
      real,    dimension(IM,JM,0:LM)    :: PKE
      real,    dimension(IM,JM   )  :: HC
      !logical, dimension(IM,JM   )  :: UNSTABLE

      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      DM       = DP*(1./MAPL_GRAV)
      TEMP     = TH*PK
      DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

      ! From BUYOANCY

       BUOY(:,:,LM) = 0.0

       HC  =  TEMP(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

       do L=LM-1,1,-1
          BUOY(:,:,L) = HC - (TEMP(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QSS(:,:,L))
          BUOY(:,:,L) = BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*TEMP(:,:,L) )
       enddo

       BUOY = MAPL_GRAV*BUOY

! New Formulation 27March2017
       DZET(:,:,1:LM) = TH(:,:,1:LM) * (PKE(:,:,1:LM) - PKE(:,:,0:LM-1)) * MAPL_CP/MAPL_GRAV

       CAPE = 0.
 
       do L=1,LM-1
          where(BUOY(:,:,L)>0.)
              CAPE = CAPE + BUOY(:,:,L)*DZET(:,:,L)
          end where
       end do

! Old formulation
!       do L=1,LM-1
!         where(BUOY(:,:,L)>0.)
!            CAPE = CAPE + BUOY(:,:,L)*DM(:,:,L)
!         end where
!      end do
!
!       UNSTABLE = CAPE > 0.0
!
!       where(.not.UNSTABLE)
!          CAPE=MAPL_UNDEF
!       end where

      
      return

      end subroutine computeCAPE
!EOC
!-----------------------------------------------------------------------
!BOP
!
      subroutine compute_ZLE_RH2 (ZLE, RH2, TH, Q, PLE, IM, JM, LM)

      use GEOS_UtilsMod

      integer,                     intent(in)  :: IM,JM,LM
      real, dimension(IM,JM,LM),   intent(in)  :: TH  ! potential temperature
      real, dimension(IM,JM,LM),   intent(in)  :: Q   ! specific humidity
      real, dimension(IM,JM,0:LM), intent(in)  :: PLE   ! pressure
      real, dimension(IM,JM,LM),   intent(out) :: RH2   ! relative humidity
      real, dimension(IM,JM,0:LM), intent(out) :: ZLE   ! geopotential height
!EOP
!-----------------------------------------------------------------------
!BOC
      integer                         :: L
      real,    dimension(IM,JM,  LM)  :: PLO, PK
      real,    dimension(IM,JM,  LM)  :: ZLO
      real,    dimension(IM,JM,0:LM)  :: CNV_PLE
      real,    dimension(IM,JM,0:LM)  :: PKE

      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)

      ZLE(:,:,LM) = 0.0
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

      RH2     = max(MIN( Q/GEOS_QSAT (TH*PK, PLO) , 1.02 ),0.0)

      return

      end subroutine compute_ZLE_RH2
!EOC
!-----------------------------------------------------------------------
      end module GEOS_ctmEnvGridComp
