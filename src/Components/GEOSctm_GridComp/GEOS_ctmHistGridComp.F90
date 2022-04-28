#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOS_ctmhistGridComp -- Produces meteorological forcings for HISTORY
!
! !INTERFACE:
!
      module GEOS_ctmHistGridComp
!
! !USES:
      use ESMF
      use MAPL

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
      type T_CTMhist_STATE
         private
         logical                    :: output_forcingData  = .FALSE. ! Export variables to HISTORY?
         character(len=ESMF_MAXSTR) :: metType                       ! MERRA2 or MERRA1 or FPIT or FP
      end type T_CTMhist_STATE

      type CTMhist_WRAP
         type (T_CTMhist_STATE), pointer :: PTR
      end type CTMhist_WRAP
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
      type (T_CTMhist_STATE), pointer :: state
      type (CTMhist_WRAP)             :: wrap

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

      call ESMF_UserCompSetInternalState ( GC, 'CTMhist', wrap, STATUS )
      VERIFY_(STATUS)

      configFile = ESMF_ConfigCreate(__RC__)

      call ESMF_ConfigLoadFile(configFile, TRIM(rcfilen), __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%output_forcingData,     &
                             Default  = .FALSE.,                       &
                             Label    = "output_forcingData:",  __RC__ )

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

!      call MAPL_AddImportSpec(GC,                              &
!           SHORT_NAME        = 'AREA',                         &
!           LONG_NAME         = 'agrid_cell_area',              &
!           UNITS             = 'm+2',                          &
!           DIMS              = MAPL_DimsHorzOnly,              &
!           VLOCATION         = MAPL_VLocationNone,    __RC__)

      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME = 'PHIS',                                &
           LONG_NAME  = 'surface_geopotential_height',         &
           UNITS      = 'm2 s-2',                              &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME = 'SLP',                                 &
           LONG_NAME  = 'sea_level_pressure',                  &
           UNITS      = 'Pa',                                  &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME = 'TROPP',                               &
           LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate', &
           UNITS      = 'Pa',                                  &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddImportSpec ( gc,                            &
           SHORT_NAME = 'QLTOT',                               &
           LONG_NAME  = 'mass_fraction_of_cloud_liquid_water', &
           UNITS      = 'kg kg-1',                             &
           DIMS       = MAPL_DimsHorzVert,                     &
           VLOCATION  = MAPL_VLocationCenter,          __RC__  )

      call MAPL_AddImportSpec ( gc,                            &
           SHORT_NAME = 'QITOT',                               &
           LONG_NAME  = 'mass_fraction_of_cloud_ice_water',    &
           UNITS      = 'kg kg-1',                             &
           DIMS       = MAPL_DimsHorzVert,                     &
           VLOCATION  = MAPL_VLocationCenter,          __RC__  )

      call MAPL_AddImportSpec(GC,                             &
          SHORT_NAME = 'TS',                                  &
          LONG_NAME  = 'surface temperature',                 &
          UNITS      = 'K',                                   &
          DIMS       = MAPL_DimsHorzOnly,                     &
          VLOCATION  = MAPL_VLocationNone,            __RC__  )

      call MAPL_AddImportSpec(GC,                             &
           SHORT_NAME = 'Q',                                  &
           LONG_NAME  = 'specific_humidity',                  &
           UNITS      = 'kg kg-1',                            &
           DIMS       = MAPL_DimsHorzVert,                    &
           VLOCATION  = MAPL_VLocationCenter,         __RC__  )

      ! This comes from CTM Env but not from external data files
      call MAPL_AddImportSpec ( gc,                                  &
           SHORT_NAME = 'PLE',                                       &
           LONG_NAME  = 'pressure_at_layer_edges',                   &
           UNITS      = 'Pa',                                        &
           DIMS       = MAPL_DimsHorzVert,                           &
           VLOCATION  = MAPL_VLocationEdge,             __RC__  )

! !EXPORT STATE:

      ! Default exports
      !----------------
      call MAPL_AddExportSpec(GC,                              &
           SHORT_NAME = 'PS',                                  &
           LONG_NAME  = 'surface_pressure',                    &
           UNITS      = 'Pa',                                  &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddExportSpec(GC,                              &
           SHORT_NAME = 'PHIS',                                &
           LONG_NAME  = 'surface_geopotential_height',         &
           UNITS      = 'm2 s-2',                              &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddExportSpec(GC,                              &
           SHORT_NAME = 'TROPP',                               &
           LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate', &
           UNITS      = 'Pa',                                  &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddExportSpec(GC,                              &
           SHORT_NAME = 'SLP',                                 &
           LONG_NAME  = 'sea_level_pressure',                  &
           UNITS      = 'Pa',                                  &
           DIMS       = MAPL_DimsHorzOnly,                     &
           VLOCATION  = MAPL_VLocationNone,              __RC__)

      call MAPL_AddExportSpec ( gc,                            &
           SHORT_NAME = 'QLTOT',                               &
           LONG_NAME  = 'mass_fraction_of_cloud_liquid_water', &
           UNITS      = 'kg kg-1',                             &
           DIMS       = MAPL_DimsHorzVert,                     &
           VLOCATION  = MAPL_VLocationCenter,          __RC__  )

      call MAPL_AddExportSpec ( gc,                            &
           SHORT_NAME = 'QITOT',                               &
           LONG_NAME  = 'mass_fraction_of_cloud_ice_water',    &
           UNITS      = 'kg kg-1',                             &
           DIMS       = MAPL_DimsHorzVert,                     &
           VLOCATION  = MAPL_VLocationCenter,          __RC__  )

      call MAPL_AddExportSpec(GC,                              &
           SHORT_NAME = 'TS',                                  &     
           LONG_NAME  = 'surface temperature',                 &     
           UNITS      = 'K',                                   &     
           DIMS       = MAPL_DimsHorzOnly,                     &     
           VLOCATION  = MAPL_VLocationNone,            __RC__  )

      call MAPL_AddExportSpec(GC,                             &
           SHORT_NAME = 'Q',                                  &      
           LONG_NAME  = 'specific_humidity',                  &      
           UNITS      = 'kg kg-1',                            &      
           DIMS       = MAPL_DimsHorzVert,                    &      
           VLOCATION  = MAPL_VLocationCenter,         __RC__  )

      call MAPL_AddExportSpec ( gc,                           &
           SHORT_NAME = 'DELP',                               &
           LONG_NAME  = 'pressure_thickness',                 &
           UNITS      = 'Pa',                                 &
           DIMS       = MAPL_DimsHorzVert,                    &
           VLOCATION  = MAPL_VLocationCenter,          __RC__ )

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
      type (T_CTMhist_STATE), pointer :: CTMhist_STATE
      type (CTMhist_WRAP)             :: WRAP
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
      call ESMF_UserCompGetInternalState(gc, 'CTMhist', WRAP, STATUS)
      VERIFY_(STATUS)

      CTMhist_STATE => WRAP%PTR

      IF (MAPL_AM_I_ROOT()) THEN
         PRINT*
         PRINT*, '<---------------------------------------------->'
         PRINT*, TRIM(Iam)//' output_forcingData:     ', CTMhist_STATE%output_forcingData
         PRINT*, '<---------------------------------------------->'
         PRINT*
      END IF


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
      type (T_CTMhist_STATE), pointer  :: CTMhist_STATE
      type (CTMhist_wrap)              :: WRAP

      real, pointer, dimension(:,:)   ::  cellArea => null()


      real,     pointer, dimension(:,:,:) ::      Qimp => null()
      real,     pointer, dimension(:,:,:) ::      Qexp => null()
      real,     pointer, dimension(:,:)   ::     PSimp => null()
      real,     pointer, dimension(:,:)   ::     PSexp => null()
      real,     pointer, dimension(:,:)   ::     TSimp => null()
      real,     pointer, dimension(:,:)   ::     TSexp => null()
      real,     pointer, dimension(:,:)   ::    SLPimp => null()
      real,     pointer, dimension(:,:)   ::    SLPexp => null()
      real,     pointer, dimension(:,:)   ::   PHISimp => null()
      real,     pointer, dimension(:,:)   ::   PHISexp => null()
      real,     pointer, dimension(:,:,:) ::   DELPexp => null()
      real,     pointer, dimension(:,:)   ::  TROPPimp => null()
      real,     pointer, dimension(:,:)   ::  TROPPexp => null()
      real,     pointer, dimension(:,:,:) ::  QITOTimp => null()
      real,     pointer, dimension(:,:,:) ::  QITOTexp => null()
      real,     pointer, dimension(:,:,:) ::  QLTOTimp => null()
      real,     pointer, dimension(:,:,:) ::  QLTOTexp => null()
      real,     pointer, dimension(:,:,:) ::       PLE => null()

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
      call ESMF_UserCompGetInternalState(gc, 'CTMhist', WRAP, STATUS)
      VERIFY_(STATUS)

      CTMhist_STATE => WRAP%PTR

      ! Get the time-step
      ! -----------------------
      call MAPL_GetResource( ggState, ndt, 'RUN_DT:', default=0, __RC__ )
      DT = ndt

      call MAPL_GetPointer ( IMPORT,    PSimp,     'PS',  __RC__ )
      call MAPL_GetPointer ( EXPORT,    PSexp,     'PS', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(PSimp).AND.ASSOCIATED(PSexp)) PSexp = PSimp

      call MAPL_GetPointer ( IMPORT,    TSimp,     'TS',  __RC__ )
      call MAPL_GetPointer ( EXPORT,    TSexp,     'TS', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(TSimp).AND.ASSOCIATED(TSexp)) TSexp = TSimp

      call MAPL_GetPointer ( IMPORT,   SLPimp,    'SLP',  __RC__ )
      call MAPL_GetPointer ( EXPORT,   SLPexp,    'SLP', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(SLPimp).AND.ASSOCIATED(SLPexp)) SLPexp = SLPimp

      call MAPL_GetPointer ( IMPORT,  PHISimp,   'PHIS', __RC__ )
      call MAPL_GetPointer ( EXPORT,  PHISexp,   'PHIS', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(PHISimp).AND.ASSOCIATED(PHISexp)) PHISexp = PHISimp

      call MAPL_GetPointer ( IMPORT,  TROPPimp,   'TROPP', __RC__ )
      call MAPL_GetPointer ( EXPORT,  TROPPexp,   'TROPP', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(TROPPimp).AND.ASSOCIATED(TROPPexp)) TROPPexp = TROPPimp

      call MAPL_GetPointer ( IMPORT,  QITOTimp,  'QITOT', __RC__ )
      call MAPL_GetPointer ( EXPORT,  QITOTexp,  'QITOT', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(QITOTimp).AND.ASSOCIATED(QITOTexp)) QITOTexp = QITOTimp

      call MAPL_GetPointer ( IMPORT,  QLTOTimp,  'QLTOT', __RC__ )
      call MAPL_GetPointer ( EXPORT,  QLTOTexp,  'QLTOT', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(QLTOTimp).AND.ASSOCIATED(QLTOTexp)) QLTOTexp = QLTOTimp

      call MAPL_GetPointer ( IMPORT,     Qimp,     'Q',  __RC__ )
      call MAPL_GetPointer ( EXPORT,     Qexp,     'Q', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(Qimp).AND.ASSOCIATED(Qexp)) Qexp = Qimp

      call MAPL_GetPointer ( IMPORT,  PLE,       'PLE', __RC__ )
      call MAPL_GetPointer ( EXPORT,  DELPexp,   'DELP', ALLOC=.TRUE., __RC__ )
      IF (ASSOCIATED(DELPexp) .AND. ASSOCIATED(PLE)) THEN
         call compute_DELP()
      ENDIF

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
         subroutine compute_PLE()
           INTEGER  :: L
           real(r8), allocatable             :: AK(:)
           real(r8), allocatable             :: BK(:)

           allocate(AK(LM+1),stat=rc)
           VERIFY_(rc)
           allocate(BK(LM+1),stat=rc)
           VERIFY_(rc)

           call ESMF_AttributeGet(esmfGrid, name="GridAK", valuelist=AK, __RC__)
           call ESMF_AttributeGet(esmfGrid, name="GridBK", valuelist=BK, __RC__)

           DO L = 1, LM+1
              PLE(:,:,L) = AK(L) + BK(L)*PSimp(:,:)
           END DO

           deallocate(AK, BK)

         end subroutine compute_PLE
         !
         !-------------------------------------------------------
         !
         subroutine compute_DELP()
           INTEGER  :: L

           DO L = 1, SIZE(PLE, 3)-1
              DELPexp(:,:,L) = PLE(:,:,L+1) - PLE(:,:,L)
           END DO

         end subroutine compute_DELP

      end subroutine Run
!------------------------------------------------------------------------------
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
      end module GEOS_ctmHistGridComp
