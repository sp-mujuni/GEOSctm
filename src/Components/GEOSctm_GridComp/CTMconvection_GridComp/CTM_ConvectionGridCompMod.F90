#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: CTM_ConvectionGridCompMod - The Convection Grid Component
!
! !INTERFACE:
!
   MODULE CTM_ConvectionGridCompMod
!
! !USES:
!
      USE ESMF
      USE MAPL_Mod
      USE GmiConvectionMethod_mod             ! GMI     Convection component
      USE GenericConvectionMethod_mod         ! Generic Convection component
      USE m_chars, ONLY : uppercase
      use CTM_rasCalculationsMod, only : INIT_RASPARAMS, DO_RAS, RASPARAM_Type
      USE Chem_UtilMod, only : pmaxmin

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices
!
! !DESCRIPTION: 
!
!  {\tt CTM\_ConvectionGridCompMod} is a ESMF gridded component for the 
!  Convection package.
!
! !REVISION HISTORY:
!
!  07Dec2012  Kouatchou  Created the Convection stub.
!
!EOP
!-------------------------------------------------------------------------
      TYPE Convection_State
         PRIVATE
         TYPE(gmiConvection_GridComp),  POINTER :: gmiCONV    => null()
         TYPE(genConvection_GridComp),  POINTER :: genCONV    => null()
      END TYPE Convection_State
    
      TYPE Convection_WRAP
         TYPE (Convection_State), pointer :: PTR => null()
      END TYPE Convection_WRAP

      integer :: convecType ! 1:  Generic Convection (only convective transport)
                            ! 2:  GMI convection

      TYPE(RASPARAM_Type) :: RASPARAMS
      logical             :: enable_rasCalculations = .FALSE.

!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for Convection Grid Component
!
! !INTERFACE:

   SUBROUTINE SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  07Dec2012  Kouatchou  First crack.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
      character(len=ESMF_MAXSTR)      :: COMP_NAME
      CHARACTER(LEN=ESMF_MAXSTR)      :: FRIENDLIES, providerName
      character(len=ESMF_MAXSTR)      :: rcfilen = 'CTM_GridComp.rc'

      type (ESMF_Config)              :: CF
      type (Convection_State), pointer:: state   ! internal, that is
      type (Convection_wrap)          :: wrap
      type (ESMF_Config)              :: convConfigFile

      integer                         :: STATUS
      integer                         :: m, n, i_XX, j_XX

      LOGICAL                         :: searchForImports
      CHARACTER(LEN=2)                :: leadChars
      CHARACTER(LEN=ESMF_MAXSTR)      :: name

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
      VERIFY_(STATUS)
      Iam = TRIM(COMP_NAME)//"::SetServices"

      !   Wrap internal state for storing in GC; rename legacyState
      !   -------------------------------------
      allocate ( state, stat=STATUS )
      VERIFY_(STATUS)
      wrap%ptr => state

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------


      IF (MAPL_AM_I_ROOT()) THEN
         PRINT *, TRIM(Iam)//': ACTIVE'
      END IF


      !   Set the Initialize, Run, Finalize entry points
      !   ----------------------------------------------

      CALL MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize_, __RC__ )
      CALL MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN,                Run_, __RC__ )
      CALL MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_FINALIZE,      Finalize_, __RC__ )

      !   Store internal state in GC
      !   --------------------------
      CALL ESMF_UserCompSetInternalState ( GC, 'Convection_state', wrap, STATUS )
      VERIFY_(STATUS)

      convConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(convConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(convConfigFile, enable_rasCalculations, &
                          Default  = .FALSE.,                          &
                          Label    = "enable_rasCalculations:", __RC__ )

! ========================== IMPORT STATE =========================

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'T',  &
           LONG_NAME          = 'air_temperature',  &
           UNITS              = 'K', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationCenter,    __RC__ )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'Q',  &
           LONG_NAME          = 'specific_humidity',  &
           UNITS              = 'kg kg-1', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationCenter,    __RC__ )


      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'PLE',  &
           LONG_NAME          = 'air_pressure',  &
           UNITS              = 'Pa', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationEdge,    __RC__ )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'MASS',  &
           LONG_NAME          = 'total_mass',  &
           UNITS              = 'kg', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationCenter,    __RC__ )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'ZLE',  &
           LONG_NAME          = 'geopotential_height',  &
           UNITS              = 'm', &
           DIMS               = MAPL_DimsHorzVert,    &
           VLOCATION          = MAPL_VLocationEdge,    __RC__ )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'ZPBL',  &
           LONG_NAME          = 'planetary_boundary_layer_height', &
           UNITS              = 'm',                               &
           DIMS               = MAPL_DimsHorzOnly,                 &
           VLOCATION          = MAPL_VLocationNone,         __RC__ )

      call MAPL_AddImportSpec(GC,                                  &
           SHORT_NAME         = 'LWI',                             &
           LONG_NAME          = 'land-ocean-ice_mask',             &
           UNITS              = '1',                               &
           DIMS               = MAPL_DimsHorzOnly,                 &
           VLOCATION          = MAPL_VLocationNone,         __RC__ )

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'AREA',  &
           LONG_NAME          = 'agrid_cell_area',  &
           UNITS              = 'm^2', &
           DIMS               = MAPL_DimsHorzOnly,    &
           VLOCATION          = MAPL_VLocationNone,    __RC__ )

      call MAPL_AddImportSpec(GC,                                  &
           SHORT_NAME         = 'ConvTR',                            &
           LONG_NAME          = 'convected_quantities',              &
           UNITS              = 'X',                                 &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
           DATATYPE           = MAPL_BundleItem,                     &
           RESTART            = MAPL_RestartOptional,                __RC__ )

      IF (enable_rasCalculations) THEN
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

         call MAPL_AddImportSpec(GC,                            &
              SHORT_NAME = 'TS',                                &
              LONG_NAME  = 'surface temperature',               &
              UNITS      = 'K',                                 &
              DIMS       = MAPL_DimsHorzOnly,                   &
              VLOCATION  = MAPL_VLocationNone,           __RC__ )
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

! ========================== EXPORT STATE =========================

      IF (enable_rasCalculations) THEN
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
      ENDIF


!#include "convTendency_ExportSpec.h"


!   Set the Profiling timers
!   ------------------------
    call MAPL_TimerAdd ( GC, name = "INITIALIZE", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_TimerAdd ( GC, name = "RUN", RC=STATUS )
    VERIFY_(STATUS)

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  END SUBROUTINE SetServices
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize GMI Convection
!
! !INTERFACE:
!
      SUBROUTINE Initialize_ ( gc, impConv, expConv, clock, rc )

! !INPUT PARAMETERS:

      type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

      type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
      type(ESMF_State), intent(inout) :: impConv     ! Import State
      type(ESMF_State), intent(inout) :: expConv     ! Export State
      integer, intent(out) ::  rc                    ! Error return code:
                                                     !  0 - all is well
                                                     !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
      character(len=ESMF_MAXSTR)      :: rcfilen = 'CTM_GridComp.rc'
      character(len=ESMF_MAXSTR)      :: COMP_NAME

      type(ESMF_Config)               :: CF
      type(ESMF_Grid)                 :: esmfGrid
      type(MAPL_MetaComp), pointer    :: ggState	   ! GEOS Generic State
      type(ESMF_State)                :: internal
      type(MAPL_VarSpec), pointer     :: InternalSpec(:)
      type (ESMF_Config)              :: convConfigFile
 
      integer                         :: STATUS
      integer                         :: nymd, nhms  ! time of day
      real                            :: cdt         ! chemistry timestep (secs)
      integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
      integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
      integer                         :: km                  ! dist grid indices
      integer                         :: dims(3), k, l, n

      type(gmiConvection_GridComp), pointer     :: gmiCONV      ! Grid Component
      type(genConvection_GridComp), pointer     :: genCONV      ! Grid Component

      rc = 0

      !  Get my name and set-up traceback handle
      !  ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
      VERIFY_(STATUS)
      Iam = TRIM(COMP_NAME)//"::Initialize"

      !  Get my internal MAPL_Generic state
      !  -----------------------------------
      call MAPL_GetObjectFromGC ( GC, ggState, RC=STATUS)
      VERIFY_(STATUS)

      !  Initialize GEOS Generic
      !  ------------------------
      call MAPL_GenericInitialize ( gc, impConv, expConv, clock,  RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_TimerOn(ggSTATE,"TOTAL")
      call MAPL_TimerOn(ggSTATE,"INITIALIZE")

      convConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(convConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      ! ------------------------------
      ! convecType
      !   1:  Generic Convection (only convective transport)
      !   2:  GMI convection
      !------ ------------------------

      call ESMF_ConfigGetAttribute(convConfigFile, convecType, &
                         label   = "convecType:",&
                         default = 1, rc=STATUS )
      VERIFY_(STATUS)


      !  Get parameters from gc and clock
      !  --------------------------------
      call extract_ ( gc, clock, genCONV, gmiCONV, nymd, nhms, cdt, STATUS )
      VERIFY_(STATUS)

      call ESMF_GridCompGet ( GC, GRID=esmfGrid, rc=STATUS)
      VERIFY_(STATUS)

      call MAPL_GridGet (esmfGrid, globalCellCountPerDim=DIMS, RC=STATUS)
      VERIFY_(STATUS)

      im = dims(1)
      jm = dims(2)

      call ESMF_GridGet(esmfGrid, localDE=0,staggerloc=ESMF_STAGGERLOC_CENTER, &
                       computationalCount=DIMS, RC=STATUS)
      VERIFY_(STATUS)

      ! Local sizes of three dimensions
      !--------------------------------
      i2 = dims(1)
      j2 = dims(2)
      km = dims(3)

      !  Call initialize
      !  ---------------
      if (convecType == 1) then
         genCONV%i1 = i1
         genCONV%i2 = i2
         genCONV%im = im
         genCONV%j1 = j1
         genCONV%j2 = j2
         genCONV%jm = jm
         genCONV%km = km

         call initializeGenericConvection ( genCONV, impConv, expConv, nymd, nhms, &
                        esmfGrid, cdt, STATUS )
         VERIFY_(STATUS)
      elseif (convecType == 2) then
         gmiCONV%i1 = i1
         gmiCONV%i2 = i2
         gmiCONV%im = im
         gmiCONV%j1 = j1
         gmiCONV%j2 = j2
         gmiCONV%jm = jm
         gmiCONV%km = km

         call initializeGmiConvection ( gmiCONV, impConv, expConv, nymd, nhms, &
                        esmfGrid, cdt, STATUS )
         VERIFY_(STATUS)
      end if

      IF (enable_rasCalculations) THEN
         IF (MAPL_AM_I_ROOT()) PRINT*, TRIM(Iam)//': Doing RAS Calculations'
         CALL INIT_RASPARAMS(RASPARAMS)
      ENDIF

      call MAPL_TimerOff(ggSTATE,"INITIALIZE")
      call MAPL_TimerOff(ggSTATE,"TOTAL")

      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE Initialize_
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs Convection
!
! !INTERFACE:
!

      SUBROUTINE Run_ ( gc, impConv, expConv, clock, rc )

! !USES:
     !USE Chem_UtilMod, only : pmaxmin

     implicit NONE

! !INPUT PARAMETERS:

      type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

      type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
      type(ESMF_State), intent(inout) :: impConv     ! Import State
      type(ESMF_State), intent(inout) :: expConv     ! Export State
      integer, intent(out) ::  rc                    ! Error return code:
                                                     !  0 - all is well
                                                     !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  07Dec2012 Kouatchou  First crack.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!  ErrLog Variables
!  ----------------
      character(len=ESMF_MAXSTR)      :: IAm = 'Run'
      integer                         :: STATUS
      character(len=ESMF_MAXSTR)      :: COMP_NAME

      type(gmiConvection_GridComp), pointer     :: gmiCONV       ! Grid Component
      type(genConvection_GridComp), pointer     :: genCONV       ! Grid Component
      integer                         :: nymd, nhms  ! time
      real                            :: cdt         ! chemistry timestep (secs)
      integer                         :: i, iOX, iT2M, k, m, n

      type(ESMF_Config)               :: CF
      type(ESMF_Grid)                 :: esmfGrid        
      type(ESMF_Time)                 :: TIME
      type (MAPL_MetaComp), pointer   :: ggState

      !REAL :: qmin, qmax
      !real, pointer, dimension(:,:,:) ::     CNV_MFC => null()
      !real, pointer, dimension(:,:,:) ::     CNV_MFD => null()

!  Get my name and set-up traceback handle
!  ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, GRID=esmfGrid, RC=STATUS )
      VERIFY_(STATUS)
      Iam = TRIM(COMP_NAME)//"::Run"


!  Get ESMF parameters from gc and clock
!  -----------------------------------------
      call extract_ ( gc, clock, genCONV, gmiCONV, nymd, nhms, & 
                      cdt, rc=status )
      VERIFY_(STATUS)

      IF (enable_rasCalculations) THEN
         CALL runRAS(impConv, expConv, esmfGrid, cdt)

         !call MAPL_GetPointer ( expConv,  CNV_MFC,  'CNV_MFC', __RC__ )
         !call MAPL_GetPointer ( expConv,  CNV_MFD,  'CNV_MFD', __RC__ )

         !CALL pmaxmin('CNV_MFC-:',CNV_MFC, qmin, qmax, size(CNV_MFD,1)*size(CNV_MFD,2), size(CNV_MFD,3)+1, 1. )
         !CALL pmaxmin('CNV_MFD-:',CNV_MFD, qmin, qmax, size(CNV_MFD,1)*size(CNV_MFD,2), size(CNV_MFD,3), 1. )

      ENDIF

!  Run
!  ---
      if (convecType == 1) then
         CALL runGenericConvection ( genCONV, impConv, expConv, nymd, nhms, &
                          cdt, enable_rasCalculations, STATUS )
         VERIFY_(STATUS)
      elseif (convecType == 2) then
         CALL runGmiConvection ( gmiCONV, impConv, expConv, nymd, nhms, &
                          cdt, enable_rasCalculations, STATUS )
         VERIFY_(STATUS)
      end if

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize CTM_ConvectionGridCompMod (ESMF)
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( gc, impConv, expConv, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impConv     ! Import State
   type(ESMF_State), intent(inout) :: expConv     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
!BOC
!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(gmiConvection_GridComp), pointer     :: gmiCONV      ! Grid Component
   type(genConvection_GridComp), pointer     :: genCONV      ! Grid Component
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)

    type(Convection_state), pointer   :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"Finalize"

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   call extract_ ( gc, clock, genCONV, gmiCONV, nymd, nhms, cdt, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Call ESMF version
!  -----------------
      if (convecType == 1) then
         call finalizeGenericConvection ( genCONV )
      elseif (convecType == 2) then
         call finalizeGmiConvection ( gmiCONV )
      end if

!  Finalize MAPL Generic.  Atanas says, "Do not deallocate foreign objects."
!  -------------------------------------------------------------------------
   call MAPL_GenericFinalize ( gc, impConv, expConv, clock,  RC=STATUS )
   VERIFY_(STATUS)

!  Destroy Legacy state
!  --------------------
      if (convecType == 1) then
         deallocate ( state%genCONV, stat = STATUS )
         VERIFY_(STATUS)
      elseif (convecType == 2) then
         deallocate ( state%gmiCONV, stat = STATUS )
         VERIFY_(STATUS)
      end if

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Finalize_
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  runRAS
!
! !INTERFACE:
!
      subroutine runRAS(impConv, expConv, esmfGrid, cdt)
!
! !USES:
      !USE Chem_UtilMod, only : pmaxmin
      use, intrinsic :: iso_fortran_env, only: REAL64

     implicit NONE
!
! !INPUT PARAMETERS:
      real,             intent(in   ) :: cdt
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Grid),  intent(inout) :: esmfGrid
      type(ESMF_State), intent(inout) :: impConv     ! Import State
      type(ESMF_State), intent(inout) :: expConv     ! Export State
!EOP
!-------------------------------------------------------------------------
!BOC

      character(len=ESMF_MAXSTR)      :: IAm = 'runRAS'
      integer                         :: RC, STATUS, IM, JM, LM
      real, pointer, dimension(:,:)   ::      FRLAND => null()
      real, POINTER, dimension(:,:)   ::        PBLH => null() 
      real, pointer, dimension(:,:)   ::          TS => null()
      real, pointer, dimension(:,:,:) ::     CNV_MFC => null()
      real, pointer, dimension(:,:,:) ::     CNV_MFD => null()
      real, pointer, dimension(:,:,:) ::           T => null()
      real, pointer, dimension(:,:,:) ::          TH => null()
      real, pointer, dimension(:,:,:) ::          KH => null()
      real, pointer, dimension(:,:,:) ::           Q => null()
      real, pointer, dimension(:,:,:) ::         PLE => null()

      real, pointer, dimension(:)     ::        PREF => null()
      real, pointer, dimension(:,:,:) ::  CNV_MFCras => null()
      real, pointer, dimension(:,:,:) ::  CNV_MFDras => null()

      real(REAL64), pointer, dimension(:,:) :: LATS    => null()
      real(REAL64), pointer, dimension(:,:) :: LONS    => null()
      real(REAL64), allocatable             :: AK(:)
      real(REAL64), allocatable             :: BK(:)
      !REAL :: qmin, qmax

         call MAPL_GetPointer ( impConv,       Q,         'Q',  __RC__ )
         call MAPL_GetPointer ( impConv,       T,         'T',  __RC__ )
         call MAPL_GetPointer ( impConv,      KH,        'KH',  __RC__ )
         call MAPL_GetPointer ( impConv,      TS,        'TS',  __RC__ )
         call MAPL_GetPointer ( impConv,     PLE,       'PLE',  __RC__ )
         call MAPL_GetPointer ( impConv,    PBLH,      'ZPBL',  __RC__ )
         call MAPL_GetPointer ( impConv,  FRLAND,    'FRLAND',  __RC__ )

         IM = SIZE(T, 1)
         JM = SIZE(T, 2)
         LM = SIZE(T, 3)

         ! Get the AK and BK
         allocate(AK(LM+1),stat=status)
         VERIFY_(STATUS)
         allocate(BK(LM+1),stat=status)
         VERIFY_(STATUS)

         call ESMF_AttributeGet(esmfGrid,name="GridAK",valuelist=AK,rc=status)
         VERIFY_(STATUS)
         call ESMF_AttributeGet(esmfGrid,name="GridBK",valuelist=BK,rc=status)
         VERIFY_(STATUS)

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
         PREF = AK + BK * MAPL_P00

         ALLOCATE(TH        (IM,JM,1:LM))
         ALLOCATE(CNV_MFDras(IM,JM,1:LM))
         ALLOCATE(CNV_MFCras(IM,JM,0:LM))

         ! Compute the potential temperature
         TH(:,:,:) = T(:,:,:)/ &
                     ((0.5*(PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) ))/100000.)**(MAPL_RGAS/MAPL_CP)

!         IF ( MAPL_am_I_root() ) THEN
!            PRINT*, "DT: ", cdt
!            PRINT*
!            PRINT*, "AK: ", AK
!            PRINT*
!            PRINT*, "BK: ", BK
!            PRINT*
!            PRINT*, "PREF: ", PREF
!            PRINT*
!         END IF
!
!         CALL pmaxmin('LONS :',REAL(LONS), qmin, qmax, IM*JM,    1, 1. )
!         CALL pmaxmin('LATS :',REAL(LATS), qmin, qmax, IM*JM,    1, 1. )
!         CALL pmaxmin('KH :',        KH, qmin, qmax, IM*JM, LM+1, 1. )
!         CALL pmaxmin('TH :',        TH, qmin, qmax, IM*JM, LM  , 1. )
!         CALL pmaxmin('TS :',        TS, qmin, qmax, IM*JM,    1, 1. )
!         CALL pmaxmin('Q  :',         Q, qmin, qmax, IM*JM, LM  , 1. )
!         CALL pmaxmin('FRLAND:', FRLAND, qmin, qmax, IM*JM,    1, 1. )
!         CALL pmaxmin('PLE :',      PLE, qmin, qmax, IM*JM, LM+1, 1. )
!         CALL pmaxmin('PBLH :',    PBLH, qmin, qmax, IM*JM,    1, 1. )
!
         ! Do the RAS calculations to obtain the values of CNV_MFD and CNV_MFC
         CALL DO_RAS(RASPARAMS, PREF, TH, PLE, KH, PBLH, Q, &
                     FRLAND, TS, CNV_MFDras, CNV_MFCras,      &
                     REAL(LATS), REAL(LONS), cdt, IM, JM, LM)

         !CALL pmaxmin('CNV_MFC :',CNV_MFCras, qmin, qmax, IM*JM, LM+1, 1. )
         !CALL pmaxmin('CNV_MFD :',CNV_MFDras, qmin, qmax, IM*JM, LM, 1. )

         ! Populate the export state
         call MAPL_GetPointer ( expConv,  CNV_MFC,  'CNV_MFC', ALLOC=.TRUE., __RC__ )
         call MAPL_GetPointer ( expConv,  CNV_MFD,  'CNV_MFD', ALLOC=.TRUE., __RC__ )

         CNV_MFC(:,:,:) = CNV_MFCras(:,:,:)
         CNV_MFD(:,:,:) = CNV_MFDras(:,:,:)

         DEALLOCATE(PREF)
         DEALLOCATE(CNV_MFDras, CNV_MFCras)
         DEALLOCATE(TH)
         deallocate(AK, BK)
      end subroutine runRAS
!EOC
!-------------------------------------------------------------------------
    SUBROUTINE extract_ ( gc, clock, genCONV, gmiCONV, nymd, nhms, cdt, &
                          rc, state )

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_Clock), intent(in)       :: clock
    type(gmiConvection_GridComp), pointer :: gmiCONV
    type(genConvection_GridComp), pointer :: genCONV
    integer, intent(out)               :: nymd, nhms
    real, intent(out)                  :: cdt
    integer, intent(out)               :: rc
    type(Convection_state), pointer, optional   :: state


    type(Convection_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(Convection_Wrap)  :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'extract_'

    rc = 0

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'Convection_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

!   This is likely to be allocated during initialize only
!   -----------------------------------------------------
      if (convecType == 1) then
         if ( .not. associated(myState%genCONV) ) then
              allocate ( myState%genCONV, stat=STATUS )
              VERIFY_(STATUS)
         end if

         genCONV  => myState%genCONV
      elseif (convecType == 2) then
         if ( .not. associated(myState%gmiCONV) ) then
              allocate ( myState%gmiCONV, stat=STATUS )
              VERIFY_(STATUS)
         end if

         gmiCONV  => myState%gmiCONV
      end if

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get time step
!   -------------
    call ESMF_ConfigGetAttribute ( CF, cdt, LABEL="RUN_DT:", RC=STATUS )
    VERIFY_(STATUS)

!   Disable Option for Modified cdt (Chemistry does not operate with modified refresh intervals)
!   --------------------------------------------------------------------------------------------
!   call ESMF_ConfigGetAttribute ( CF, cdt, LABEL="CHEMISTRY_DT:", DEFAULT=cdt, RC=STATUS )
!   VERIFY_(STATUS)

!   Extract nymd, nhms, day of year from clock
!   ------------------------------------------
    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   END SUBROUTINE extract_

 END MODULE CTM_ConvectionGridCompMod
