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
      USE MAPL
      !USE GmiConvectionMethod_mod             ! GMI     Convection component
      USE convectiveTransport_mod
      USE m_chars, ONLY : uppercase
      !USE Chem_UtilMod, only : pmaxmin
      use GmiArrayBundlePointer_mod
      USE GmiESMFrcFileReading_mod

   IMPLICIT NONE

      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
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
      TYPE T_Convection_State
         PRIVATE
         logical :: det_ent      = .FALSE. ! flag for doing detrainment then entrainment
         logical :: do_downdraft = .FALSE. ! flag for doing downdrafts
         integer :: convecType   = 1       ! 1:  Generic Convection (only convective transport)
                                           ! 2:  GMI convection
         integer :: numSpecies
         logical, pointer :: isFixedConcentration(:) => null()
         REAL(KIND=DBL), ALLOCATABLE :: cellArea(:,:)
         logical :: FIRST = .TRUE.
      END TYPE T_Convection_State
    
      TYPE Convection_WRAP
         TYPE (T_Convection_State), pointer :: PTR => null()
      END TYPE Convection_WRAP

      integer :: i1=1, i2, ig=0, im  ! dist grid indices
      integer :: j1=1, j2, jg=0, jm  ! dist grid indices
      integer :: km                  ! dist grid indices
      integer :: k1=1, k2, ivert, ilong

      REAL, PARAMETER :: mwtAir    = 28.9
      REAL, PARAMETER :: rStar     = 8.314E+03
      REAL, PARAMETER :: Pa2hPa    = 0.01
      REAL, PARAMETER :: ToGrPerKg = 1000.00
      REAL, PARAMETER :: secPerDay = 86400.00

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
      type (T_Convection_State), pointer:: state   ! internal, that is
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

      call ESMF_ConfigGetAttribute(convConfigFile, state%det_ent, &
                                   Label = "det_ent:",      &
                                   default=.FALSE., __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, state%do_downdraft, &
                                   Label = "do_downdraft:", &
                                   default=.FALSE., __RC__ )
      ! ------------------------------
      ! convecType
      !   1:  Generic Convection (only convective transport)
      !   2:  GMI convection
      !------ ------------------------

      call ESMF_ConfigGetAttribute(convConfigFile, state%convecType, &
                                   label   = "convecType:",&
                                   default = 1, __RC__ )
      state%FIRST = .TRUE.

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT*," -----> det_ent      = ", state%det_ent
         PRINT*," -----> do_downdraft = ", state%do_downdraft
         PRINT*," -----> convecType   = ", state%convecType
         PRINT *,"Done Reading the Convection Resource File"
      END IF

! ========================== IMPORT STATE =========================

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'T',  &
           LONG_NAME          = 'air_temperature',  &
           UNITS              = 'K', &
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

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME         = 'AREA',  &
           LONG_NAME          = 'agrid_cell_area',  &
           UNITS              = 'm^2', &
           DIMS               = MAPL_DimsHorzOnly,    &
           VLOCATION          = MAPL_VLocationNone,    __RC__ )

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

      call MAPL_AddImportSpec(GC,                                  &
           SHORT_NAME         = 'ConvTR',                            &
           LONG_NAME          = 'convected_quantities',              &
           UNITS              = 'X',                                 &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
           DATATYPE           = MAPL_BundleItem,                     &
           RESTART            = MAPL_RestartOptional,                __RC__ )

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
      character(len=ESMF_MAXSTR)      :: COMP_NAME

      type(ESMF_Config)               :: CF
      type(ESMF_Grid)                 :: esmfGrid
      type(MAPL_MetaComp), pointer    :: ggState	   ! GEOS Generic State
      type(ESMF_State)                :: internal
      type(MAPL_VarSpec), pointer     :: InternalSpec(:)
      type (ESMF_Config)              :: convConfigFile
      type (T_Convection_State), pointer:: conv_state   ! internal, that is
      type (Convection_wrap)          :: wrap
 
      integer                         :: STATUS
      integer                         :: nymd, nhms  ! time of day
      real                            :: cdt         ! chemistry timestep (secs)
      integer                         :: dims(3), k, l, n
      type (ESMF_Field)               :: FIELD
      type (ESMF_Array)               :: ARRAY
      type (ESMF_FieldBundle)         :: ConvTR
      REAL, POINTER, DIMENSION(:,:,:) :: S
      character(len=ESMF_MAXSTR)      :: field_name
      integer                         :: ic
      REAL, POINTER, DIMENSION(:,:) :: gridBoxArea

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

      ! Get my private state from the component
      !----------------------------------------

      call ESMF_UserCompGetInternalState(gc, 'Convection_state', WRAP, STATUS)
      VERIFY_(STATUS)

      conv_state => WRAP%PTR

      !  Get parameters from gc and clock
      !  --------------------------------
      call extract_ ( gc, clock, nymd, nhms, cdt, STATUS )
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
      i2    = dims(1)
      j2    = dims(2)
      km    = dims(3)
      k2    = km
      ivert = km
      ilong = i2-i1+1

      ! Grid box surface area, m^{2}
      ! ----------------------------
      CALL MAPL_GetPointer(impConv, gridBoxArea,    'AREA', __RC__)
      allocate(conv_state%cellArea(i1:i2,j1:j2), STAT=STATUS); VERIFY_(STATUS)
      conv_state%cellArea = gridBoxArea

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

      integer                         :: nymd, nhms  ! time
      real                            :: cdt         ! chemistry timestep (secs)
      integer                         :: i, iOX, iT2M, k, m, n

      type(ESMF_Config)               :: CF
      type(ESMF_Grid)                 :: esmfGrid        
      type(ESMF_Time)                 :: TIME
      type (MAPL_MetaComp), pointer   :: ggState
      type (T_Convection_State), pointer  :: conv_state
      type (Convection_wrap)              :: WRAP

      REAL, POINTER, DIMENSION(:,:)   :: zpbl
      REAL, POINTER, DIMENSION(:,:,:) :: ple, zle, totalMass
      REAL, POINTER, DIMENSION(:,:,:) :: CNV_MFC, CNV_MFD, T

      integer :: ic, kR, ik, is
      REAL, ALLOCATABLE :: pl(:,:,:)

      REAL(KIND=DBL), ALLOCATABLE :: pbl(:,:)
      REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: gridBoxHeight(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: cmf(:,:,:)
      REAL(KIND=DBL), allocatable :: dtrain       (:, :, :)
      REAL(KIND=DBL), allocatable :: eu         (:, :, :)
      REAL(KIND=DBL), allocatable :: ed         (:, :, :)
      REAL(KIND=DBL), allocatable :: md         (:, :, :)
      REAL(KIND=DBL), allocatable :: kel        (:, :, :)
      REAL(KIND=DBL)              :: tdt

      type (ESMF_Field)                   :: FIELD
      type (ESMF_Array)                   :: ARRAY
      type (ESMF_FieldBundle)             :: ConvTR
      REAL, POINTER, DIMENSION(:,:,:)     :: S

      type (t_GmiArrayBundle), pointer :: concentration(:)
      character(len=ESMF_MAXSTR) :: NAME, speciesName
      REAL :: qmin, qmax
!EOP
!-------------------------------------------------------------------------
!BOC

!  Get my name and set-up traceback handle
!  ---------------------------------------
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, GRID=esmfGrid, RC=STATUS )
      VERIFY_(STATUS)
      Iam = TRIM(COMP_NAME)//"::Run"

      ! Get my internal MAPL_Generic state
      !-----------------------------------
      call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )

      call MAPL_TimerOn(ggState,"TOTAL")
      call MAPL_TimerOn(ggState,"RUN")

      ! Get my private state from the component
      !----------------------------------------
      call ESMF_UserCompGetInternalState(gc, 'Convection_state', WRAP, STATUS)
      VERIFY_(STATUS)

      conv_state => WRAP%PTR

!  Get ESMF parameters from gc and clock
!  -----------------------------------------
      call extract_ ( gc, clock, nymd, nhms, cdt, rc=status )
      VERIFY_(STATUS)

!  Run
!  ---
      !---------------------
      ! Convective Transport
      !---------------------

      ! Get the bundles containing the quantities to be diffused, 
      !----------------------------------------------------------

      call ESMF_StateGet(impConv, 'ConvTR' ,    ConvTR,     RC=STATUS)
      VERIFY_(STATUS)

      ! Identify the "fixed" tracers
      !-----------------------------
      IF (conv_state%FIRST) THEN
         conv_state%FIRST = .FALSE.

         call ESMF_FieldBundleGet(ConvTR, fieldCOUNT=conv_state%numSpecies, RC=STATUS)
         VERIFY_(STATUS)

         allocate(conv_state%isFixedConcentration(conv_state%numSpecies), STAT=STATUS)
         VERIFY_(STATUS)
         conv_state%isFixedConcentration(:) = .FALSE.

         DO ic = 1, conv_state%numSpecies
            ! Get field and name from tracer bundle
            !--------------------------------------
            call ESMF_FieldBundleGet(ConvTR, ic, FIELD, RC=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
            VERIFY_(STATUS)

            !Identify fixed species such as O2, N2, ad.
            if (TRIM(NAME) == 'ACET' .OR. TRIM(NAME) == 'N2'   .OR. &
                TRIM(NAME) == 'O2'   .OR. TRIM(NAME) == 'NUMDENS')  THEN
               conv_state%isFixedConcentration(ic) = .TRUE.
            end if

         END DO
      END IF

      ! Get the tracers from the ESMF Bundle
      !-------------------------------------
      ALLOCATE(concentration(conv_state%numSpecies), STAT=STATUS)
      VERIFY_(STATUS)

      DO ic = 1, conv_state%numSpecies
         ! Get field and name from tracer bundle
         !--------------------------------------
         call ESMF_FieldBundleGet(ConvTR, ic, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

         ! Get pointer to the quantity
         !----------------------------
         call ESMFL_BundleGetPointerToData(ConvTR, NAME, S, RC=STATUS)
         VERIFY_(STATUS)

         ! The quantity must exist; others are optional.
         !----------------------------------------------
         ASSERT_(associated(S ))

         ALLOCATE(concentration(ic)%pArray3d(i1:i2,j1:j2,km), STAT=STATUS)
         VERIFY_(STATUS)

         concentration(ic)%pArray3d(:,:,km:1:-1) = S(:,:,:)
      END DO

      ! Satisfy the imports
      !--------------------
      CALL MAPL_GetPointer(impConv,          T,        'T', __RC__)
      CALL MAPL_GetPointer(impConv,        zpbl,    'ZPBL', __RC__)
      CALL MAPL_GetPointer(impConv,         ple,     'PLE', __RC__)
      CALL MAPL_GetPointer(impConv,   totalMass,    'MASS', __RC__)
      CALL MAPL_GetPointer(impConv,         zle,     'ZLE', __RC__)
      CALL MAPL_GetPointer(impConv,     CNV_MFD, 'CNV_MFD', __RC__)
      CALL MAPL_GetPointer(impConv,     CNV_MFC, 'CNV_MFC', __RC__)

      allocate(pbl(i1:i2,j1:j2),                 STAT=STATUS); VERIFY_(STATUS)
      allocate(press3c(i1:i2,j1:j2,k1:k2),       STAT=STATUS); VERIFY_(STATUS)
      allocate(press3e(i1:i2,j1:j2,k1-1:k2),     STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(pl(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(mass(i1:i2,j1:j2,k1:k2),          STAT=STATUS); VERIFY_(STATUS)
      allocate(gridBoxHeight(i1:i2,j1:j2,k1:k2), STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(cmf(i1:i2,j1:j2,k1:k2),           STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(dtrain(i1:i2,j1:j2,k1:k2),        STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(eu(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(ed(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(md(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(kel(i1:i2,j1:j2,k1:k2),           STAT=STATUS); VERIFY_(STATUS)

      pbl     (:,:)       = zpbl(:,:)
      kel     (:,:,k1:k2) = T        (:,:,km:1:-1)
      cmf     (:,:,k1:k2) = CNV_MFC  (:,:,km-1:0:-1)
      dtrain  (:,:,k1:k2) = CNV_MFD  (:,:,km:1:-1)
      mass    (:,:,k1:k2) = totalMass(:,:,km:1:-1)
      pl      (:,:,k1:k2) = (ple(:,:,0:km-1)+ple(:,:,1:km))*0.50

      press3c(i1:i2,j1:j2,:)       = pl (:,:,km:1:-1)*Pa2hPa
      press3e(i1:i2,j1:j2,k1-1:k2) = ple(:,:,km:0:-1)*Pa2hPa

      ! This formulation was suggested by Steve Steentod (25Nov2014)
      DO ik=1,km
         kR = km-ik+1
         eu(:,:,kR) = CNV_MFC(:,:,ik-1) - CNV_MFC(:,:,ik) + dtrain(:,:,kR)
      END DO

      ed = 0.0
      md = 0.0

      DO ik=1,km
         kR = km-ik+1
         gridBoxHeight(:,:,kR) = zle(:,:,ik-1)-zle(:,:,ik)         ! m
      END DO

      tdt = cdt

      call doConvectiveTransport (conv_state%det_ent, conv_state%do_downdraft, pbl, cmf, &
                      dtrain, eu, ed, md, gridBoxHeight, mass, kel, press3e,       &
                      concentration, conv_state%isFixedConcentration, conv_state%cellArea, tdt,  &
                      i1, i2, j1, j2, k1, k2, ilong, ivert, conv_state%numSpecies)

      deallocate(pbl )
      deallocate(dtrain, cmf, kel, eu, ed, md)
      deallocate(press3e, press3c, mass, gridBoxHeight, pl)

      ! Pass back the tracers to the ESMF Bundle
      !------------------------------------------
      DO ic = 1, conv_state%numSpecies
         call ESMF_FieldBundleGet(ConvTR, ic, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

         call ESMFL_BundleGetPointerToData(ConvTR, NAME, S, RC=STATUS)
         VERIFY_(STATUS)

         ! Do not forget to re-order the vertical levels
         S(:,:,:) = concentration(ic)%pArray3d(:,:,km:1:-1)
      END DO

      CALL CleanArrayPointer(concentration, STATUS)
      VERIFY_(STATUS)

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

   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"Finalize"

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   call extract_ ( gc, clock,  nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Finalize MAPL Generic.  Atanas says, "Do not deallocate foreign objects."
!  -------------------------------------------------------------------------
   call MAPL_GenericFinalize ( gc, impConv, expConv, clock,  RC=STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Finalize_
!EOC
!-------------------------------------------------------------------------
    SUBROUTINE extract_ ( gc, clock, nymd, nhms, cdt, rc )

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_Clock), intent(in)       :: clock
    integer, intent(out)               :: nymd, nhms
    real, intent(out)                  :: cdt
    integer, intent(out)               :: rc

    type(T_Convection_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'extract_'

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get time step
!   -------------
    call ESMF_ConfigGetAttribute ( CF, cdt, LABEL="RUN_DT:", RC=STATUS )
    VERIFY_(STATUS)

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
