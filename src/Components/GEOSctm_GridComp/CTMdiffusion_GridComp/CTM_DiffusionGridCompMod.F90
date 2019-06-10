#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: CTM_DiffusionGridCompMod - The Diffusion Grid Component
!
! !INTERFACE:
!
   MODULE CTM_DiffusionGridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL_Mod
   USE GmiDiffusionMethod_mod                   ! ESMF parent component
   USE m_chars, ONLY : uppercase

   IMPLICIT NONE
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices
!
! !DESCRIPTION: 
!
!  {\tt CTM\_DiffusionGridCompMod} is a ESMF gridded component for the 
!  Diffusion package.
!
! !REVISION HISTORY:
!
!  07Dec2012  Kouatchou  Created the Diffusion stub.
!
!EOP
!-------------------------------------------------------------------------
  TYPE Diffusion_State
     PRIVATE
     TYPE(Diffusion_GridComp),  POINTER :: gcDIFF    => null()
  END TYPE Diffusion_State

  TYPE Diffusion_WRAP
     TYPE (Diffusion_State), pointer :: PTR => null()
  END TYPE Diffusion_WRAP

!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for Diffusion Grid Component
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

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Local derived type aliases
!   --------------------------
    type (ESMF_Config)              :: CF
    type (Diffusion_State), pointer   :: state   ! internal, that is
    type (Diffusion_wrap)             :: wrap

    integer                         :: m, n, i_XX, j_XX
    CHARACTER(LEN=ESMF_MAXSTR)      :: FRIENDLIES, providerName

    LOGICAL :: searchForImports
    CHARACTER(LEN=2) :: leadChars
    CHARACTER(LEN=ESMF_MAXSTR) :: name

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


    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *, TRIM(Iam)//': ACTIVE'
    END IF


!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------

    CALL MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, &
                                      RC=STATUS)
    VERIFY_(STATUS)

    CALL MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_RUN,  Run_,        &
                                      RC=STATUS)
    VERIFY_(STATUS)

    CALL MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_FINALIZE,  Finalize_,  &
                                      RC=STATUS)
    VERIFY_(STATUS)

!   Store internal state in GC
!   --------------------------
    CALL ESMF_UserCompSetInternalState ( GC, 'Diffusion_state', wrap, STATUS )
    VERIFY_(STATUS)

! ========================== IMPORT STATE =========================

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'PLE',  &
        LONG_NAME          = 'air_pressure',  &
        UNITS              = 'Pa', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationEdge,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'MASS',  &
        LONG_NAME          = 'total_mass',  &
        UNITS              = 'kg', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'KH',  &
        LONG_NAME          = 'total_scalar_diffusivity',  &
        UNITS              = 'm+2 s-1', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationEdge,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                               &
        SHORT_NAME         = 'TROPP',                                          &
        LONG_NAME          = 'tropopause_pressure_based_on_blended_estimate',  &
        UNITS              = 'Pa',                                             &
        DIMS               = MAPL_DimsHorzOnly,                                &
        VLOCATION          = MAPL_VLocationNone,                               &
                                                                    RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'ZLE',  &
        LONG_NAME          = 'geopotential_height',  &
        UNITS              = 'm', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationEdge,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'ZPBL',  &
        LONG_NAME          = 'planetary_boundary_layer_height',  &
        UNITS              = 'm', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationNone,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC, &
        SHORT_NAME         = 'AREA',  &
        LONG_NAME          = 'agrid_cell_area',  &
        UNITS              = 'm^2', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationNone,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'DiffTR',                            &
        LONG_NAME          = 'diffused_quantities',               &
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART            = MAPL_RestartOptional,                &

                                                       RC=STATUS  )
     VERIFY_(STATUS)


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


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize Diffusion
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( gc, impDiff, expDiff, clock, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impDiff     ! Import State
   type(ESMF_State), intent(inout) :: expDiff     ! Export State
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
!
!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Diffusion_GridComp), pointer     :: gcDIFF      ! Grid Component
   integer                         :: nymd, nhms  ! time of day
   real                            :: cdt         ! chemistry timestep (secs)

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km                  ! dist grid indices
   integer                         :: dims(3), k, l, n

   type(MAPL_MetaComp), pointer    :: ggState	   ! GEOS Generic State
   type (ESMF_State)               :: internal
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)

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
   call MAPL_GenericInitialize ( gc, impDiff, expDiff, clock,  RC=STATUS )
   VERIFY_(STATUS)

   call MAPL_TimerOn(ggSTATE,"TOTAL")
   call MAPL_TimerOn(ggSTATE,"INITIALIZE")

!  Get parameters from gc and clock
!  --------------------------------
   call extract_ ( gc, clock, gcDIFF, nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

   call ESMF_GridCompGet ( GC, GRID=grid, rc=STATUS)
   VERIFY_(STATUS)

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, RC=STATUS)
   VERIFY_(STATUS)

   im = dims(1)
   jm = dims(2)

   call ESMF_GridGet(GRID, localDE=0,staggerloc=ESMF_STAGGERLOC_CENTER, &
   		     computationalCount=DIMS, RC=STATUS)
   VERIFY_(STATUS)

! Local sizes of three dimensions
!--------------------------------
   i2 = dims(1)
   j2 = dims(2)
   km = dims(3)

! Broadcast necessary information to individual GCs
! -------------------------------------------------
   gcDIFF%i1 = i1
   gcDIFF%i2 = i2
   gcDIFF%im = im
   gcDIFF%j1 = j1
   gcDIFF%j2 = j2
   gcDIFF%jm = jm
   gcDIFF%km = km

!  Call initialize
!  ---------------
   call initializeDiffusion ( gcDIFF, impDiff, expDiff, nymd, nhms, &
                              grid, cdt, STATUS )
   VERIFY_(STATUS)

   call MAPL_TimerOff(ggSTATE,"INITIALIZE")
   call MAPL_TimerOff(ggSTATE,"TOTAL")

   RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Initialize_
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs Diffusion
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( gc, impDiff, expDiff, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impDiff     ! Import State
   type(ESMF_State), intent(inout) :: expDiff     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!  07Dec2012 Kouatchou  First crack.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Run'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Diffusion_GridComp), pointer     :: gcDIFF       ! Grid Component
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   integer                         :: i, im, iOX, iT2M, jm, k, km, m, n

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid        
   type(ESMF_Time)                 :: TIME

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, GRID=grid, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"::Run"

!  Get ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, gcDIFF, nymd, nhms, & 
                   cdt, rc=status )
   VERIFY_(STATUS)

!  Run
!  ---
   CALL runDiffusion ( gcDIFF, impDiff, expDiff, nymd, nhms, &
                          cdt, STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize CTM_DiffusionGridCompMod (ESMF)
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( gc, impDiff, expDiff, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impDiff     ! Import State
   type(ESMF_State), intent(inout) :: expDiff     ! Export State
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

   type(Diffusion_GridComp), pointer     :: gcDIFF      ! Grid Component
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)

    type(Diffusion_state), pointer   :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = TRIM(COMP_NAME)//"Finalize"

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   call extract_ ( gc, clock, gcDIFF, nymd, nhms, cdt, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Call ESMF version
!  -----------------
   call FinalizeDiffusion ( gcDIFF )

!  Finalize MAPL Generic.  Atanas says, "Do not deallocate foreign objects."
!  -------------------------------------------------------------------------
   call MAPL_GenericFinalize ( gc, impDiff, expDiff, clock,  RC=STATUS )
   VERIFY_(STATUS)

!  Destroy Legacy state
!  --------------------
   deallocate ( state%gcDIFF, stat = STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Finalize_
!EOC
!-------------------------------------------------------------------------
    SUBROUTINE extract_ ( gc, clock, gcDIFF, nymd, nhms, cdt, &
                          rc, state )

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_Clock), intent(in)       :: clock
    type(Diffusion_GridComp), pointer  :: gcDIFF
    integer, intent(out)               :: nymd, nhms
    real, intent(out)                  :: cdt
    integer, intent(out)               :: rc
    type(Diffusion_state), pointer, optional   :: state


    type(Diffusion_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(Diffusion_Wrap)  :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'extract_'

    rc = 0

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'Diffusion_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

!   This is likely to be allocated during initialize only
!   -----------------------------------------------------
    if ( .not. associated(myState%gcDIFF) ) then
         allocate ( myState%gcDIFF, stat=STATUS )
         VERIFY_(STATUS)
    end if

    gcDIFF  => myState%gcDIFF

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

 END MODULE CTM_DiffusionGridCompMod
