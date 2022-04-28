#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP

! !MODULE:  GEOS\_ctmGridCompMod -- A Module to combine Chemistry, 
!                                   advCore (Transport), Convection and 
!                                   Diffusion Gridded Components
!
! !INTERFACE:
!
      module GEOS_ctmGridCompMod
!
! !USES:
      use ESMF
      use MAPL
      use Bundle_IncrementMod

      use GEOS_ctmHistGridComp,      only : HctmSetServices  => SetServices
      use GEOS_ctmEnvGridComp,       only : EctmSetServices  => SetServices
      use CTM_ConvectionGridCompMod, only : ConvSetServices  => SetServices
      use CTM_DiffusionGridCompMod,  only : DiffSetServices  => SetServices
      use GEOS_ChemGridCompMod,      only : ChemSetServices  => SetServices
      use AdvCore_GridCompMod,       only : AdvCSetServices  => SetServices
      use CTM_pTracersGridCompMod,   only : pTraSetServices  => SetServices
      use m_set_eta, only: set_eta
      use, intrinsic :: iso_fortran_env, only: REAL64
      use Chem_GroupMod
!
      implicit none
      private
!
! !PUBLIC MEMBER FUNCTIONS:

      public SetServices
!
! !DESCRIPTION: 
!   This gridded component (GC) combines Transport (AdvCore), 
!   Convection, GEOSchem, and Diffusion GCs into a new 
!   composite GEOSctm GC.
!   The friendly tracers are variables from GEOSchem.  
! 
! \paragraph{Runnng the Code:}
!
!   The code can be run in two main configurations:
!   \begin{enumerate}
!   \item \textbf{Passive Tracer Run:} This experiment is done to verify how well 
!         AdvCore transports the tracers. We want to find out if the advection 
!         module conserves the mass of each tracer over time. To carry out this 
!         experiment, you need to initialize the tracers. This can be done in
!         two possibe ways:
!         \begin{enumerate} 
!         \item Idealized tracers which concentrations are computed in the
!               initialize method of the pTracer component. To select this
!               option, set \emph{do\_AdvColdStart: T} in the resource file
!               \texttt{pTracer\_GridComp.rc}.
!         \item User provided restart file (on the cubed sphere grid) with
!               concentration of each tracer.
!         \end{enumerate} 
!   \item \textbf{ Chemistry Run:} Here we exercise any Chemistry configuration
!         (GOCART, GMI, GEOS CHEM).
!          We have the option to determine if we want to do Convection and/or
!          Diffusion by setting the variables \emph{do\_ctmConvection} and
!          \emph{do\_ctmDiffusion} in the resource file \texttt{CTM\_GridComp.rc}.
!   \end{enumerate}
!
! \paragraph{External Data Files:}
!
!   To run this code, you need to have data files that have at least
!   the following variables:
!
!  \begin{description}
!  \item[PLE]:  edge pressure 
!  \item[U]:    eastward wind
!  \item[V]:    northward wind
!  \item[ZLE]:  geopotential height
!  \item[Q]:    specific humidity
!  \item[T]:    temperatute
!  \end{description}
!
!  The above variables are required to for passive tracer experiments.
!  More variables should be provided for other experiments.
!  Regarless of the type of run you choose to carry out, you will need
!  to edit a resource file (named \texttt{MAPL_ExtData_state.rc}) that
!  lists (among other information) each variable and the external file
!  which contains the variable.
!
! 
!EOP
!------------------------------------------------------------------------------
      !Derived types for the internal state
      type T_CTM_STATE
         private
         integer                    :: CHEM = -1
         integer                    :: CONV = -1
         integer                    :: DIFF = -1
         integer                    :: ADV3 = -1
         integer                    :: ECTM = -1
         integer                    :: HCTM = -1
         integer                    :: PTRA = -1
         logical                    :: do_ctmConvection = .FALSE. ! do Convection?
         logical                    :: do_ctmDiffusion  = .FALSE. ! do Diffusion?
         logical                    :: do_ctmAdvection  = .TRUE.  ! do Advection?
         logical                    :: output_forcingData  = .FALSE. ! Export variables to HISTORY?
         logical                    :: enable_pTracers  = .FALSE. ! do idealized Passive Tracers?
         character(len=ESMF_MAXSTR) :: metType                    ! MERRA2 or MERRA1 or FPIT or FP
      end type T_CTM_STATE

      type CTM_WRAP
         type (T_CTM_STATE), pointer :: PTR
      end type CTM_WRAP

!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
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
    integer,             intent(  OUT) :: RC  ! return code
!
! !DESCRIPTION:  
!   The SetServices for the GEOSctm GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (CHEM, Diffusion, Convection, advCore) and runs their
!   respective SetServices.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer                       :: STATUS
      integer                       :: I
      type (ESMF_Config)            :: CF
      type (ESMF_Config)            :: configFile
      character(len=ESMF_MAXSTR)    :: COMP_NAME
      CHARACTER(LEN=ESMF_MAXSTR)    :: rcfilen = 'CTM_GridComp.rc'
      character(len=ESMF_MAXSTR)    :: IAm = 'SetServices'
      type (T_CTM_STATE), pointer   :: state
      type (CTM_WRAP)               :: wrap

      ! Get my name and set-up traceback handle
      ! ---------------------------------------

      Iam = 'SetServices'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // "::" // Iam

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
      call ESMF_UserCompSetInternalState ( GC, 'CTM_GridComp', wrap, STATUS )
      VERIFY_(STATUS)

      ! Choose children to birth and which children not to conceive
      ! -----------------------------------------------------------
      configFile = ESMF_ConfigCreate( __RC__ )

      call ESMF_ConfigLoadFile(configFile, TRIM(rcfilen),  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%output_forcingData,     &
                             Default  = .FALSE.,                       &
                             Label    = "output_forcingData:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%enable_pTracers,       &
                                     Default  = .FALSE.,                    &
                                     Label    = "ENABLE_pTracers:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%do_ctmConvection,      &
                                     Default  = .FALSE.,                    &
                                     Label    = "do_ctmConvection:", __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%do_ctmAdvection,       &
                                     Default  = .TRUE.,                     &
                                     Label    = "do_ctmAdvection:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, state%do_ctmDiffusion,       &
                                     Default  = .FALSE.,                    &
                                     Label    = "do_ctmDiffusion:",  __RC__ )

      ! Type of meteological fields (MERRA2 or MERRA1 or FPIT or FP)
      call ESMF_ConfigGetAttribute(configFile, state%metType,               &
                                     Default  = 'MERRA2',                   &
                                     Label    = "metType:",          __RC__ )

      IF ((TRIM(state%metType) == "F515_516") .OR. &
          (TRIM(state%metType) == "F5131"))            state%metType = "FP"

      ! Turn Convection on for any Chemistry configuration
      IF (.NOT. state%enable_pTracers) THEN
         state%do_ctmConvection = .TRUE.
      ENDIF

      IF ( MAPL_am_I_root() ) THEN
         PRINT*
         PRINT*, "---------------------------------------------------"
         PRINT*, "-----             GEOS CTM Settings           -----"
         PRINT*, "---------------------------------------------------"
         PRINT*,'   Doing Passive Tracer?: ', state%enable_pTracers
         PRINT*,'    Output forcing data?: ', state%output_forcingData
         PRINT*,'               Advection: ', state%do_ctmAdvection
         PRINT*,'              Convection: ', state%do_ctmConvection
         PRINT*,'               Diffusion: ', state%do_ctmDiffusion
         PRINT*,'     Meteological Fields: ', TRIM(state%metType)
         PRINT*, "---------------------------------------------------"
         PRINT*
      END IF

      ! -----------------------------------------------------------------
      ! Create children`s gridded components and invoke their SetServices
      ! -----------------------------------------------------------------

      IF (state%enable_pTracers) THEN
         ! Doing passive tracer experiment
         !--------------------------------
         state%ADV3 = MAPL_AddChild(GC, NAME='DYNAMICS',   SS=AdvCSetServices, __RC__)
         state%ECTM = MAPL_AddChild(GC, NAME='CTMenv',     SS=EctmSetServices, __RC__)
         state%PTRA = MAPL_AddChild(GC, NAME='PTRACERS',   SS=pTraSetServices, __RC__)

         ! Do you want to provide forcing data to HISTOTY?
         if (state%output_forcingData) then
            state%HCTM = MAPL_AddChild(GC, NAME='CTMhist',    SS=HctmSetServices, __RC__)
         endif

         ! Are you doing Convection?
         if (state%do_ctmConvection) then
            state%CONV = MAPL_AddChild(GC, NAME='CONVECTION', SS=ConvSetServices, __RC__)
         end if

         ! Are you doing Diffusion?
         if (state%do_ctmDiffusion) then
            state%DIFF = MAPL_AddChild(GC, NAME='DIFFUSION',  SS=DiffSetServices, __RC__ )
         end if
      ELSE
         state%ADV3 = MAPL_AddChild(GC, NAME='DYNAMICS',   SS=AdvCSetServices, __RC__)
         state%CHEM = MAPL_AddChild(GC, NAME='CHEMISTRY',  SS=ChemSetServices, __RC__)
         state%ECTM = MAPL_AddChild(GC, NAME='CTMenv',     SS=EctmSetServices, __RC__)

         ! Do you want to provide forcing data to HISTOTY?
         if (state%output_forcingData) then
            state%HCTM = MAPL_AddChild(GC, NAME='CTMhist',    SS=HctmSetServices, __RC__)
         endif

         ! Are you doing Convection?
         if (state%do_ctmConvection) then
            state%CONV = MAPL_AddChild(GC, NAME='CONVECTION', SS=ConvSetServices, __RC__)
         end if

         ! Are you doing Diffusion?
         if (state%do_ctmDiffusion) then
            state%DIFF = MAPL_AddChild(GC, NAME='DIFFUSION',  SS=DiffSetServices, __RC__)
         end if
      END IF

      call MAPL_TimerAdd(GC, name="INITIALIZE", __RC__)
      call MAPL_TimerAdd(GC, name="RUN"       , __RC__)

      ! -------------------------------
      ! Connectivities between Children
      ! -------------------------------
      CALL MAPL_AddConnectivity ( GC, &
               SHORT_NAME  = (/'AREA'/), &
               DST_ID = state%ECTM, SRC_ID = state%ADV3, __RC__  )

      CALL MAPL_AddConnectivity ( GC, &
               SRC_NAME  = (/ 'CXr8  ', 'CYr8  ', 'MFXr8 ', 'MFYr8 ', 'PLE0r8', 'PLE1r8' /), &
               DST_NAME  = (/ 'CX    ', 'CY    ', 'MFX   ', 'MFY   ', 'PLE0  ', 'PLE1  ' /), &
               DST_ID = state%ADV3, SRC_ID = state%ECTM, __RC__  )
      CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'TRADV'/),          &
               CHILD = state%ADV3,    __RC__  )

      IF (state%output_forcingData) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'PLE'/), &
                 DST_ID = state%HCTM, SRC_ID = state%ECTM, __RC__  )
      ENDIF

      ! Doing Convection
      IF (state%do_ctmConvection) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/ 'CNV_MFC', 'CNV_MFD' /), &
                 DST_ID = state%CONV, SRC_ID = state%ECTM, __RC__  )

         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = state%CONV, SRC_ID = state%ADV3, __RC__  )

         IF ( (TRIM(state%metType) == 'MERRA2') .OR. &
              (TRIM(state%metType) == 'FPIT')   .OR. &
              (TRIM(state%metType) == 'FP') ) THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'LWI '/), &
                    DST_ID = state%CONV, SRC_ID = state%ECTM, __RC__  )
         ELSEIF ( TRIM(state%metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'ZLE ', 'LWI '/), &
                    DST_ID = state%CONV, SRC_ID = state%ECTM, __RC__  )
         END IF

         CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'ConvTR'/),          &
               CHILD = state%CONV,    __RC__  )
      END IF

      ! Doing Diffusion
      IF (state%do_ctmDiffusion) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = state%DIFF, SRC_ID = state%ADV3, __RC__  )

         IF ( (TRIM(state%metType) == 'MERRA2') .OR. &
              (TRIM(state%metType) == 'FPIT')   .OR. &
              (TRIM(state%metType) == 'FP') ) THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS'/), &
                    DST_ID = state%DIFF, SRC_ID = state%ECTM, __RC__  )
         ELSEIF ( TRIM(state%metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'ZLE '/), &
                    DST_ID = state%DIFF, SRC_ID = state%ECTM, __RC__  )
         END IF

         CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'DiffTR'/),          &
               CHILD = state%DIFF,    __RC__  )
      END IF

      IF (state%enable_pTracers) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = state%PTRA, SRC_ID = state%ADV3, __RC__  )

         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'PLE'/), &
                 DST_ID = state%PTRA, SRC_ID = state%ECTM, __RC__  )
      ELSE
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = state%CHEM, SRC_ID = state%ADV3, __RC__  )

         IF (state%do_ctmConvection) THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/ 'CNV_MFC', 'CNV_MFD' /), &
                    DST_ID = state%CHEM, SRC_ID = state%ECTM, __RC__  )
         ENDIF

         ! This is done for MERRA1, MERRA2, FPIT, FP
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/ 'LFR   ', 'QCTOT ', 'TH    ', 'PLE   ', &
                                  'LWI   ', 'CNV_QC', 'U     ', 'V     ', &
                                  'BYNCY ', 'ITY   ', 'QICN  ', 'QLCN  ' /), &
                 DST_ID = state%CHEM, SRC_ID = state%ECTM, __RC__  )

         ! Additional Connectivity if using MERRA1
         IF ( TRIM(state%metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/ 'ZLE   ', 'RH2   ' /), &
                    DST_ID = state%CHEM, SRC_ID = state%ECTM, __RC__  )
         END IF
      END IF

      ! EXPORT States for Increment Bundles:
      !---------------
      call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME  = 'TRADVI',                                 &
           LONG_NAME   = 'advected_quantities_tendencies',         &
           units       = 'UNITS s-1',                              &
           DATATYPE    = MAPL_BundleItem,         __RC__ )

      call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME  = 'MTRI',                                   &
           LONG_NAME   = 'moist_quantities_tendencies',            &
           units       = 'UNITS s-1',                              &
           DATATYPE    = MAPL_BundleItem,          __RC__ )

      call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME  = 'TRI',                                    &
           LONG_NAME   = 'turbulence_quantities_tendencies',       &
           units       = 'UNITS s-1',                              &
           DATATYPE    = MAPL_BundleItem,        __RC__ )

      ! Create grid for this GC
      !------------------------
      !call MAPL_GridCreate  (GC, __RC__ )

      call MAPL_GenericSetServices ( GC, __RC__ )

      RETURN_(ESMF_SUCCESS)
  
      end subroutine SetServices
!EOC
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: Initialize -- Initialize method for the GEOS CTM Gridded Component
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
! !INPUT/OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
!  The Initialize method of the GEOS CTM Gridded Component.
!  It acts as a driver for the initializtion of the four children: 
!  GEOSchem, advCore, Diffusion and Convection. 
!  It also sets up the frieldly connections between the children.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer                             :: STATUS
      type (MAPL_MetaComp),       pointer :: STATE
      type (ESMF_GridComp),       pointer :: GCS(:)
      type (ESMF_State),          pointer :: GIM(:)
      type (ESMF_State),          pointer :: GEX(:)
      type (ESMF_FieldBundle)             :: BUNDLE, iBUNDLE
      type (ESMF_Field)                   :: FIELD
      type (ESMF_State)                   :: DUMMY
      type (ESMF_Grid)                    :: GRID

      integer                             :: NUM_TRACERS
      integer                             :: I
      integer                             :: NA
      character(len=ESMF_MAXSTR), pointer :: NAMES(:)
      character(len=ESMF_MAXSTR)          :: myNAME
      character(len=ESMF_MAXSTR)          :: varName
      character(len=ESMF_MAXSTR)          ::  iNAME
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      character(len=ESMF_MAXSTR)          :: IAm = "Initialize"
      real(REAL64), allocatable           :: ak(:),bk(:)
      real(REAL64)                        :: ptop, pint
      integer                             :: counts(3),lm,ls, ib
      type (T_CTM_STATE), pointer         :: CTM_STATE
      type (CTM_WRAP)                     :: WRAP


      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------

      !call ESMF_GridCompGet ( GC, name=COMP_NAME, GRID=GRID, RC=STATUS )
      call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
      Iam = trim(COMP_NAME) // "::" // TRIM(Iam)

      ! Create grid for this GC
      !------------------------
      call MAPL_GridCreate  (GC, __RC__ )

      ! Add the AK and BK information
      ! -----------------------------
      call ESMF_GridCompGet(GC,grid=grid,__RC__)
      call MAPL_GridGet(grid, globalCellCountPerDim=counts,__RC__)
      lm = counts(3)
      allocate(ak(lm+1),stat=status)
      VERIFY_(STATUS)
      allocate(bk(lm+1),stat=status)
      VERIFY_(STATUS)
      call set_eta(lm,ls,ptop,pint,ak,bk)
      call ESMF_AttributeSet(grid, name='GridAK', itemCount=LM+1, &
                             valuelist=ak, __RC__)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(grid, name='GridBK', itemCount=LM+1, &
                             valuelist=bk, __RC__)
      VERIFY_(STATUS)
      deallocate(ak,bk)

      ! Get my MAPL_Generic state
      !--------------------------

      call MAPL_GetObjectFromGC ( GC, STATE, __RC__ )

      ! Call Initialize for every Child
      !--------------------------------
      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  __RC__ )

      call MAPL_TimerOn(STATE,"TOTAL")
      call MAPL_TimerOn(STATE,"INITIALIZE")

      ! Get my private state from the component
      !----------------------------------------
      call ESMF_UserCompGetInternalState(GC, 'CTM_GridComp', WRAP, STATUS )
      VERIFY_(STATUS)

      CTM_STATE => WRAP%PTR

#ifdef PRINT_STATES
      call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
      if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, __RC__ )
      call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
      if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, __RC__ )
#endif


      ! Get children and their im/ex states from my generic state.
      !----------------------------------------------------------

      call MAPL_Get ( STATE, GCS=GCS, GIM=GIM, GEX=GEX, __RC__ )

      ! Extract the friendly tracers
      !-----------------------------
      IF (CTM_STATE%enable_pTracers) THEN
         if (CTM_STATE%do_ctmDiffusion) then
            !------------------
            ! Diffusion Tracers
            !------------------
            call ESMF_StateGet   (GIM(CTM_STATE%DIFF),  'DiffTR' , BUNDLE, __RC__  )

            call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%PTRA), "TURBULENCE", BUNDLE, __RC__  )

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Diffusion Tracer Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

            ! Fill the diffusion increments bundle
            !---------------------------------
            call Initialize_IncBundle_init(GC, GIM(CTM_STATE%DIFF), EXPORT, TRIincCTM, __RC__)

            ! Count tracers
            !--------------

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )

            ! Get the names of all tracers to fill other turbulence bundles.
            !---------------------------------------------------------------

            allocate(NAMES(NUM_TRACERS),STAT=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  __RC__ )
         end if

         if (CTM_STATE%do_ctmConvection) then
            !-------------------
            ! Convection Tracers
            !-------------------
            call ESMF_StateGet       (GIM(CTM_STATE%CONV), 'ConvTR', BUNDLE, __RC__  )

            call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%PTRA), "MOIST", BUNDLE, __RC__  )

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Convective Transport Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

            ! Fill the moist increments bundle
            !---------------------------------
            call Initialize_IncBundle_init(GC, GIM(CTM_STATE%CONV), EXPORT, MTRIincCTM, __RC__)

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )
         end if

         !----------------
         ! AdvCore Tracers
         !----------------
         call ESMF_StateGet       (GIM(CTM_STATE%ADV3), 'TRADV', BUNDLE, __RC__  )
  
         call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%PTRA), "DYNAMICS", BUNDLE, __RC__  )

#ifdef PRINT_STATES
         call WRITE_PARALLEL ( trim(Iam)//": AdvCore Tracer Bundle" )
         if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

         call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )

         ! Initialize the advection increments bundle (TRADVI)
         ! with tracer increment names 
         !--------------------------------
         call Initialize_IncBundle_init(GC, GIM(CTM_STATE%ADV3), EXPORT, DYNinc, __RC__)


         ! Get the names of all tracers to fill other turbulence bundles.
         !---------------------------------------------------------------

         allocate(NAMES(NUM_TRACERS),STAT=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  __RC__ )

      ELSE 
         if (CTM_STATE%do_ctmDiffusion) then
            !------------------
            ! Diffusion Tracers
            !------------------
            call ESMF_StateGet   (GIM(CTM_STATE%DIFF),  'DiffTR' , BUNDLE, __RC__  )
  
            call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%CHEM), "TURBULENCE", BUNDLE, __RC__  )

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Diffusion Tracer Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

            ! Fill the diffusion increments bundle
            !---------------------------------
            call Initialize_IncBundle_init(GC, GIM(CTM_STATE%DIFF), EXPORT, TRIincCTM, __RC__)

            ! Count tracers
            !--------------

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )

            ! Get the names of all tracers to fill other turbulence bundles.
            !---------------------------------------------------------------

            allocate(NAMES(NUM_TRACERS),STAT=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  __RC__ )
         end if

         if (CTM_STATE%do_ctmConvection) then
            !-------------------
            ! Convection Tracers
            !-------------------
            call ESMF_StateGet       (GIM(CTM_STATE%CONV), 'ConvTR', BUNDLE, __RC__  )
  
            call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%CHEM), "MOIST", BUNDLE, __RC__  )

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Convective Transport Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

            ! Fill the moist increments bundle
            !---------------------------------
            call Initialize_IncBundle_init(GC, GIM(CTM_STATE%CONV), EXPORT, MTRIincCTM, __RC__)

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )

            IF (NUM_TRACERS .EQ. 0) THEN
               IF ( MAPL_am_I_root() ) THEN
                  PRINT*, '======================================='
                  PRINT*, '-----> No tracer was friendly to MOIST'
                  PRINT*, '-----> Convection will not be performed'
                  PRINT*, '======================================='
               END IF
               CTM_STATE%do_ctmConvection = .FALSE.
            ELSE
               IF ( MAPL_am_I_root() ) THEN
                  PRINT*, '======================================='
                  PRINT*, '------ List of Convected Tracers ------'
                  PRINT*, '======================================='
                  DO ib = 1, NUM_TRACERS
                     call ESMF_FieldBundleGet(BUNDLE, ib, FIELD,  __RC__)
                     call ESMF_FieldGet(FIELD, name=varName,  __RC__)
                     PRINT '(i4,a5,a20)', ib, '-->', TRIM(varName)
                  ENDDO
                  PRINT*, '======================================='
               END IF
            END IF
         end if

         !----------------
         ! AdvCore Tracers
         !----------------
         call ESMF_StateGet       (GIM(CTM_STATE%ADV3), 'TRADV', BUNDLE, __RC__  )

         call MAPL_GridCompGetFriendlies(GCS(CTM_STATE%CHEM), "DYNAMICS", BUNDLE, __RC__  )

#ifdef PRINT_STATES
         call WRITE_PARALLEL ( trim(Iam)//": AdvCore Tracer Bundle" )
         if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, __RC__ )
#endif

         ! Initialize the advection increments bundle (TRADVI)
         ! with tracer increment names
         !--------------------------------
         call Initialize_IncBundle_init(GC, GIM(CTM_STATE%ADV3), EXPORT, DYNinc, __RC__)

         call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, __RC__ )
      END IF

      call MAPL_TimerOff(STATE,"INITIALIZE")
      call MAPL_TimerOff(STATE,"TOTAL")

      ! All Done
      !---------
      RETURN_(ESMF_SUCCESS)

      end subroutine Initialize
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run -- Run method for the GEOS CTM Gridded Component
!
! !INTERFACE:

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
!  The run method for the GEOSctm calls the children`s
!   run methods. 
!
!EOP  
!=============================================================================
!BOC
!
! !LOCAL VARIABLES:
      integer                             :: STATUS
      type (MAPL_MetaComp),      pointer  :: STATE
      type (ESMF_GridComp),      pointer  :: GCS(:)
      type (ESMF_State),         pointer  :: GIM(:)
      type (ESMF_State),         pointer  :: GEX(:)
      type (ESMF_State)                   :: INTERNAL
      type (ESMF_Config)                  :: CF
      character(len=ESMF_MAXSTR),pointer  :: GCNames(:)
      integer                             :: I, L, K
      integer                             :: IM, JM, LM
      real                                :: DT
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      character(len=ESMF_MAXSTR)          :: IAm = "Run"

      type (ESMF_FieldBundle)             :: Bundle

      character(len=ESMF_MAXSTR), pointer :: Names(:)
      character(len=ESMF_MAXSTR)          :: STRING

      CHARACTER(LEN=ESMF_MAXSTR)          :: fieldName

      real, pointer, dimension(:,:)       :: AREA   => null()
      real, pointer, dimension(:,:,:)     :: PLE    => null()
      real, pointer, dimension(:,:,:)     :: Q      => null()

      integer                             :: NQ
      type (ESMF_FieldBundle)             :: Bundletest
      type (ESMF_Field)                   :: FIELD
      type (T_CTM_STATE), pointer         :: CTM_STATE
      type (CTM_WRAP)                     :: WRAP

      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------

      call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, __RC__ )
      Iam = trim(COMP_NAME) // "::" // TRIM(Iam)

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(STATE,"TOTAL")
      call MAPL_TimerOn(STATE,"RUN")

      ! Get my private state from the component
      !----------------------------------------
      call ESMF_UserCompGetInternalState(GC, 'CTM_GridComp', WRAP, STATUS )
      VERIFY_(STATUS)

      CTM_STATE => WRAP%PTR

      ! Get the children`s states from the generic state
      !-------------------------------------------------

      call MAPL_Get ( STATE,   &
          GCS=GCS, GIM=GIM, GEX=GEX,       &
          IM = IM, JM = JM, LM = LM,       &
          GCNames = GCNames,               &
          INTERNAL_ESMF_STATE = INTERNAL,  __RC__ )

      call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:" , __RC__ )

      !---------------------
      ! Cinderella Component: to derive variables for other components
      !---------------------
      I = CTM_STATE%ECTM

      call MAPL_TimerOn (STATE,GCNames(I))
      call ESMF_GridCompRun (GCS(I),               &
                             importState = GIM(I), &
                             exportState = GEX(I), &
                             clock       = CLOCK,  &
                             userRC      = STATUS  )
      VERIFY_(STATUS)
      call MAPL_TimerOff(STATE,GCNames(I))

      !-----------------------------------------------
      ! CTM History to provide forcing data to HISTORY
      !-----------------------------------------------
      if (CTM_STATE%output_forcingData) then
         I = CTM_STATE%HCTM

          call MAPL_TimerOn (STATE,GCNames(I))
          call ESMF_GridCompRun (GCS(I),               &
                                 importState = GIM(I), &
                                 exportState = GEX(I), &
                                 clock       = CLOCK,  &
                                 userRC      = STATUS  )
          VERIFY_(STATUS)
          call MAPL_TimerOff(STATE,GCNames(I))
      end if

      !--------
      ! advCore
      !--------
      
      ! Initialize Dynamics increment bundle
      !--------------------------------------------
      call Initialize_IncBundle_run(GIM(CTM_STATE%ADV3), EXPORT, DYNinc, __RC__)

      IF (CTM_STATE%do_ctmAdvection) THEN
         I = CTM_STATE%ADV3

         call Pack_Chem_Groups( GIM(I) )  ! Prepare to transport chemical families

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                              importState = GIM(I), &
                              exportState = GEX(I), &
                              clock       = CLOCK,  &
                              userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))

         ! Compute Dynamics increments and fill bundle
         !--------------------------------------------
         call Compute_IncBundle(GIM(CTM_STATE%ADV3), EXPORT, DYNinc, STATE, __RC__)

         call MAPL_GetPointer( GEX(CTM_STATE%ADV3), AREA, 'AREA', __RC__ )
         call MAPL_GetPointer( GEX(CTM_STATE%ECTM), PLE,  'PLE',  __RC__ )
         call MAPL_GetPointer( GIM(CTM_STATE%ECTM), Q,      'Q',  __RC__ )
         call Unpack_Chem_Groups( GIM(CTM_STATE%ADV3), PLE, AREA, Q )   ! Finish transporting chemical families
      END IF

      !-----------
      ! Convection
      !-----------
      IF (CTM_STATE%do_ctmConvection) THEN
         I = CTM_STATE%CONV

         ! Initialize Moist increment bundle
         !----------------------------------
         call Initialize_IncBundle_run(GIM(I), EXPORT, MTRIincCTM, __RC__)

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                exportState = GEX(I), &
                                   clock       = CLOCK,  &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))

         ! Compute Moist increments and fill bundle
         !--------------------------------------------
         call Compute_IncBundle(GIM(I), EXPORT, MTRIincCTM, STATE, __RC__)
      END IF

      IF (.NOT. CTM_STATE%enable_pTracers) THEN
         !----------
         ! Chemistry: Phase 1
         !----------
         I = CTM_STATE%CHEM   

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                exportState = GEX(I), &
                                clock       = CLOCK,  &
                                phase       = 1,      &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      END IF


      !----------
      ! Diffusion
      !----------
      IF (CTM_STATE%do_ctmDiffusion) THEN
         I = CTM_STATE%DIFF

         ! Initialize Diffusion increment bundle
         !----------------------------------
         call Initialize_IncBundle_run(GIM(I), EXPORT, TRIincCTM, __RC__)

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                   exportState = GEX(I), &
                                clock       = CLOCK,  &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))

         ! Compute Diffusion increments and fill bundle
         !--------------------------------------------
         call Compute_IncBundle(GIM(I), EXPORT, TRIincCTM, STATE, __RC__)
      END IF

      IF (CTM_STATE%enable_pTracers) THEN
         !---------------
         ! Passive Tracer
         !---------------
         I = CTM_STATE%PTRA

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                exportState = GEX(I), &
                                clock       = CLOCK,  &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      ELSE
         !----------
         ! Chemistry: Phase 2
         !----------
         I = CTM_STATE%CHEM   

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                exportState = GEX(I), &
                                clock       = CLOCK,  &
                                phase       = 2,      &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      END IF

      call MAPL_TimerOff(STATE,"RUN")
      call MAPL_TimerOff(STATE,"TOTAL")

      RETURN_(ESMF_SUCCESS)

      end subroutine Run
!
!EOC
!------------------------------------------------------------------------------
      end module GEOS_ctmGridCompMod
