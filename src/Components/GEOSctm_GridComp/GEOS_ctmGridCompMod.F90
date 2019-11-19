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
      use MAPL_Mod

      use GEOS_ctmEnvGridComp,       only : EctmSetServices  => SetServices
      use CTM_ConvectionGridCompMod, only : ConvSetServices  => SetServices
      use CTM_DiffusionGridCompMod,  only : DiffSetServices  => SetServices
      use GEOS_ChemGridCompMod,      only : ChemSetServices  => SetServices
      use AdvCore_GridCompMod,       only : AdvCSetServices  => SetServices
      use CTM_pTracersGridCompMod,   only : pTraSetServices  => SetServices
      use m_set_eta, only: set_eta
      use, intrinsic :: iso_fortran_env, only: REAL64
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
!   The code acn be run in two main configurations:
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
      integer :: CHEM = -1
      integer :: CONV = -1
      integer :: DIFF = -1
      integer :: ADV3 = -1
      integer :: ECTM = -1
      integer :: PTRA = -1

      logical :: enable_pTracers  = .FALSE.
      logical :: do_ctmAdvection  = .FALSE.
      logical :: do_ctmConvection = .FALSE.
      logical :: do_ctmDiffusion  = .FALSE.
      character(len=ESMF_MAXSTR) :: metType ! MERRA2 or MERRA1 or FPIT or FP
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

      ! Get my name and set-up traceback handle
      ! ---------------------------------------

      Iam = 'SetServices'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // "::" // Iam

      ! Register services for this component
      ! ------------------------------------

      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        RC=STATUS )
      VERIFY_(STATUS)

      ! Choose children to birth and which children not to conceive
      ! -----------------------------------------------------------
      configFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(configFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(configFile, enable_pTracers,             &
                                     Default  = .FALSE.,                    &
                                     Label    = "ENABLE_pTracers:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, do_ctmConvection,            &
                                     Default  = .FALSE.,                    &
                                     Label    = "do_ctmConvection:", __RC__ )

      call ESMF_ConfigGetAttribute(configFile, do_ctmAdvection,             &
                                     Default  = .TRUE.,                     &
                                     Label    = "do_ctmAdvection:",  __RC__ )

      call ESMF_ConfigGetAttribute(configFile, do_ctmDiffusion,             &
                                     Default  = .FALSE.,                    &
                                     Label    = "do_ctmDiffusion:",  __RC__ )

      ! Type of meteological fields (MERRA2 or MERRA1 or FPIT or FP)
      call ESMF_ConfigGetAttribute(configFile, metType,                     &
                                     Default  = 'MERRA2',                   &
                                     Label    = "metType:",          __RC__ )

      IF ((TRIM(metType) == "F515_516") .OR. &
          (TRIM(metType) == "F5131"))            metType = "FP"

      IF ( MAPL_am_I_root() ) THEN
         PRINT*
         PRINT*, "---------------------------------------------------"
         PRINT*, "-----             GEOS CTM Settings           -----"
         PRINT*, "---------------------------------------------------"
         PRINT*,'   Doing Passive Tracer?: ', enable_pTracers
         PRINT*,'               Advection: ', do_ctmAdvection
         PRINT*,'              Convection: ', do_ctmConvection
         PRINT*,'               Diffusion: ', do_ctmDiffusion
         PRINT*,'     Meteological Fields: ', TRIM(metType)
         PRINT*, "---------------------------------------------------"
         PRINT*
      END IF

      ! -----------------------------------------------------------------
      ! Create children`s gridded components and invoke their SetServices
      ! -----------------------------------------------------------------

      IF (enable_pTracers) THEN
         ! Doing passive tracer experiment
         !--------------------------------
         ADV3 = MAPL_AddChild(GC, NAME='DYNAMICS',   SS=AdvCSetServices, __RC__)
         ECTM = MAPL_AddChild(GC, NAME='CTMenv',     SS=EctmSetServices, __RC__)
         PTRA = MAPL_AddChild(GC, NAME='PTRACERS',   SS=pTraSetServices, __RC__)

         ! Are you doing Convection?
         if (do_ctmConvection) then
            CONV = MAPL_AddChild(GC, NAME='CONVECTION', SS=ConvSetServices, __RC__)
         end if

         ! Are you doing Diffusion?
         if (do_ctmDiffusion) then
            DIFF = MAPL_AddChild(GC, NAME='DIFFUSION',  SS=DiffSetServices, __RC__)
         end if
      ELSE
         ADV3 = MAPL_AddChild(GC, NAME='DYNAMICS',   SS=AdvCSetServices, __RC__)
         CHEM = MAPL_AddChild(GC, NAME='CHEMISTRY',  SS=ChemSetServices, __RC__)
         ECTM = MAPL_AddChild(GC, NAME='CTMenv',     SS=EctmSetServices, __RC__)

         ! Are you doing Convection?
         if (do_ctmConvection) then
            CONV = MAPL_AddChild(GC, NAME='CONVECTION', SS=ConvSetServices, __RC__)
         end if

         ! Are you doing Diffusion?
         if (do_ctmDiffusion) then
            DIFF = MAPL_AddChild(GC, NAME='DIFFUSION',  SS=DiffSetServices, __RC__)
         end if
      END IF

      call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
      VERIFY_(STATUS)

      ! -------------------------------
      ! Connectivities between Children
      ! -------------------------------
      CALL MAPL_AddConnectivity ( GC, &
               SHORT_NAME  = (/'AREA'/), &
               DST_ID = ECTM, SRC_ID = ADV3, __RC__  )

      CALL MAPL_AddConnectivity ( GC, &
               SRC_NAME  = (/ 'CXr8  ', 'CYr8  ', 'MFXr8 ', 'MFYr8 ', 'PLE0r8', 'PLE1r8' /), &
               DST_NAME  = (/ 'CX    ', 'CY    ', 'MFX   ', 'MFY   ', 'PLE0  ', 'PLE1  ' /), &
               DST_ID = ADV3, SRC_ID = ECTM, __RC__  )
      CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'TRADV'/),          &
               CHILD = ADV3,    __RC__  )

      ! Doing Convection
      IF (do_ctmConvection) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = CONV, SRC_ID = ADV3, __RC__  )

         IF ( (TRIM(metType) == 'MERRA2') .OR. &
              (TRIM(metType) == 'FPIT')   .OR. &
              (TRIM(metType) == 'FP') ) THEN
               CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'LWI '/), &
                    DST_ID = CONV, SRC_ID = ECTM, __RC__  )
         ELSEIF ( TRIM(metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'ZLE ', 'LWI '/), &
                    DST_ID = CONV, SRC_ID = ECTM, __RC__  )
         END IF

         CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'ConvTR'/),          &
               CHILD = CONV,    __RC__  )
      END IF

      ! Doing Diffusion
      IF (do_ctmDiffusion) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = DIFF, SRC_ID = ADV3, __RC__  )

         IF ( (TRIM(metType) == 'MERRA2') .OR. &
              (TRIM(metType) == 'FPIT')   .OR. &
              (TRIM(metType) == 'FP') ) THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS'/), &
                    DST_ID = DIFF, SRC_ID = ECTM, __RC__  )
         ELSEIF ( TRIM(metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/'PLE ', 'MASS', 'ZLE '/), &
                    DST_ID = DIFF, SRC_ID = ECTM, __RC__  )
         END IF

         CALL MAPL_TerminateImport    ( GC,    &
               SHORT_NAME = (/'DiffTR'/),          &
               CHILD = DIFF,    __RC__  )
      END IF

      IF (enable_pTracers) THEN
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = PTRA, SRC_ID = ADV3, __RC__  )

         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'PLE'/), &
                 DST_ID = PTRA, SRC_ID = ECTM, __RC__  )
      ELSE
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/'AREA'/), &
                 DST_ID = CHEM, SRC_ID = ADV3, __RC__  )

         ! This is done for MERRA1, MERRA2, FPIT, FP
         CALL MAPL_AddConnectivity ( GC, &
                 SHORT_NAME  = (/ 'LFR   ', 'QCTOT ', 'TH    ', 'PLE   ', &
                                  'LWI   ', 'CNV_QC', &
                                  'BYNCY ', 'ITY   ', 'QICN  ', 'QLCN  ' /), &
                 DST_ID = CHEM, SRC_ID = ECTM, __RC__  )

         ! Additional Connectivity if using MERRA1
         IF ( TRIM(metType) == 'MERRA1') THEN
            CALL MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/ 'ZLE   ', 'RH2   ' /), &
                    DST_ID = CHEM, SRC_ID = ECTM, __RC__  )
         END IF
      END IF

      ! Create grid for this GC
      !------------------------
      !call MAPL_GridCreate  (GC, __RC__ )

      call MAPL_GenericSetServices ( GC, RC=STATUS )
      VERIFY_(STATUS)

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
      character(len=ESMF_MAXSTR)          ::  iNAME
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      character(len=ESMF_MAXSTR)          :: IAm = "Initialize"
      real(REAL64), allocatable           :: ak(:),bk(:)
      real(REAL64)                        :: ptop, pint
      integer                             :: counts(3),lm,ls


      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------

      !call ESMF_GridCompGet ( GC, name=COMP_NAME, GRID=GRID, RC=STATUS )
      call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
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
      call ESMF_AttributeSet(grid,name='GridAK', itemCount=LM+1, &
           valuelist=ak,rc=status)
      VERIFY_(STATUS)
      call ESMF_AttributeSet(grid,name='GridBK', itemCount=LM+1, &
           valuelist=bk,rc=status)
      VERIFY_(STATUS)
      deallocate(ak,bk)

      ! Get my MAPL_Generic state
      !--------------------------

      call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
      VERIFY_(STATUS)

      ! Call Initialize for every Child
      !--------------------------------
      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(STATE,"TOTAL")
      call MAPL_TimerOn(STATE,"INITIALIZE")

#ifdef PRINT_STATES
      call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
      if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
      call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
      if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif


      ! Get children and their im/ex states from my generic state.
      !----------------------------------------------------------

      call MAPL_Get ( STATE, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
      VERIFY_(STATUS)

      ! Extract the friendly tracers
      !-----------------------------
      IF (enable_pTracers) THEN
         if (do_ctmDiffusion) then
            !------------------
            ! Diffusion Tracers
            !------------------
            call ESMF_StateGet   (GIM(DIFF),  'DiffTR' , BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

            call MAPL_GridCompGetFriendlies(GCS(PTRA), "TURBULENCE", BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Diffusion Tracer Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

            ! Count tracers
            !--------------

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
            VERIFY_(STATUS)

            ! Get the names of all tracers to fill other turbulence bundles.
            !---------------------------------------------------------------

            allocate(NAMES(NUM_TRACERS),STAT=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  RC=STATUS)
            VERIFY_(STATUS)
         end if

         if (do_ctmConvection) then
            !-------------------
            ! Convection Tracers
            !-------------------
            call ESMF_StateGet       (GIM(CONV), 'ConvTR', BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

            call MAPL_GridCompGetFriendlies(GCS(PTRA), "MOIST", BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Convective Transport Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
            VERIFY_(STATUS)
         end if

         !----------------
         ! AdvCore Tracers
         !----------------
         call ESMF_StateGet       (GIM(ADV3), 'TRADV', BUNDLE, RC=STATUS )
         VERIFY_(STATUS)
  
         call MAPL_GridCompGetFriendlies(GCS(PTRA), "DYNAMICS", BUNDLE, RC=STATUS )
         VERIFY_(STATUS)

#ifdef PRINT_STATES
         call WRITE_PARALLEL ( trim(Iam)//": AdvCore Tracer Bundle" )
         if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

         call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
         VERIFY_(STATUS)

         ! Get the names of all tracers to fill other turbulence bundles.
         !---------------------------------------------------------------

         allocate(NAMES(NUM_TRACERS),STAT=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  RC=STATUS)
         VERIFY_(STATUS)

      ELSE 
         if (do_ctmDiffusion) then
            !------------------
            ! Diffusion Tracers
            !------------------
            call ESMF_StateGet   (GIM(DIFF),  'DiffTR' , BUNDLE, RC=STATUS )
            VERIFY_(STATUS)
  
            call MAPL_GridCompGetFriendlies(GCS(CHEM), "TURBULENCE", BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Diffusion Tracer Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

            ! Count tracers
            !--------------

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
            VERIFY_(STATUS)

            ! Get the names of all tracers to fill other turbulence bundles.
            !---------------------------------------------------------------

            allocate(NAMES(NUM_TRACERS),STAT=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  RC=STATUS)
            VERIFY_(STATUS)
         end if

         if (do_ctmConvection) then
            !-------------------
            ! Convection Tracers
            !-------------------
            call ESMF_StateGet       (GIM(CONV), 'ConvTR', BUNDLE, RC=STATUS )
            VERIFY_(STATUS)
  
            call MAPL_GridCompGetFriendlies(GCS(CHEM), "MOIST", BUNDLE, RC=STATUS )
            VERIFY_(STATUS)

#ifdef PRINT_STATES
            call WRITE_PARALLEL ( trim(Iam)//": Convective Transport Bundle" )
            if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

            call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
            VERIFY_(STATUS)

            IF (NUM_TRACERS .EQ. 0) THEN
               IF ( MAPL_am_I_root() ) THEN
                  PRINT*, '======================================='
                  PRINT*, '-----> No tracer friendly to MOIST'
                  PRINT*, '-----> Convection will not be performed'
                  PRINT*, '======================================='
               END IF
               do_ctmConvection = .FALSE.
            END IF
         end if

         !----------------
         ! AdvCore Tracers
         !----------------
         call ESMF_StateGet       (GIM(ADV3), 'TRADV', BUNDLE, RC=STATUS )
         VERIFY_(STATUS)

         call MAPL_GridCompGetFriendlies(GCS(CHEM), "DYNAMICS", BUNDLE, RC=STATUS )
         VERIFY_(STATUS)

#ifdef PRINT_STATES
         call WRITE_PARALLEL ( trim(Iam)//": AdvCore Tracer Bundle" )
         if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

         call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
         VERIFY_(STATUS)
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


      ! Get the target components name and set-up traceback handle.
      ! -----------------------------------------------------------

      call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // "::" // TRIM(Iam)

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(STATE,"TOTAL")
      call MAPL_TimerOn(STATE,"RUN")

      ! Get the children`s states from the generic state
      !-------------------------------------------------

      call MAPL_Get ( STATE,   &
          GCS=GCS, GIM=GIM, GEX=GEX,       &
          IM = IM, JM = JM, LM = LM,       &
          GCNames = GCNames,               &
          INTERNAL_ESMF_STATE = INTERNAL,  &
                               RC=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:" , RC=STATUS)
      VERIFY_(STATUS)

      !---------------------
      ! Cinderella Component: to derive variables for other components
      !---------------------
      I=ECTM

      call MAPL_TimerOn (STATE,GCNames(I))
      call ESMF_GridCompRun (GCS(I),               &
                             importState = GIM(I), &
                             exportState = GEX(I), &
                             clock       = CLOCK,  &
                             userRC      = STATUS  )
      VERIFY_(STATUS)
      call MAPL_TimerOff(STATE,GCNames(I))

      !--------
      ! advCore
      !--------
      IF (do_ctmAdvection) THEN
         I=ADV3

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                              importState = GIM(I), &
                              exportState = GEX(I), &
                              clock       = CLOCK,  &
                              userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      END IF

      !-----------
      ! Convection
      !-----------
      IF (do_ctmConvection) THEN
         I=CONV

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                exportState = GEX(I), &
                                   clock       = CLOCK,  &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      END IF

      IF (.NOT. enable_pTracers) THEN
         !----------
         ! Chemistry: Phase 1
         !----------
         I=CHEM   

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
      IF (do_ctmDiffusion) THEN
         I=DIFF

         call MAPL_TimerOn (STATE,GCNames(I))
         call ESMF_GridCompRun (GCS(I),               &
                                importState = GIM(I), &
                                   exportState = GEX(I), &
                                clock       = CLOCK,  &
                                userRC      = STATUS  )
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,GCNames(I))
      END IF

      IF (enable_pTracers) THEN
         !---------------
         ! Passive Tracer
         !---------------
         I=PTRA

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
         I=CHEM   

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
