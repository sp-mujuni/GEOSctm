#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CTM_rasCalculationsMod
!
! !INTERFACE:
!
      MODULE CTM_rasCalculationsMod
!
! !USES:
!
      use ESMF
      use MAPL_Mod
      use GEOS_Mod
      use GEOS_UtilsMod

      IMPLICIT NONE

      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: INIT_RASPARAMS
      PUBLIC  :: DO_RAS
!
! !PUBLIC DATA MEMBERS:
!
      !=================================================================
      ! Derived type for RASPARAMS
      !=================================================================
      TYPE, PUBLIC  :: RASPARAM_TYPE
            REAL    :: CUFRICFAC
            REAL    :: SHR_LAMBDA_FAC
            REAL    :: QC_CRIT_CN
            REAL    :: RASAL1
            REAL    :: RASAL2
            REAL    :: RASNCL
            REAL    :: LAMBDA_FAC
            REAL    :: LAMBMX_FAC
            REAL    :: MIN_DIAMETER
            REAL    :: CUFRICLAMBDA
            REAL    :: RDTLEXPON
            REAL    :: STRAPPING
            REAL    :: SDQV2
            REAL    :: SDQV3
            REAL    :: SDQVT1
            REAL    :: ACRITFAC
            REAL    :: HMINTRIGGER
            REAL    :: LLDISAGGXP
            REAL    :: PBLFRAC
            REAL    :: RASAUTORAMPB
            REAL    :: AUTOC_CN_ZDEP
            REAL    :: MAXDALLOWED_S
            REAL    :: MAXDALLOWED_D
            INTEGER :: RASAL_EXP
            REAL    :: RAS_RHMIN
            REAL    :: RAS_RHFULL
      END TYPE RASPARAM_TYPE
!

   character(len=32), parameter, public :: RasCodes(-2:7) = (/ &
         "Default                         ", & ! -2
         "Invalid Code                    ", & ! -1
         "Cloud detraining here is active ", & ! 0
         "PBL h < layer's h*              ", & ! 1
         "No Valid lambda                 ", & ! 2
         "Lambda out of bounds            ", & ! 3
         "A < Acrit                       ", & ! 4
         "Negative Kernel                 ", & ! 5
         "Invalid Code                    ", & ! 6
         "RH Trigger not met              "  & ! 7
         /)
!
!EOP
!-------------------------------------------------------------------------
      CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Initializes RASPARAMS with default values.
!\\
!\\
!
      SUBROUTINE INIT_RASPARAMS(RASPARAMS )
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(RASPARAM_Type), INTENT(INOUT) :: RASPARAMS
! 
! !REVISION HISTORY:
! 28 Feb 2017 - K. Yu - Initial Version
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Set Default values based on strippedinfo.txt

      RASPARAMS%CUFRICFAC      =      1.000
      RASPARAMS%SHR_LAMBDA_FAC =      0.05
      RASPARAMS%QC_CRIT_CN     =      8.0e-4
      RASPARAMS%RASAL1         =   1800.0
      RASPARAMS%RASAL2         = -43200.0
      RASPARAMS%RASNCL         =   -300.0
      RASPARAMS%LAMBDA_FAC     =      4.0
      RASPARAMS%LAMBMX_FAC     =      0.0
      RASPARAMS%MIN_DIAMETER   =    200.0
      RASPARAMS%CUFRICLAMBDA   =      7.5e-4
      RASPARAMS%RDTLEXPON      =      1.0
      RASPARAMS%STRAPPING      =     -1.0
      RASPARAMS%SDQV2          =      1.3
      RASPARAMS%SDQV3          =      1.3
      RASPARAMS%SDQVT1         =    263.0
      RASPARAMS%ACRITFAC       =      0.5
      RASPARAMS%HMINTRIGGER    =      1.0
      RASPARAMS%LLDISAGGXP     =      0.0
      RASPARAMS%PBLFRAC        =      0.1
      RASPARAMS%RASAUTORAMPB   =      0.8
      RASPARAMS%MAXDALLOWED_S  =   4000.0
      RASPARAMS%MAXDALLOWED_D  =       RASPARAMS%MAXDALLOWED_S
      RASPARAMS%RASAL_EXP      =      1
      RASPARAMS%RAS_RHMIN      =      0.5
      RASPARAMS%RAS_RHFULL     =      0.65
      RASPARAMS%AUTOC_CN_ZDEP  =      1.0

      END SUBROUTINE INIT_RASPARAMS
!EOC
!-------------------------------------------------------------------------
!BOP
!
      SUBROUTINE DO_RAS( RASPARAMS, PREF, TH, PLE, KH, PBLH, Q, FRLAND, TS, &
                         CNV_FLXD, CNV_FLXC, LATS, LONS, DT, IM, JM, LM)
!
! !INPUT PARAMETERS:
!
      INTEGER,             INTENT(IN) :: IM, JM, LM
      REAL,                INTENT(IN) :: DT                 ! model time step (s)
      real,                INTENT(IN) :: LONS  (IM,JM)      ! array of longitudes (rad)
      real,                INTENT(IN) :: LATS  (IM,JM)      ! array of latitudes  (rad)
      real,                INTENT(IN) :: FRLAND(IM,JM)      ! fraction of land
      real,                INTENT(IN) :: TS    (IM,JM)      ! surface temperature
      real,                INTENT(IN) :: PBLH  (IM,JM)      ! PBL height [m]
      real,                INTENT(IN) :: PREF  (      0:LM) ! reference_air_pressure
      real,                INTENT(IN) :: TH    (IM,JM,  LM) ! potential temperature
      real,                INTENT(IN) :: Q     (IM,JM,  LM) ! Specific humidity
      real,                INTENT(IN) :: PLE   (IM,JM,0:LM) ! air pressure
      real,                INTENT(IN) :: KH    (IM,JM,0:LM) ! scalar_diffusivity
      TYPE(RASPARAM_Type), INTENT(IN) :: RASPARAMS          ! parameters
!
! !OUTPUT PARAMETERS:
      real,               INTENT(OUT) :: CNV_FLXD(IM,JM,  LM) ! detraining mass flux
      real,               INTENT(OUT) :: CNV_FLXC(IM,JM,0:LM) ! cumulative mass flux 
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER                         :: IDIM, I, L, J, RC
      INTEGER                         :: IRUN
      INTEGER                         :: K0
      INTEGER                         :: ICMIN
      INTEGER                         :: KCBLMIN
      REAL                            :: PMIN_DET
      REAL                            :: PMIN_CBL
      integer, dimension(IM,JM,2)     :: SEEDRAS
      real,    dimension(IM,JM,  LM)  :: TEMP, PLO, PK, QSS, DQS
      real,    dimension(IM,JM,0:LM)  :: CNV_PLE
      real,    dimension(IM,JM,  LM)  :: ZLO
      real,    dimension(IM,JM,0:LM)  :: ZLE
      real,    dimension(IM,JM,0:LM)  :: PKE
      integer, dimension(IM,JM)       :: IRAS, JRAS, KCBL
      real,    dimension(      0:LM)  :: SIGE
      INTEGER, dimension(IM,JM)       :: KPBL
      REAL,    dimension(IM,JM)       :: KPBLIN
      real,    dimension(IM,JM,  LM)  :: WGT0, WGT1
      real,    dimension(IM,JM)       :: ZCBLx
      real,    dimension(IM,JM)       :: TPERT, QPERT, RASAL2_2d
      real,    dimension(IM,JM)       :: MXDIAMx
      real,    dimension(IM,JM,  LM)  :: CNV_FLX
      integer, dimension(IM,JM,  LM)  :: IRC
      ! Convective Fraction
      real   , dimension(IM,JM)       :: QV600, QSSFC
      real   , dimension(IM,JM)       :: CNV_FRACTION
      real                            :: CNV_Q600_MIN
      real                            :: CNV_Q600_MAX
      real,    dimension(IM,JM,  LM)  :: Q1, TH1
      INTEGER                         :: STATUS, CBL_METHOD
      !INTEGER                         :: levs600
      REAL                            :: CBL_QPERT, CBL_TPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND
      character(len=ESMF_MAXSTR)      :: IAm = "DO_RAS"
      REAL,     dimension(IM,JM,0:LM) :: tempZLE
      REAL :: qmin, qmax

      !<------ Beginning
      !-----------------

      IDIM     = IM*JM
      IRUN     = IM*JM
      K0       = LM

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% Set local variables to zero for safety's sake; at least the ones
      !%%% that aren't defined further down in the code. (bmy, 9/5/18)

      ! Integers
      CBL_METHOD   = 0
      ICMIN        = 0
      KCBLMIN      = 0
      SEEDRAS      = 0
      IRAS         = 0
      JRAS         = 0
      KCBL         = 0
      KPBL         = 0
      IRC          = 0
      STATUS       =-1

      ! Floating-point reals
      TH1          = 0.0
      Q1           = 0.0
      CNV_PLE      = 0.0
      PLO          = 0.0
      PKE          = 0.0
      PK           = 0.0
      TEMP         = 0.0
      ZLE          = 0.0
      ZLO          = 0.0
      SIGE         = 0.0
      ZCBLx        = 0.0
      WGT0         = 0.0
      WGT1         = 0.0
      TPERT        = 0.0
      QPERT        = 0.0
      QSSFC        = 0.0
      QSS          = 0.0
      DQS          = 0.0
      CNV_FRACTION = 0.0
      CNV_FLX      = 0.0
      CNV_FLXD     = 0.0
      CNV_FLXC     = 0.0
      RASAL2_2d    = 0.0
      MXDIAMx      = 0.0

      TH1      = TH                         ! local potential temperature
      Q1       = Q                          ! local specific humidity
      CNV_PLE  = PLE*.01                    ! convert to mb
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      TEMP     = TH*PK                      ! temperature

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

      PMIN_DET   =  3000.0
      PMIN_CBL   = 50000.0

      !--------------
      ! Compute ICMIN
      !--------------
      ICMIN    = max(1,count(PREF < PMIN_DET))
      KCBLMIN  =       count(PREF < PMIN_CBL)

      !----------------
      ! Compute SEEDRAS
      !----------------
      SEEDRAS(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
      SEEDRAS(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

      !----------------------
      ! Compute IRAS and JRAS
      !----------------------
      IRAS       = nint(LONS*100)
      JRAS       = nint(LATS*100)

      !-------------
      ! Compute SIGE
      !-------------
      SIGE(0:LM) = PREF(0:LM)/PREF(LM)

      !-------------------------------------------
      ! Compute the Planetary Boundary Layer Level
      !-------------------------------------------
      ! This piece of code was based on email exchanges between 
      ! Jules and Andrea on January 9, 2018.
      DO L=0,LM
         tempZLE(:,:,L) = ZLE(:,:,L) - ZLE(:,:,LM)
      ENDDO

      DO J = 1, JM
         DO I = 1, IM
            DO L = LM-1,0,-1
               IF ((PBLH(I,J) .LE. tempZLE(I,J,L))      .AND. &
                   (PBLH(I,J) .GT. tempZLE(I,J,L+1)))   THEN
                  KPBLIN(I,J) = L
                  CYCLE
               ENDIF
            ENDDO
         ENDDO
      ENDDO

           !---> This is the GEOS_MoistGridComp formulation
      KPBL = FINDPBL( KH, IM, JM, LM )

           !---> This is the Harvard formulation
!      DO J = 1,JM
!         DO I = 1,IM
!            IF (NINT(KPBLIN(I,J)).NE.0) THEN
!               KPBL(I,J) = MIN(NINT(KPBLIN(I,J))+1,LM-1)
!            ELSE
!               KPBL(I,J) = LM-1
!            ENDIF
!         ENDDO
!      ENDDO

      !-------------------------
      ! Compute KCBL, WGT0, WGT1
      !-------------------------

      CBL_METHOD = 6

      do J=1,JM
         do I=1,IM
            KCBL(I,J)              = KPBL(I,J)
            KCBL(I,J)              = MAX( KCBL(I,J), KCBLMIN )

            WGT0(I,J,:)            = 0.
            WGT0(I,J,KCBL(I,J):K0) = 1.0
            WGT1(I,J,:)            = 0.
            WGT1(I,J,KCBL(I,J):K0) = 1.0
 
            !ZCBLx(I,J)             = ZLE( I, J, KCBL(I,J)-1 )
            ZCBLx(I,J)             = ZLE( I, J, KCBL(I,J) )    ! Harvard formulation
         end do
      end do

      !---------------------
      ! Compute CNV_FRACTION
      !---------------------
      IF ( .FALSE. ) THEN
         ! Find QV at 600mb level
              !---> This is the GEOS_MoistGridComp formulation
         !call VertInterp(QV600,Q,log(PLE),log(60000.),STATUS)
         !VERIFY_(STATUS)
         !! Fill undefs (600mb below the surface) with surface QV values L=LM
         !levs600  = max(1,count(PREF < 60000.))
         !WHERE (QV600 == MAPL_UNDEF)
         !   QV600 = Q(:,:,levs600)
         !END WHERE

              !---> This is the Harvard formulation
         DO J = 1, JM
            DO I = 1, IM
               DO L = 1, LM
                  QV600(I,J) = Q(I,J,L)
                  IF (CNV_PLE(I,J,L) > 600.0) EXIT
               ENDDO
           ENDDO
         ENDDO

         !    CNV_FRACTION = 0.0  --->  Large-scale
         !    CNV_FRACTION = 1.0  --->  Deep-Convection
         CNV_FRACTION = 0.0

         CNV_Q600_MIN = 0.00250
         CNV_Q600_MAX = 0.00600

         if ( CNV_Q600_MAX > CNV_Q600_MIN ) then
            DO J=1, JM
               DO I=1, IM
                  CNV_FRACTION(I,J) =MAX(0.0,MIN(1.0,(QV600(I,J)-CNV_Q600_MIN)/(CNV_Q600_MAX-CNV_Q600_MIN)))
               END DO
            END DO
         endif
      ENDIF

           !---> This is the Harvard formulation
      CNV_FRACTION = 1.0 ! KNJR This is was set in the Harvard case

      !------------------------
      ! Compute TPERT and QPERT
      !------------------------

      CBL_QPERT       = 1.0 ! 1.0
      CBL_TPERT       =-2.0 !-2.0
      CBL_TPERT_MXOCN = 3.0 ! 0.035 ! 3.0
      CBL_TPERT_MXLND = 3.0 ! 0.035 ! 3.0

           !---> This is the GEOS_MoistGridComp formulation
      !QSSFC    = GEOS_QSAT( TS , CNV_PLE(:,:,LM) )
           !---> This is the Harvard formulation
      QSSFC = 0.0   ! for now

      TPERT  = ABS(CBL_TPERT) * ( TS - ( TEMP(:,:,LM)+ MAPL_GRAV*ZLO(:,:,LM)/MAPL_CP )  )

           !---> This is the Harvard formulation 
           !---> (section commented out when CNV_FRACTION = 1.0)
      !IF (CBL_TPERT .LT. 0) THEN
      !   ! Make TPERT 0 in areas of deep convection
      !   TPERT = TPERT*(1.0-CNV_FRACTION)
      !ENDIF

      QPERT  = CBL_QPERT * ( QSSFC - Q(:,:,LM) )
      TPERT  = MAX( TPERT , 0.0 )
      QPERT  = MAX( QPERT , 0.0 )

      where (FRLAND<0.1)
         TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
      elsewhere
         TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
      end where

      !------------------
      ! Compute RASAL2_2d
      !------------------

      if ( RASPARAMS%RASAL2 > 0.0) then
         RASAL2_2d(:,:) =  RASPARAMS%RASAL2
      else
         ! include CNV dependence
         RASAL2_2d(:,:) = CNV_FRACTION(:,:)*ABS(RASPARAMS%RASAL2) + (1.0-CNV_FRACTION(:,:))*RASPARAMS%RASAL1
      endif
      
      !---------------------
      ! Compute QSS and DSS
      !---------------------
           !---> This is the GEOS_MoistGridComp formulation
      !DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)
           !---> This is the Harvard formulation
      DO L = 1, LM
         DO J = 1, JM
            DO I = 1, IM
               CALL QSAT(TEMP(I,J,L), PLO(I,J,L), QSS(I,J,L), DQS(I,J,L), .TRUE.)
            END DO
         ENDDO
      ENDDO


      !-------------------------
      ! Call the main subroutine
      !-------------------------

      CALL RASE(                         &
                IDIM,                    &
                IRUN,                    &
                K0,                      &
                ICMIN,                   &
                DT ,                     &
                MAPL_CP,                 &   ! CPO
                MAPL_ALHL,               &   ! ALHLO
                MAPL_GRAV,               &   ! GRAVO
                SEEDRAS,                 &
                IRAS,                    &
                JRAS,                    &
                SIGE,                    &
                ! Inputs for CBL
                KCBL,                    &
                WGT0,                    &
                WGT1,                    &
                ZCBLx,                   &
                MXDIAMx,                 &
                TPERT,                   &
                QPERT,                   &
                ! Inputs
                TH1,                     &
                Q1,                      &
                QSS,                     &
                DQS,                     &
                CNV_FRACTION,            &
                RASAL2_2d,               &
                CNV_PLE,                 &
                PKE,                     &
                ! Outputs
                CNV_FLX,                 &    ! -> progno_clod
                CNV_FLXD,                &    ! -> diag
                CNV_FLXC,                &    ! -> progno_clod
                RASPARAMS,               &
                IRC                      )

      RETURN

      END SUBROUTINE DO_RAS
!EOC
!-------------------------------------------------------------------------
!BOP


   SUBROUTINE RASE(IDIM, IRUN, K0, ICMIN, DT ,            &
         CPO,ALHLO,GRAVO,                                 &
         SEEDRAS,IRAS,JRAS,SIGE,                          &
         KCBL,WGT0,WGT1,ZCBL,MXDIAM,TPERT,QPERT,          &
         THO, QHO,                                        & 
         QSS, DQS, CNV_FRACTION, RASAL2_2d,               &
         !pko, plo, phio, phie, qlo, qio,                  &
         PLE, PKE, FLX, FLXD, FLXC,                       &
         !ENTLAM,                                          &
         RASPARAMS,                                       &
         IRC                                              )


      !*********************************************************************
      !*********************************************************************
      !******************** Relaxed Arakawa-Schubert ***********************
      !************************ Parameterization ***************************
      !********************** SCALAR RAS-1 VERSION  ************************
      !************************* 31 DECEMBER 1999 **************************
      !*********************************************************************
      !************************** Developed By *****************************
      !*********************************************************************
      !************************ Shrinivas Moorthi **************************
      !******************************* and *********************************
      !************************** Max J. Suarez ****************************
      !*********************************************************************
      !******************** Laboratory for Atmospheres *********************
      !****************** NASA/GSFC, Greenbelt, MD 20771 *******************
      !*********************************************************************
      !*********************************************************************

      !  Input:
      !  ------
      ! 
      !     K0      : Number of vertical levels (increasing downwards)
      !
      !     DT      : Time step in seconds
      !
      !     RASAL   : Array of dimension K-1 containing relaxation parameters
      !               for cloud-types detraining at those levels
      !
      !     CPO     : Specific heat at constant pressure (J/kg/K)
      !
      !     ALHLO   : Latent Heat of condensation (J/kg)
      !
      !     GRAVO   : Acceleration due to gravity (m/s^2)
      !
      !     PLE     : 2D array of dimension (IDIM,K0+1) containing pressure
      !               in hPa at the interfaces of K-layers from top of the 
      !               atmosphere to the bottom  (mb)
      !
      !     PKE     : 2D array of dimension (IDIM,K0+1) containing (PRS/P00) **
      !               RKAP.  i.e. Exner function at layer edges.
      !
      !     PKL     : 2D array of dimension (IDIM,K0) ) containing the
      !               Exner function at the layers.
      !
      !     QSS     : 2D array of dimension (IDIM,K0  ) containing the
      !               saturation specific humidity at the layers. (kg/kg)
      !
      !     DQS     : 2D array of dimension (IDIM,K0  ) containing
      !               d(qss)/dt at the layers.  (1/K)
      !   
      !     CNV_FRACTION    : 1D array of dimension (IDIM) containing
      !               fraction of grid cell considered to be convective
      !   
      !  Update:
      !  -------
      !
      !     THO     : 2D array of dimension (IDIM,K0) containing potential
      !               temperature (K)
      !
      !     QHO     : 2D array of dimension (IDIM,K0) containing specific
      !               humidity (kg/kg)
      !
      !  Output:
      !  -------
      !!
      !     FLX     : 2D array of dimension (IDIM,K0) containing the
      !               cloud-base mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXD    : 2D array of dimension (IDIM,K0) containing the
      !               detrained  mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXC    : 2D array of dimension (IDIM,K0+1) containing the
      !               total cloud mass flux for all cloud types through
      !               the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
      !               and  FLXD(L) = FLXC(L+1)-FLSD(L) )
      !                (kg/m^2/s) 
      !
      !
      !************************************************************************

      !-->----------------
      !--> INPUT VARIABLES
      !-->----------------
      INTEGER,                     INTENT(IN   ) ::  IDIM
      INTEGER,                     INTENT(IN   ) ::  IRUN
      INTEGER,                     INTENT(IN   ) ::  K0     ! Number of vertical levels (increasing downwards)
      INTEGER,                     INTENT(IN   ) ::  ICMIN
      REAL,                        INTENT(IN   ) ::  DT     ! Time step in seconds
      REAL,                        INTENT(IN   ) ::  CPO    ! Specific heat at constant pressure (J/kg/K)
      REAL,                        INTENT(IN   ) ::  ALHLO  ! Latent Heat of condensation (J/kg)
      REAL,                        INTENT(IN   ) ::  GRAVO  ! Acceleration due to gravity (m/s^2)
      REAL, DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  PLE    ! pressure in hPa at the interfaces of K-layers
                                                            ! from top of the atmosphere to the bottom  (mb)
      REAL, DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  PKE    ! (PRS/P00) ** RKAP.  i.e. Exner function at layer edges.
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  QSS    ! saturation specific humidity at the layers (kg/kg)
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  DQS    ! d(qss)/dt at the layers.  (1/K)
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  CNV_FRACTION ! fraction of grid cell considered to be convective
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  RASAL2_2d
      REAL, DIMENSION (     K0+1), INTENT(IN   ) ::  SIGE
      INTEGER, DIMENSION (IDIM,2), INTENT(IN   ) ::  SEEDRAS
      INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  IRAS
      INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  JRAS
      INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  KCBL
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  ZCBL
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  TPERT
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  QPERT
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT0
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT1
      type (RASPARAM_TYPE),        INTENT(IN   ) ::  RASPARAMS

      !-->-----------------------
      !--> INPUT/OUTPUT VARIABLES
      !-->-----------------------
      REAL, DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  THO    ! potential temperature (K)
      REAL, DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  QHO    ! specific humidity (kg/kg)

      !-->-----------------
      !--> OUTPUT VARIABLES
      !-->-----------------
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  FLX    ! cloud-base mass flux for each cloud type ordered by
                                                            ! detrainment level (kg/m^2/s)
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  FLXD   ! detrained  mass flux for each cloud type ordered by
                                                            ! detrainment level (kg/m^2/s)
      REAL, DIMENSION (IDIM,K0+1), INTENT(  OUT) ::  FLXC   ! total cloud mass flux for all cloud types through
                                                            ! the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
                                                            ! and  FLXD(L) = FLXC(L+1)-FLSD(L) ) (kg/m^2/s)
      REAL, DIMENSION (IDIM     ), INTENT(  OUT) ::  MXDIAM
      INTEGER, DIMENSION(IDIM,K0), INTENT(  OUT) ::  IRC

      !-->----------------
      !--> LOCAL VARIABLES
      !-->----------------

      REAL,  DIMENSION(K0) :: POI_SV, QOI_SV
      REAL,  DIMENSION(K0) :: POI, QOI, DQQ, BET, GAM
      REAL,  DIMENSION(K0) :: POI_c, QOI_c
      REAL,  DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI
      REAL,  DIMENSION(K0) :: TCU, QCU, UCU, VCU,  CLN, RNS, POL,DM
      REAL,  DIMENSION(K0) :: QST, SSL,  RMF, RNN, RN1, RMFC, RMFP
      REAL,  DIMENSION(K0) :: GMS, ETA, GMH, EHT,  HCC, RMFD
      REAL,  DIMENSION(K0) :: HOL, HST, QOL, ZOL, HCLD
      REAL,  DIMENSION(K0) :: LAMBDSV, BKE , CVW, UPDFRC
      REAL,  DIMENSION(K0) :: RASAL, MTKWI, UPDFRP,BK2,BK3

      REAL,  DIMENSION(K0+1) :: PRJ, PRS, QHT, SHT ,ZET, XYD, XYD0

      INTEGER,  DIMENSION(K0-1) :: RC

      INTEGER :: K,MY_PE

      REAL, DIMENSION(IDIM,K0) :: LAMBDSV2


      REAL TX2, TX3, UHT, VHT, AKM, ACR, ALM, TTH, QQH, SHTRG, DQX
      REAL WFN, TEM, TRG, TRGEXP, EVP, WLQ, QCC,MTKW_MAX 
      REAL SHTRG_FAC, SIGE_MINHOL, WFNOG

      INTEGER I, IC, L, ICL , ITR , ICL_C, N_DTL
      INTEGER NDTLEXPON

      INTEGER,  ALLOCATABLE, DIMENSION(:) :: ICL_V

      !  RASE GLOBAL CONSTANTS


      REAL GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP, OBG, AFC 

!!!!!!!!!
      REAL PBLFRAC,AUTORAMPB,CO_ZDEP
      REAL RASAL1, RASAL2, RASAL2i, CO_T, RASNCL,FRICLAMBDA
      REAL LAMBDA_FAC,ACRITFAC
      REAL LAMBMX_FAC, DIAMMN_MIN,RDTLEXPON, CLI_CRIT, MAXDALLOWED_D, MAXDALLOWED_S
      REAL RHMN, RHMX
      INTEGER RASAL_EXP

      real cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!

      REAL, PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
      !      REAL, PARAMETER :: PBLFRAC = 0.5
      REAL, PARAMETER :: RHMAX   = 0.9999

      !  LAMBDA LIMITS
      REAL            :: LAMBDA_MIN
      REAL            :: LAMBDA_MAX

      !  TRIGGER PARAMETERS
      real, parameter :: RHO_W  =  1.0e3  ! Density of liquid water in kg/m^3

      !character(len=ESMF_MAXSTR)          :: CBL_STYLE

      real*8, dimension(k0) :: tcu8, qcu8, pcu, flx8
      real*8    :: cup  !, dpd, tla
      logical :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft

      real*8, dimension(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, trcfac
      real*8, dimension(k0) :: ALFIND,ALFINT,ALFINQ,RHC_LS
      real*8, dimension(k0+1) :: prs8, phih8
      real*8 :: FRACBL, dt8, rasalf
      integer :: KPBL

      !real*8 :: MAX_NEG_BOUY = 1.0 ! no inhibition for =1.0
      !!real*8 :: ALFINT = 0.5
      !!real*8 :: ALFINQ = 0.5
      !real*8 :: RHFACL = 0.0 ! not used
      !real*8 :: RHFACS = 0.0 ! no inhibition
      !!real*8 :: ALFIND = 1.0
      !!real*8 :: RHC_LS = 0.80
      !real*8 :: dsfc   = 0.001
      !real*8 :: cd     = 1.e-3
      !real*8 :: wfnc   = 0.0
      !real*8 :: tla    = -1.0
      !real*8 :: dpd    = 300.

      ! *********************************************************************

      IF(IRUN <= 0) RETURN

      !ENTLAM    =0.
      IRC = -2

      ! Zero out local varaiables for safety's sake 
      ! INTEGER
      RC            = 0
      K             = 0
      MY_PE         = 0
      I             = 0
      IC            = 0
      L             = 0
      ICL           = 0
      ITR           = 0
      ICL_C         = 0
      N_DTL         = 0
      NDTLEXPON     = 0
      KPBL          = 0

      !REAL
      POI_SV        = 0.0
      QOI_SV        = 0.0
      POI           = 0.0
      QOI           = 0.0
      DQQ           = 0.0
      BET           = 0.0
      GAM           = 0.0
      POI_c         = 0.0
      QOI_c         = 0.0
      PRH           = 0.0
      PRI           = 0.0
      GHT           = 0.0
      DPT           = 0.0
      DPB           = 0.0
      PKI           = 0.0
      TCU           = 0.0
      QCU           = 0.0
      UCU           = 0.0   !<---- not used
      VCU           = 0.0   !<---- not used
      CLN           = 0.0
      RNS           = 0.0
      POL           = 0.0
      DM            = 0.0
      QST           = 0.0
      SSL           = 0.0
      RMF           = 0.0   !<---- This init was necessary
      RNN           = 0.0   !<---- not used
      RN1           = 0.0   !<---- not used
      RMFC          = 0.0
      RMFP          = 0.0
      GMS           = 0.0
      ETA           = 0.0
      GMH           = 0.0
      EHT           = 0.0
      HCC           = 0.0
      RMFD          = 0.0   !<---- This init was necessary
      HOL           = 0.0
      HST           = 0.0
      QOL           = 0.0
      ZOL           = 0.0
      HCLD          = 0.0   !<---- not used
      LAMBDSV       = 0.0
      BKE           = 0.0   !<---- not used
      CVW           = 0.0   !<---- not used
      UPDFRC        = 0.0
      RASAL         = 0.0
      MTKWI         = 0.0   !<---- not used
      UPDFRP        = 0.0   !<---- not used
      BK2           = 0.0   !<---- not used
      BK3           = 0.0   !<---- not used
      PRJ           = 0.0   !<---- This init was necessary
      PRS           = 0.0   !<---- This init was necessary
      QHT           = 0.0
      SHT           = 0.0
      ZET           = 0.0
      XYD           = 0.0   !<---- not used
      XYD0          = 0.0   !<---- not used
      LAMBDSV2      = 0.0
      TX2           = 0.0
      TX3           = 0.0   !<---- not used
      UHT           = 0.0   !<---- not used
      VHT           = 0.0   !<---- not used
      AKM           = 0.0
      ACR           = 0.0
      ALM           = 0.0
      TTH           = 0.0   !<---- not used
      QQH           = 0.0   !<---- not used
      SHTRG         = 0.0   !<---- not used
      DQX           = 0.0   !<---- not used?
      WFN           = 0.0
      TEM           = 0.0
      TRG           = 0.0
      TRGEXP        = 0.0   !<---- not used
      EVP           = 0.0   !<---- not used
      WLQ           = 0.0   !<---- only place set
      QCC           = 0.0
      MTKW_MAX      = 0.0
      SHTRG_FAC     = 0.0   !<---- not used
      SIGE_MINHOL   = 0.0
      WFNOG         = 0.0
      cld_radius    = 0.0   !<---- not used
      areal_frac    = 0.0   !<---- not used
      spect_mflx    = 0.0   !<---- not used
      cvw_cbase     = 0.0   !<---- not used
      LAMBDA_MIN    = 0.0
      LAMBDA_MAX    = 0.0

      ! REAL*8
      tcu8          = 0.0d0
      qcu8          = 0.0d0
      pcu           = 0.0d0
      flx8          = 0.0d0
      cup           = 0.0d0
      toi8          = 0.0d0
      qoi8          = 0.0d0
      prsm8         = 0.0d0
      phil8         = 0.0d0
      qli8          = 0.0d0
      qii8          = 0.0d0
      trcfac        = 0.0d0
      ALFIND        = 0.0d0
      ALFINT        = 0.0d0
      ALFINQ        = 0.0d0
      RHC_LS        = 0.0d0
      prs8          = 0.0d0
      phih8         = 0.0d0
      FRACBL        = 0.0d0
      dt8           = 0.0d0
      rasalf        = 0.0d0

      ! LOGICALS
      revap         = .FALSE.
      wrkfun        = .FALSE.
      calkpb        = .FALSE.
      crtfun        = .FALSE.
      lprnt         = .FALSE.
      dndrft        = .FALSE.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      RASAL1        = RASPARAMS%RASAL1
      RASAL2        = RASPARAMS%RASAL2
      RASNCL        = RASPARAMS%RASNCL
      LAMBDA_FAC    = RASPARAMS%LAMBDA_FAC
      LAMBMX_FAC    = RASPARAMS%LAMBMX_FAC
      DIAMMN_MIN    = RASPARAMS%MIN_DIAMETER
      RDTLEXPON     = RASPARAMS%RDTLEXPON
      ACRITFAC      = RASPARAMS%ACRITFAC
      PBLFRAC       = RASPARAMS%PBLFRAC
      AUTORAMPB     = RASPARAMS%RASAUTORAMPB
      CO_ZDEP       = RASPARAMS%AUTOC_CN_ZDEP
      MAXDALLOWED_S = RASPARAMS%MAXDALLOWED_S
      MAXDALLOWED_D = RASPARAMS%MAXDALLOWED_D
      RASAL_EXP     = RASPARAMS%RASAL_EXP
      RHMN          = RASPARAMS%RAS_RHMIN
      RHMX          = RASPARAMS%RAS_RHFULL

      GRAV  = GRAVO
      ALHL  = ALHLO
      CP    = CPO
      CPI   = 1.0/CP      
      ALHI  = 1.0/ALHL
      GRAVI = 1.0/GRAV
      CPBG  = CP*GRAVI
      DDT   = DAYLEN/DT
      AFC   = -1.04E-4*SQRT(DT*113.84)
      LBCP  = ALHL*CPI
      OBG   = 100.*GRAVI

      DO I=1,IRUN

         K = KCBL(I)

         rc(icmin) = 0

         CALL FINDDTLS

         IF ( K > 0 ) THEN 
            CALL STRAP( FINAL=0 )
            call HTEST
           
            DO ICL_C = 1,N_DTL
               ICL   = ICL_V( ICL_C )

               IF( ICL > ICMIN ) THEN
                  CALL CLOUDE(ICL)
               ENDIF

            ENDDO

            IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

               CALL STRAP( FINAL=1 )

            ELSE

               CALL STRAP( FINAL=2 )

            ENDIF

         ELSE 

            CALL STRAP( FINAL=2 )

         ENDIF

      ENDDO

      IF(ALLOCATED(ICL_V)) DEALLOCATE(ICL_V)

      RETURN

   CONTAINS

      !*********************************************************************

      SUBROUTINE CLOUDE(IC)

         INTEGER, INTENT(IN ) :: IC
         REAL :: DEEP_FACT,WSCALE

         REAL :: CLI , TE_A, C00_X, CLI_CRIT_X, PETE, TOKI, GMHx, HSTx  !, dQx
         REAL :: DT_LYR, RATE, CLOSS, F2, F3, F4, F5
         INTEGER :: K700

         !=============================AER_CLOUD local variables ====================
         REAL :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
	      alph_e, beta_e, RH_AMB, ECRIT

         INTEGER :: INX, naux

         ALM   = 0.
         TRG   = AMIN1(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))

         F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  ! F4 should ramp from 0 at SIG=AUTORAMPB
         ! to 1 at SIG=AUTORAMPB-0.2

         if ( SIGE(IC) >= 0.5 ) then
            F5 = 1.0
         else
            F5 = 1.0 - 2.*CO_ZDEP *( 0.5 - SIGE(IC) )
            F5 = MAX( F5 , 0.0 )
         endif


         IF(TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>
            RC(IC) = 7
            RETURN
         ENDIF

         !  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL

         POI_c = POI
         QOI_c = QOI
         POI_c(K) =  POI_c(K) + TPERT(I)
         QOI_c(K) =  QOI_c(K) + QPERT(I)

         ZET(K+1) = 0.
         SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
         DO L=K,IC,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI_c(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
         ENDDO

         DO L=IC+1,K
            TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
            SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
            QHT(L)  = .5*(QOL(L)+QOL(L-1))
         ENDDO

         !  CALCULATE LAMBDA, ETA, AND WORKFUNCTION

         LAMBDA_MIN = .2/MXDIAM(I)
         LAMBDA_MAX = .2/  200. 


         IF (HOL(K) <= HST(IC)) THEN   ! CANNOT REACH IC LEVEL  ======>>
            RC(IC) = 1
            RETURN
         ENDIF

         !  LAMBDA CALCULATION: MS-A18

         TEM  =       (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
         DO L=IC+1,K-1
            TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
         ENDDO

         IF(TEM <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
            RC(IC) = 2
            RETURN
         ENDIF

         ALM     = (HOL(K)-HST(IC)) / TEM

         IF(ALM > LAMBDA_MAX) THEN
            RC(IC) = 3
            RETURN
         ENDIF

         TOKI=1.0

         IF(ALM < LAMBDA_MIN) THEN
       TOKI = ( ALM/LAMBDA_MIN )**2  !we can probably replace this by a actual distribution based on grid cell size

            !RC(IC) = 6
            !RETURN
         ENDIF

         !  ETA CALCULATION: MS-A2

         DO L=IC+1,K
            ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
         ENDDO
         ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))

         !  WORKFUNCTION CALCULATION:  MS-A22

         WFN     = 0.0
         HCC(K)  = HOL(K)
         DO L=K-1,IC+1,-1
            HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
            TEM    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
            EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
            WFN    = WFN + (TEM - EHT(L)*HST(L))*GAM(L)
         ENDDO
         HCC(IC) = HST(IC)*ETA(IC)
         WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)

         !   ALPHA CALCULATION 
         RASAL2i = RASAL2_2d(I)
         IF ( ZET(IC) <  2000. ) RASAL(IC) = RASAL1
         IF ( ZET(IC)  >= 2000. ) THEN 
   !WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
            RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*MIN(1.0,(ZET(IC) - 2000.)/8000.)**RASAL_EXP
         ENDIF
   !WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
         RASAL(IC) = DT / RASAL(IC)

         !  TEST FOR CRITICAL WORK FUNCTION

         CALL ACRITN(POL(IC), PRS(K), ACR)

         IF(WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
            RC(IC) = 4
            RETURN
         ENDIF


         !     CALCULATE GAMMAS AND KERNEL

         GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O GRAV)
         GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL     ! MS-A31 (W/O GRAV)
         AKM    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O GRAV)

         TX2     = GMH(K)
         DO L=K-1,IC+1,-1
            GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  ))   & 
                  + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
            GMH(L) = GMS(L)                         &
                  + ( ETA(L  )*(QHT(L)-QOL(L  ))   &
                  + ETA(L+1)*(QOL(L)-QHT(L+1)) )*ALHL*PRI(L)
            TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
            AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)

         ENDDO

         GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
         AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

         GMH(IC) =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL &
               + ETA(IC  )*(HST(IC)-HOL(IC  ))     )*PRI(IC)

         !    CLOUD BASE MASS FLUX

         IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN  !  =========>
            RC(IC) = 5
            RETURN
         ENDIF

         WFN = - (WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
         WFN = MIN( ( RASAL(IC)*TRG*TOKI )*WFN  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))



         !    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT

         WFNOG    = WFN*GRAVI
         TEM      = WFN*GRAVI
         RMF (IC) = RMF (IC) +     TEM           ! (kg/m^2/step)
         RMFD(IC) = RMFD(IC) +     TEM * ETA(IC) ! (kg/m^2/step)

         DO L=IC+1,K
            RMFP(L) = TEM * ETA(L)                 ! (kg/m^2/step)
            RMFC(L) = RMFC(L)  +  RMFP(L)          ! (kg/m^2/step)

         ENDDO

         !    THETA AND Q CHANGE DUE TO CLOUD TYPE IC

         DO L=IC,K
            GMH(L) = GMH(L) * WFN
            GMS(L) = GMS(L) * WFN
            QOI(L) = QOI(L) + (GMH(L) - GMS(L)) * ALHI
            POI(L) = POI(L) + GMS(L)*PKI(L)*CPI
            QST(L) = QST(L) + GMS(L)*BET(L)*CPI
         ENDDO

         RC(IC) = 0

         RETURN

      END SUBROUTINE CLOUDE

      SUBROUTINE ACRITN( PL, PLB, ACR)

         REAL, INTENT(IN ) :: PL, PLB
         REAL, INTENT(OUT) :: ACR

         INTEGER IWK

         !!REAL, PARAMETER :: FACM=0.5

         REAL, PARAMETER :: &
               PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
               550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

         REAL, PARAMETER :: & 
               A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
               0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
               0.0553, 0.0445, 0.0633     /)   !!*FACM

         IWK = PL * 0.02 - 0.999999999

         IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
            ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
         ELSEIF(IWK > 15) THEN
            ACR = A(15)
         ELSE
            ACR = A(1)
         ENDIF

         ACR = ACRITFAC  * ACR * (PLB - PL)

         RETURN

      END SUBROUTINE ACRITN



      SUBROUTINE HTEST
         REAL,  DIMENSION(K0) :: HOL1
         integer  :: lminhol
         real     :: minhol


         hol=0.            ! HOL initialized here in order not to confuse Valgrind debugger
         lminhol  = K+1
         MINHOL   = -999999.
         ZET(K+1) = 0
         SHT(K+1) = CP*POI(K)*PRJ(K+1)
         DO L=K,ICMIN,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
         ENDDO

         HOL1=HOL
         DO L=K-1,ICMIN+1,-1
            HOL1(L)  = (0.25*HOL(L+1)+0.50*HOL(L)+0.25*HOL(L-1))
            if ( ( MINHOL>=HOL1(L) ) .OR. (MINHOL<0.) ) THEN
               MINHOL   =  HOL1(L)
               LMINHOL  =  L
            end if
         ENDDO

         SIGE_MINHOL = SIGE(LMINHOL)

      end subroutine HTEST


      subroutine FINDDTLS
         real :: SIGDT0,sigmax,sigmin
         integer :: LL

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% GCST UPDATE:
         !%%% Dimension The_Seed with a large number of integers (100)
         !%%% even though we probably won't use all elements.
         !%%% Also remove platform-specific code. (bmy, 8/30/18)
         INTEGER :: The_Seed(100)
         INTEGER :: Seed_Size

         !%%% Zero out all elements of The_Seed
         The_Seed = 0


         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         if(THE_SEED(1) == 0) THE_SEED(1) =  5
         if(THE_SEED(2) == 0) THE_SEED(2) = -5

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% GCST UPDATE:
         !%%% First call Random_Seed to get the size,
         !%%% then only use that many elements of The_Seed.
         !%%% Otherwise this throws a compile error with gfortran 8.2.0.
         !%%% (bmy, 8/30/18)
         CALL Random_Seed( SIZE = Seed_Size              )
         call Random_Seed( PUT  = The_Seed(1:Seed_Size)  )
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         SIGMAX=SIGE(K)
         SIGMIN=SIGE(ICMIN)

         if ( RASNCL < 0.0 ) then 
            !! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
            N_DTL = K - ICMIN 
         else 
            N_DTL = min( int( RASNCL ) , K-ICMIN )
         endif

         if(allocated(ICL_V)) deallocate(ICL_V)
         allocate(ICL_V(N_DTL))

         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
            do L=1,N_DTL
               ICL_V(L) = ICMIN + L - 1
            enddo
         else if ( RASNCL < -100.0 ) then 
            do L=1,N_DTL
               ICL_V(L) = K - L
               !! NO SHALLOW CONV           ICL_V(L) = 56 - L
            enddo
         else
            do L=1,N_DTL
               call random_number ( SIGDT0 )
               SIGDT0 = 1.00 - ( SIGDT0**RDTLEXPON )
               SIGDT0 = SIGMIN+SIGDT0*(SIGMAX-SIGMIN)

               do LL=ICMIN,K
                  if ( (SIGE(LL+1)>=SIGDT0) .and. (SIGE(LL)<SIGDT0 ) ) ICL_V(L)=LL
               enddo
            end do
         endif

      end subroutine FINDDTLS

      SUBROUTINE STRAP(FINAL)

         INTEGER :: FINAL
         REAL , DIMENSION(K0)  :: WGHT, MASSF

         REAL :: WGHT0, PRCBL

         integer, parameter :: nrands=1
         real ::    rndu(nrands)
         integer :: seedcbl(nrands)

         ! !DESCRIPTION: 
         !   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
         !   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
         !   to redistribute convective tendencies within CBL

         integer :: KK

         !  LOCAL VARIABLES FOR USE IN CLOUDE

         !!IF (.NOT. PRESENT(FINAL)) THEN
         IF (FINAL==0) THEN

            do kk=icmin,k+1
               PRJ(kk) = PKE(I,kk)
            enddo

            poi=0.        ! These initialized here in order not to confuse Valgrind debugger
            qoi=0.        ! Do not believe it actually makes any difference.

            PRS(ICMIN:K0+1) = PLE(I,ICMIN:K0+1)
            POI(ICMIN:K)   = THO(I,ICMIN:K)
            QOI(ICMIN:K)   = QHO(I,ICMIN:K)

            QST(ICMIN:K) = QSS(I,ICMIN:K)
            DQQ(ICMIN:K) = DQS(I,ICMIN:K)

!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
            MASSF(:) = WGT0(I,:)

!!! RESET PRESSURE at bottom edge of CBL 
            PRCBL = PRS(K)
            do l= K,K0
               PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
            end do
            PRS(K+1) = PRCBL
            PRJ(K+1) = (PRS(K+1)/1000.)**(MAPL_RGAS/CP)

            DO L=K,ICMIN,-1
               POL(L)  = 0.5*(PRS(L)+PRS(L+1))
               PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) &
                     / (ONEPKAP*(PRS(L+1)-PRS(L)))
               PKI(L)  = 1.0 / PRH(L)
               DPT(L)  = PRH(L  ) - PRJ(L)
               DPB(L)  = PRJ(L+1) - PRH(L)
               PRI(L)  = .01 / (PRS(L+1)-PRS(L))
            ENDDO

!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
            if( K<=K0) then
               POI(K) = 0.
               QOI(K) = 0.

               !! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
               WGHT = 0.
               DO L=K,K0
                  WGHT(L)   = MASSF(L) *                &
                        ( PLE(I,L+1) - PLE(I,L) ) &
                        /( PRS(K+1)   - PRS(K)  )
               END DO

               DO L=K,K0
                  POI(K) = POI(K) + WGHT(L)*THO(I,L)
                  QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
               ENDDO

               CALL QSAT( POI(K)*PRH(K) , POL(K), QST(K), DQQ(K), .FALSE. )

            endif

            rndu(:) = max( seedras(I,1)/1000000., 1e-6 )

            MXDIAM(I) = CNV_FRACTION(I)*ABS(MAXDALLOWED_D) + (1-CNV_FRACTION(I))*ABS(MAXDALLOWED_S)
            if (MAXDALLOWED_D > 0) then
               ! Make MXDIAM stochastic
               MXDIAM(I) = MXDIAM(I)*( rndu(1)**(-1./2.) )
            endif

            DO L=K,ICMIN,-1
               BET(L)  = DQQ(L)*PKI(L)  !*
               GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
               IF(L<K) THEN
                  GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
               ENDIF
            ENDDO

            TCU(ICMIN:K) = -POI(ICMIN:K)*PRH(ICMIN:K)
            QCU(ICMIN:K) = -QOI(ICMIN:K)

            RMF  = 0.
            RMFD = 0.
            RMFC = 0.
            RMFP = 0.

            POI_SV = POI
            QOI_SV = QOI

         END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (FINAL==1) THEN 

            THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
            QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)

            !! De-strap tendencies from RAS
            !! specify weighting "SHAPE"
            WGHT   = WGT1(I,:)

            !! Scale properly by layer masses
            wght0 = 0.
            DO L=K,K0 
               wght0 = wght0 + WGHT(L)* ( PLE(I,L+1) - PLE(I,L) )
            END DO

            wght0 = ( PRS(K+1)   - PRS(K)  )/wght0

            WGHT  = wght0 * WGHT


            DO L=K,K0 
               THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
               QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
            END DO

            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
            FLXD(I,ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)

            FLX (I,1:ICMIN-1) = 0.
            FLXD(I,1:ICMIN-1) = 0.
            FLXC(I,1:ICMIN-1) = 0.


            IF ( K < K0 ) THEN 
               FLX (I,K:K0) = 0.
               FLXD(I,K:K0) = 0.
               FLXC(I,K:K0) = 0.
            END IF
          
            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

            IF(ALLOCATED(ICL_V)) DEALLOCATE( ICL_V )

         ENDIF

         IF (FINAL==2) THEN 


            FLX (I,:) = 0.
            FLXD(I,:) = 0.
            FLXC(I,:) = 0.

            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

         ENDIF

         RETURN

      END SUBROUTINE STRAP


   END SUBROUTINE RASE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function FINDPBL( KH, IM, JM, LM ) result( KPBL )
    ! !DESCRIPTION: 
    integer                    , intent(in) :: IM,JM,LM
    real, dimension(IM,JM,0:LM), intent(in) :: KH

    integer, dimension(IM,JM)               :: KPBL

    integer                                 :: I, J, L
    real                                    :: KHCRIT

    KHCRIT = 2.0  ! m+2 s-1

    do I = 1, IM
       do J = 1, JM

          KPBL(I,J) = LM
          do L = LM-1,1,-1
             IF( ( KH(I,J,L) >= KHCRIT ).AND.( KH(I,J,L-1) < KHCRIT ) ) THEN ! "top" is between L and L-1
                KPBL(I,J) = L+1   ! returned index for CBL q,t etc. is just below PBL top
                EXIT
             ENDIF
          enddo

          KPBL(I,J)=MIN( LM-1, KPBL(I,J) )

       end do
    end do


  end function FINDPBL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

      subroutine qsat (TT,P,Q,DQDT,LDQDT)
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Compute Saturation Specific Humidity
!
!  INPUT:
!  ======
!    TT ......... Temperature (Kelvin)
!    P .......... Pressure (mb)
!    LDQDT ...... Logical Flag to compute QSAT Derivative
!
!  OUTPUT:
!  =======
!    Q .......... Saturation Specific Humidity
!    DQDT ....... Saturation Specific Humidity Derivative wrt Temperature
!
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      IMPLICIT NONE
      REAL TT, P, Q, DQDT
      LOGICAL LDQDT
      REAL AIRMW, H2OMW
      
      PARAMETER ( AIRMW  = 28.97      )                                         
      PARAMETER ( H2OMW  = 18.01      )                                         

      REAL ESFAC, ERFAC
      PARAMETER ( ESFAC = H2OMW/AIRMW       )
      PARAMETER ( ERFAC = (1.0-ESFAC)/ESFAC )

      real aw0, aw1, aw2, aw3, aw4, aw5, aw6
      real bw0, bw1, bw2, bw3, bw4, bw5, bw6
      real ai0, ai1, ai2, ai3, ai4, ai5, ai6
      real bi0, bi1, bi2, bi3, bi4, bi5, bi6

      real d0, d1, d2, d3, d4, d5, d6
      real e0, e1, e2, e3, e4, e5, e6
      real f0, f1, f2, f3, f4, f5, f6
      real g0, g1, g2, g3, g4, g5, g6

! ********************************************************
! ***  Polynomial Coefficients WRT Water (Lowe, 1977) ****
! ***              (Valid +50 C to -50 C)             ****
! ********************************************************

      parameter ( aw0 =  6.107799961e+00 * esfac )
      parameter ( aw1 =  4.436518521e-01 * esfac )
      parameter ( aw2 =  1.428945805e-02 * esfac )
      parameter ( aw3 =  2.650648471e-04 * esfac )
      parameter ( aw4 =  3.031240396e-06 * esfac )
      parameter ( aw5 =  2.034080948e-08 * esfac )
      parameter ( aw6 =  6.136820929e-11 * esfac )

      parameter ( bw0 = +4.438099984e-01 * esfac )
      parameter ( bw1 = +2.857002636e-02 * esfac )
      parameter ( bw2 = +7.938054040e-04 * esfac )
      parameter ( bw3 = +1.215215065e-05 * esfac )
      parameter ( bw4 = +1.036561403e-07 * esfac )
      parameter ( bw5 = +3.532421810e-10 * esfac )
      parameter ( bw6 = -7.090244804e-13 * esfac )


! ********************************************************
! ***   Polynomial Coefficients WRT Ice  (Lowe, 1977) ****
! ***              (Valid  +0 C to -50 C)             ****
! ********************************************************

      parameter ( ai0 = +6.109177956e+00 * esfac )
      parameter ( ai1 = +5.034698970e-01 * esfac )
      parameter ( ai2 = +1.886013408e-02 * esfac )
      parameter ( ai3 = +4.176223716e-04 * esfac )
      parameter ( ai4 = +5.824720280e-06 * esfac )
      parameter ( ai5 = +4.838803174e-08 * esfac )
      parameter ( ai6 = +1.838826904e-10 * esfac )

      parameter ( bi0 = +5.030305237e-01 * esfac )
      parameter ( bi1 = +3.773255020e-02 * esfac )
      parameter ( bi2 = +1.267995369e-03 * esfac )
      parameter ( bi3 = +2.477563108e-05 * esfac )
      parameter ( bi4 = +3.005693132e-07 * esfac )
      parameter ( bi5 = +2.158542548e-09 * esfac )
      parameter ( bi6 = +7.131097725e-12 * esfac )


! ********************************************************
! ***         Polynomial Coefficients WRT Ice         ****
! ***   Starr and Cox (1985) (Valid -40 C to -70 C)   ****
! ********************************************************


      parameter ( d0 = 0.535098336e+01 * esfac )
      parameter ( d1 = 0.401390832e+00 * esfac )
      parameter ( d2 = 0.129690326e-01 * esfac )
      parameter ( d3 = 0.230325039e-03 * esfac )
      parameter ( d4 = 0.236279781e-05 * esfac )
      parameter ( d5 = 0.132243858e-07 * esfac )
      parameter ( d6 = 0.314296723e-10 * esfac )

      parameter ( e0 = 0.469290530e+00 * esfac )
      parameter ( e1 = 0.333092511e-01 * esfac )
      parameter ( e2 = 0.102164528e-02 * esfac )
      parameter ( e3 = 0.172979242e-04 * esfac )
      parameter ( e4 = 0.170017544e-06 * esfac )
      parameter ( e5 = 0.916466531e-09 * esfac )
      parameter ( e6 = 0.210844486e-11 * esfac )


! ********************************************************
! ***         Polynomial Coefficients WRT Ice         ****
! ***   Starr and Cox (1985) (Valid -65 C to -95 C)   ****
! ********************************************************

      parameter ( f0 = 0.298152339e+01 * esfac )
      parameter ( f1 = 0.191372282e+00 * esfac )
      parameter ( f2 = 0.517609116e-02 * esfac )
      parameter ( f3 = 0.754129933e-04 * esfac )
      parameter ( f4 = 0.623439266e-06 * esfac )
      parameter ( f5 = 0.276961083e-08 * esfac )
      parameter ( f6 = 0.516000335e-11 * esfac )

      parameter ( g0 = 0.312654072e+00 * esfac )
      parameter ( g1 = 0.195789002e-01 * esfac )
      parameter ( g2 = 0.517837908e-03 * esfac )
      parameter ( g3 = 0.739410547e-05 * esfac )
      parameter ( g4 = 0.600331350e-07 * esfac )
      parameter ( g5 = 0.262430726e-09 * esfac )
      parameter ( g6 = 0.481960676e-12 * esfac )

      REAL        TMAX, TICE
      PARAMETER ( TMAX=323.15, TICE=273.16)
      
      REAL T, D, W, QX, DQX
      T = MIN(TT,TMAX) - TICE
      DQX = 0.
      QX  = 0.

! Fitting for temperatures above 0 degrees centigrade
! ---------------------------------------------------
      if(t.gt.0.) then
       qx = aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6)))))
      if (ldqdt)  then
      dqx = bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6)))))
      endif
      endif

! Fitting for temperatures between 0 and -40
! ------------------------------------------
      if( t.le.0. .and. t.gt.-40.0 ) then
        w = (40.0 + t)/40.0
       qx =     w *(aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6))))))  &
          + (1.-w)*(ai0+T*(ai1+T*(ai2+T*(ai3+T*(ai4+T*(ai5+T*ai6))))))
      if (ldqdt)  then
      dqx =     w *(bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6))))))  &
          + (1.-w)*(bi0+T*(bi1+T*(bi2+T*(bi3+T*(bi4+T*(bi5+T*bi6))))))
      endif
      endif

! Fitting for temperatures between -40 and -70
! --------------------------------------------
      if( t.le.-40.0 .and. t.ge.-70.0 ) then
       qx = d0+T*(d1+T*(d2+T*(d3+T*(d4+T*(d5+T*d6)))))
      if (ldqdt) then
      dqx = e0+T*(e1+T*(e2+T*(e3+T*(e4+T*(e5+T*e6)))))
      endif
      endif

! Fitting for temperatures less than -70
! --------------------------------------
      if(t.lt.-70.0) then
       qx = f0+t*(f1+t*(f2+t*(f3+t*(f4+t*(f5+t*f6)))))
      if (ldqdt) then
      dqx = g0+t*(g1+t*(g2+t*(g3+t*(g4+t*(g5+t*g6)))))
      endif
      endif

! Compute Saturation Specific Humidity
! ------------------------------------
      D = (P-ERFAC*QX)
      IF(D.LT.0.) THEN
       Q = 1.0
       IF (LDQDT)  DQDT = 0.
      ELSE
       D = 1.0 / D
       Q = MIN(QX * D,1.0)
       IF (LDQDT)  DQDT = (1.0 + ERFAC*Q) * D * DQX
      ENDIF
      RETURN

      end subroutine qsat
!EOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine VertInterp(v2,v3,ple,pp,rc)

    real    , intent(OUT) :: v2(:,:)
    real    , intent(IN ) :: v3(:,:,:)
    real    , intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer k,km
    logical edge

    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    ASSERT_(edge .or. size(v3,3)==km)

    v2   = MAPL_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
          where( (pp>pb.and.pp<=ple(:,:,km+1)) )
             v2 = v3(:,:,km)
          end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

!-------------------------------------------------------------------------

END MODULE CTM_rasCalculationsMod
