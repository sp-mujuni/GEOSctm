!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GenericConvectionMethod_mod 
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
      module GenericConvectionMethod_mod
!
! !USES:
      use ESMF
      use MAPL_Mod
      !USE Chem_UtilMod
      USE convectiveTransport_mod
      use GmiArrayBundlePointer_mod
      USE GmiESMFrcFileReading_mod
!
      implicit none
!
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
!
! !PUBLIC MEMBER FUNCTIONS:
      private 
      public  :: initializeGenericConvection
      public  :: runGenericConvection
      public  :: finalizeGenericConvection
!
! !PUBLIC DATA MEMBERS:
      public  :: genConvection_GridComp
!
      TYPE genConvection_GridComp
        character(len=ESMF_MAXSTR) :: name = "Generic Convection"

        ! Dimension
        INTEGER :: i1, i2, im, j1, j2, jm, km

        ! Surface area of grid cells
        REAL(KIND=DBL), POINTER :: cellArea(:,:)

        logical :: det_ent      ! flag for doing detrainment then entrainment
        logical :: do_downdraft ! flag for doing downdrafts

        logical, pointer :: isFixedConcentration(:)
        logical :: FIRST
      END TYPE genConvection_GridComp
!
      character(len=ESMF_MAXSTR) :: rcfilen = 'CTM_GridComp.rc' 
!
      REAL, PARAMETER :: mwtAir    = 28.9
      REAL, PARAMETER :: rStar     = 8.314E+03
      REAL, PARAMETER :: Pa2hPa    = 0.01
      REAL, PARAMETER :: ToGrPerKg = 1000.00
      REAL, PARAMETER :: secPerDay = 86400.00
!
! !DESCRIPTION:
! Generic convection for the convective transport only.
!
! !AUTHOR:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeGenericConvection
!
! !INTERFACE:
!
      subroutine initializeGenericConvection (self, impConv, expConv, nymd, nhms, &
                                      grid, tdt, rc)
!
! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: nymd, nhms     ! Time from AGCM
      REAL,    INTENT(IN) :: tdt            ! Chemistry time step (secs)
      TYPE(ESMF_GRID), INTENT(IN)  :: grid
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(out) ::  rc           ! Error return code:
                                            !  0 - all is well
                                            !  1 - 
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ESMF_State),   INTENT(INOUT)  :: impConv    ! Import State
      TYPE(ESMF_State),   INTENT(INOUT)  :: expConv    ! Export State
      type (genConvection_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Convection related variables from the resource file.
!
! !LOCAL VARIABLES:
      ! Grid cell area can be set by initialize
      ! ---------------------------------------
      REAL, POINTER, DIMENSION(:,:) :: cellArea
      integer ::  i1, i2, j1, j2, km, im, jm
      integer :: k1, k2
      integer                    :: numSpecies, ic, is
      integer                    :: STATUS
      type (ESMF_Config)         :: convConfigFile
      REAL, ALLOCATABLE :: var2D(:,:)
      character(len=ESMF_MAXSTR) :: varName, fileName
      character(len=ESMF_MAXSTR) :: IAm = "initializeGenericConvection"
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Get the GMI grid information
      i1   = self%i1
      i2   = self%i2
      im   = self%im
      j1   = self%j1
      j2   = self%j2
      jm   = self%jm
      km   = self%km

      k1   = 1
      k2   = km

      self%FIRST = .TRUE.

      !################################
      ! Begin reading the resource file
      !################################

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Starting Reading the Convection Resource File"
      ENDIF

      convConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(convConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)


      call rcEsmfReadLogical(convConfigFile, self%det_ent, &
     &               "det_ent:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(convConfigFile, self%do_downdraft, &
     &               "do_downdraft:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT*," -----> det_ent      = ", self%det_ent
         PRINT*," -----> do_downdraft = ", self%do_downdraft
         PRINT *,"Done Reading the Convection Resource File"
      END IF

      !##############################
      ! End reading the resource file
      !##############################

      ! Grid box surface area, m^{2}
      ! ----------------------------
      CALL MAPL_GetPointer(impConv, cellArea, 'AREA', rc=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(self%cellArea(i1:i2,j1:j2), STAT=STATUS)
      VERIFY_(STATUS)
      self%cellArea(i1:i2,j1:j2)=cellArea(:,:)

      return

      end subroutine initializeGenericConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! IROUTINE: runGenericConvection
!
! !INTERFACE:
!
      subroutine runGenericConvection (self, impConv, expConv, nymd, nhms,  &
                   tdt, enable_rasCalculations, rc)
!
! !USES:
!
! !INPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: impConv ! Import State
      INTEGER, INTENT(IN) :: nymd, nhms          ! time
      REAL,    INTENT(IN) :: tdt                 ! chemical timestep (secs)
      LOGICAL, INTENT(IN) :: enable_rasCalculations ! do RAS calculations?
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(OUT) ::  rc                ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: expConv   ! Export State
      type (genConvection_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Runs the convection component.
!
! !LOCAL VARIABLES:
      REAL, POINTER, DIMENSION(:,:)   :: zpbl
      REAL, POINTER, DIMENSION(:,:,:) :: ple, zle, totalMass
      REAL, POINTER, DIMENSION(:,:,:) :: CNV_MFC, CNV_MFD, T

      integer                    :: STATUS
      integer :: numSpecies, ic, kR, ik, is
      integer ::  i1, i2, j1, j2, km, im, jm
      integer :: k1, k2, ivert, ilong
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
      REAL(KIND=DBL), allocatable :: zmdu       (:, :, :)
      REAL(KIND=DBL), allocatable :: zmeu       (:, :, :)
      REAL(KIND=DBL), allocatable :: zmed       (:, :, :)
      REAL(KIND=DBL), allocatable :: zmmd       (:, :, :)
      REAL(KIND=DBL), allocatable :: zmmu       (:, :, :)
      REAL(KIND=DBL), allocatable :: hkdu       (:, :, :)
      REAL(KIND=DBL), allocatable :: hkeu       (:, :, :)
      REAL(KIND=DBL), allocatable :: hkmu       (:, :, :)
      REAL(KIND=DBL) :: cdt

      type (ESMF_Field)                   :: FIELD
      type (ESMF_Array)                   :: ARRAY
      type (ESMF_FieldBundle)             :: ConvTR
      REAL, POINTER, DIMENSION(:,:,:)     :: S

      type (ESMF_Config)         :: convConfigFile
      type (t_GmiArrayBundle), pointer :: concentration(:)
      character(len=ESMF_MAXSTR) :: IAm = "runGenericConvection"
      character(len=ESMF_MAXSTR) :: NAME, speciesName 
!EOP
!-------------------------------------------------------------------------
!BOC

      ! Get the GMI grid information
      i1   = self%i1
      i2   = self%i2
      im   = self%im
      j1   = self%j1
      j2   = self%j2
      jm   = self%jm
      km   = self%km

      k1   = 1
      k2   = km
      ivert= km
      ilong= i2-i1+1

      ! Get the bundles containing the quantities to be diffused, 
      !----------------------------------------------------------

      call ESMF_StateGet(impConv, 'ConvTR' ,    ConvTR,     RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(ConvTR, fieldCOUNT=numSpecies, RC=STATUS)
      VERIFY_(STATUS)

      ! Get the molecular weights of species to be convected
      ! Identify the "fixed" tracers
      ! Map the convected tracer indices into the GMI tracer indices
      !-------------------------------------------------------------

      IF (self%FIRST) THEN
         self%FIRST = .FALSE.

         allocate(self%isFixedConcentration(numSpecies), STAT=STATUS)
         VERIFY_(STATUS)
         self%isFixedConcentration(:) = .FALSE.

         DO ic = 1, numSpecies
            ! Get field and name from tracer bundle
            !--------------------------------------
            call ESMF_FieldBundleGet(ConvTR, ic, FIELD, RC=STATUS)
            VERIFY_(STATUS)
   
            call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
            VERIFY_(STATUS)
   
            !Identify fixed species such as O2, N2, ad.
            if (TRIM(NAME) == 'ACET' .OR. TRIM(NAME) == 'N2'   .OR. &
                TRIM(NAME) == 'O2'   .OR. TRIM(NAME) == 'NUMDENS')  THEN
               self%isFixedConcentration(ic) = .TRUE.
            end if

         END DO

         convConfigFile = ESMF_ConfigCreate(rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigLoadFile(convConfigFile, TRIM(rcfilen), rc=STATUS )
         VERIFY_(STATUS)

      END IF

      ! Get the tracers from the ESMF Bundle
      !-------------------------------------
      ALLOCATE(concentration(numSpecies), STAT=STATUS)
      VERIFY_(STATUS)

      DO ic = 1, numSpecies
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
      CALL MAPL_GetPointer(impConv,         T,       'T', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,      zpbl,    'ZPBL', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,       ple,     'PLE', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv, totalMass,    'MASS', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,       zle,     'ZLE', RC=STATUS); VERIFY_(STATUS)
      IF (enable_rasCalculations) THEN
         CALL MAPL_GetPointer(expConv,   CNV_MFD, 'CNV_MFD', RC=STATUS); VERIFY_(STATUS)
         CALL MAPL_GetPointer(expConv,   CNV_MFC, 'CNV_MFC', RC=STATUS); VERIFY_(STATUS)
      ELSE
         CALL MAPL_GetPointer(impConv,   CNV_MFD, 'CNV_MFD', RC=STATUS); VERIFY_(STATUS)
         CALL MAPL_GetPointer(impConv,   CNV_MFC, 'CNV_MFC', RC=STATUS); VERIFY_(STATUS)
      ENDIF

      allocate(pbl(i1:i2,j1:j2),                 STAT=STATUS); VERIFY_(STATUS)
      allocate(press3c(i1:i2,j1:j2,k1:k2),  STAT=STATUS); VERIFY_(STATUS)
      allocate(press3e(i1:i2,j1:j2,k1-1:k2),STAT=STATUS); VERIFY_(STATUS)
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

      cdt = tdt

      call doGenericConvectiveTransport (self%det_ent, self%do_downdraft, pbl, cmf, &
                      dtrain, eu, ed, md, gridBoxHeight, mass, kel, press3e,       &
                      concentration, self%isFixedConcentration, self%cellArea, cdt,  &
                      i1, i2, j1, j2, k1, k2, ilong, ivert, numSpecies)

      deallocate(pbl )
      deallocate(dtrain, cmf, kel, eu, ed, md)
      deallocate(press3e, press3c, mass, gridBoxHeight, pl)

      ! Pass bacck the tracers to the ESMF Bundle
      !------------------------------------------
      DO ic = 1, numSpecies
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


      return

      end subroutine RunGenericConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FinalizeGenericConvection
!
! !INTERFACE:
!
      subroutine FinalizeGenericConvection (self)
!
! !INPUT/OUTPUT PARAMETERS:
      type (genConvection_GridComp), intent(inout) :: self
!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR) :: IAm = "finalizeGenericConvection"
!EOP
!-------------------------------------------------------------------------
!BOC
      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Completed Convection"
      ENDIF

      return

      end subroutine FinalizeGenericConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doGenericConvectiveTransport
!
! !INTERFACE:
!
      subroutine doGenericConvectiveTransport (det_ent, do_downdraft, pbl, cldmas, &
                      dtrn, eu, ed, md, grid_height, mass, kel, press3e,       &
                      concentration, isFixedConcentration, mcor, tdt,          &
                      i1, i2, j1, j2, k1, k2, ilong, ivert, num_species)
!
! !INPUT PARAMETER:
      integer, intent(in) :: i1, i2, j1, j2, k1, k2
      integer, intent(in) :: ilong, ivert, num_species
      logical, intent(in) :: isFixedConcentration(:)
      logical, intent(in) :: det_ent ! flag for doing detrainment then entrainment
      logical, intent(in) :: do_downdraft ! flag for doing downdrafts
      real*8 , intent(in) :: tdt          ! model time step  (s)
      real*8 , intent(in) :: mcor       (i1:i2, j1:j2) ! area of each grid box (m^2)
      real*8 , intent(in) :: pbl        (i1:i2,j1:j2) ! planetary boundary layer thickness (m)
      real*8 , intent(in) :: cldmas     (i1:i2,j1:j2,k1:k2) ! convective mass flux in     updraft (kg/m^2/s)
      real*8 , intent(in) :: dtrn       (i1:i2,j1:j2,k1:k2) ! detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
      real*8 , intent(in) :: eu         (i1:i2,j1:j2,k1:k2) ! ntrainment into convective updraft (s^-1)
      real*8 , intent(in) :: ed         (i1:i2,j1:j2,k1:k2) ! entrainment into convective downdraft (s^-1)
      real*8 , intent(in) :: md         (i1:i2,j1:j2,k1:k2) ! convective mass flux in downdraft (kg/m^2/s)
      real*8 , intent(in) :: grid_height(i1:i2,j1:j2,k1:k2) ! grid box height  (m)
      real*8 , intent(in) :: mass       (i1:i2,j1:j2,k1:k2) ! mass of air in each grid box (kg)
      real*8 , intent(in) :: kel        (i1:i2,j1:j2,k1:k2) ! temperature      (degK)
      real*8 , intent(in) :: press3e    (i1:i2,j1:j2,k1-1:k2) ! atmospheric pressure at the edge of each grid box (mb)
!
! !INPUT/OUTPUT PARAMETERS:
                             ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DEFINED PARAMETERS:
      real*8, parameter :: MBSTH = 1.0d-15 ! threshold below which we treat 
                                           ! mass fluxes as zero (mb/s)
      real*8, parameter :: GMI_G  =   9.81d0 ! mean surface gravity accel. (m/s^2)
      real*8, parameter :: PASPMB = 100.00d0 ! pascals  per millibar
!
! !DESCRIPTION:
! This is the interface routine to Convective Transport.  
! It formats the gem variables to satisfy the Convective Transport routine.
!
!EOP
!------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: iku
      integer :: il, ij, ik, ic
      integer :: il2g                      ! gathered index to operate over
      integer :: itdt_conv
      integer :: num_conv_steps
      real*8  :: rnum_conv_steps
      real*8  :: tdt_conv
      real*8  :: xmbsth
      integer :: ideep(ilong)              ! gathering array
      integer :: pbli (ilong)              ! index of pbl height
      real*8  :: updraft_velocity(i1:i2)   ! velocity in convective updraft
                                           ! (m/s)
      real*8  :: dpi(ilong,  k1:k2)        ! delta pressure between interfaces
      real*8  :: dui(ilong,  k1:k2)        ! mass detraining from updraft
      real*8  :: eui(ilong,  k1:k2)        ! mass entraining into updraft
      real*8  :: mui(ilong,  k1:k2)        ! mass flux up
      real*8  :: mdi(ilong,  k1:k2)        ! mass flux down
      real*8  :: fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  :: qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture

      ideep(:) = 0
      pbli (:) = 0

      dpi(:,:) = 0.0d0
      dui(:,:) = 0.0d0
      eui(:,:) = 0.0d0
      mui(:,:) = 0.0d0
      mdi(:,:) = 0.0d0

      fracis(:,:,:) = 0.0d0

      xmbsth = MBSTH

      updraft_velocity(:) = 0.0d0

      ! -----------------------------------------------------------
      ! Calculate any needed sub-cycling of the convection operator
      ! by comparing the mass flux over the full time step to the
      ! mass within the grid box.
      ! -----------------------------------------------------------

      num_conv_steps = Maxval (tdt * cldmas(:,:,:) * Spread (mcor(:,:), 3, k2-k1+1) /  &
                               mass(:,:,:)) + 1.0d0

      rnum_conv_steps = num_conv_steps
      tdt_conv        = tdt / rnum_conv_steps

      IJLOOP: do ij = j1, j2
         il2g = 0

         ILLOOP: do il = i1, i2
            ! ----------------------------------------------------
            ! Verify that there is convection in the current cell.
            ! ----------------------------------------------------
            if (Maxval (cldmas(il,ij,:)) >= 0.0d0) then
               il2g = il2g + 1
               ideep(il2g) = il

               dpi(il2g,:) = (press3e(il,ij,k2-1:k1-1:-1) -  &
                              press3e(il,ij,k2:k1:-1)) * PASPMB

               mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G

               eui(il2g,k1+1:k2-1) = Max (0.0d0, cldmas(il,ij,k2-1:k1+1:-1) -  &
                                                 cldmas(il,ij,k2-2:k1  :-1) +  &
                                                 dtrn  (il,ij,k2-1:k1+1:-1))

               eui(il2g,:) = eui(il2g,:) * GMI_G

               if (det_ent) dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * GMI_G
            end if
         end do ILLOOP

         IL2GIF: if (il2g /= 0) then
            fracis(:,:,:) = 1.0d0

            ! ----------------------------------------------------------------
            ! Find the index of the top of the planetary boundary layer (pbl).
            ! Convection will assume well mixed tracers below that level.
            ! ----------------------------------------------------------------
            do il = 1, il2g
               pbli(il) = 0

               IKLOOP: do ik = k1, k2
                  if (pbl(ideep(il),ij) < Sum (grid_height(ideep(il),ij,k1:ik))) then
                     pbli(il) = ik
                     exit IKLOOP
                  end if
               end do IKLOOP

               if (pbli(il) == 0) then
                  err_msg = 'Could not find pbl in doGenericConvectiveTransport.'
                  PRINT*, err_msg
                  stop
               end if

               pbli(il) = k2 - pbli(il)
            end do

            ITDTCLOOP: do itdt_conv = 1, num_conv_steps
               do ic = 1, num_species
                  qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
               end do

               call convectiveTransport (il2g, tdt_conv, xmbsth, ideep,  &
                           pbli, dui, eui, mui, mdi, dpi, fracis, qq, &
                           isFixedConcentration, i1, i2, k1, k2, ilong, num_species)

               do ic = 1, num_species
                  concentration(ic)%pArray3D(:,ij,k1:k2) = qq(:,k2:k1:-1,ic)
               end do
            end do ITDTCLOOP
         end if IL2GIF
      end do IJLOOP

      return

      end subroutine doGenericConvectiveTransport
!EOC
!-----------------------------------------------------------------------------
  end module GenericConvectionMethod_mod
