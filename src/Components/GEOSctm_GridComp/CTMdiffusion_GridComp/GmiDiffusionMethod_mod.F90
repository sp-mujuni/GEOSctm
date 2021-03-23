!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDiffusionMethod_mod 
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
      module GmiDiffusionMethod_mod
!
! !USES:
      use ESMF
      use MAPL
      USE Chem_UtilMod
      USE updateDiffusion_mod
      use GmiArrayBundlePointer_mod
!
      implicit none
!
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
!
! !PUBLIC MEMBER FUNCTIONS:
      private 
      public  :: initializeDiffusion
      public  :: runDiffusion
      public  :: finalizeDiffusion
!
! !PUBLIC DATA MEMBERS:
      public  :: Diffusion_GridComp
!
      TYPE Diffusion_GridComp
        character(len=ESMF_MAXSTR) :: name = "GMI Diffusion"
        character(len=ESMF_MAXSTR) :: metType

        ! Dimension
        INTEGER :: i1, i2, im, j1, j2, jm, km

        ! Surface area of grid cells
        REAL(KIND=DBL), POINTER :: cellArea(:,:)

        integer :: diffu_opt        ! diffusion option
        real*8  :: vert_diffu_coef  ! Scalar vertical diffusion coefficient  (m^2/s)
        real*8  :: pbl_mixing_tau   ! PBL mixing tau
        logical, pointer :: isFixedConcentration(:)
        logical :: FIRST
      END TYPE Diffusion_GridComp
!
      REAL, PARAMETER :: mwtAir    = 28.9
      REAL, PARAMETER :: rStar     = 8.314E+03
      REAL, PARAMETER :: Pa2hPa    = 0.01
      REAL, PARAMETER :: ToGrPerKg = 1000.00
      REAL, PARAMETER :: secPerDay = 86400.00
!
! !DESCRIPTION:
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
! !IROUTINE: initializeDiffusion
!
! !INTERFACE:
!
      subroutine initializeDiffusion (self, impDiff, expDiff, nymd, nhms, &
                                      grid, tdt, rc)
!
! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: nymd, nhms     ! Time from AGCM
      REAL,    INTENT(IN) :: tdt            ! Chemistry time step (secs)
      type(ESMF_Grid)          :: grid
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(out) ::  rc           ! Error return code:
                                            !  0 - all is well
                                            !  1 - 
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ESMF_State),   INTENT(INOUT)  :: impDiff    ! Import State
      TYPE(ESMF_State),   INTENT(INOUT)  :: expDiff    ! Export State
      type (Diffusion_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Diffusion related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer ::  i1, i2, j1, j2, km, im, jm
      integer :: k1, k2, ilo, ihi, julo, jhi
      integer                    :: STATUS
      type (ESMF_Config)         :: diffConfigFile
      REAL, ALLOCATABLE :: var2D(:,:)
      ! Grid cell area can be set by initialize
      ! ---------------------------------------
      REAL, POINTER, DIMENSION(:,:) :: cellArea
      character(len=ESMF_MAXSTR) :: varName, fileName
      character(len=ESMF_MAXSTR) :: IAm = "initializeDiffusion"
      character(len=ESMF_MAXSTR) :: rcfilen = 'CTM_GridComp.rc' 
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
      ilo  = i1
      ihi  = i2
      julo = j1
      jhi  = j2

      !################################
      ! Begin reading the resource file
      !################################

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Starting Reading the Diffusion Resource File"
      ENDIF

      diffConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(diffConfigFile, TRIM(rcfilen), rc=STATUS )
      VERIFY_(STATUS)

      ! --------------------------------
      ! diffu_opt
      !   0:  no diffusion
      !   1:  do DAO2 vertical diffusion
      ! --------------------------------

      call ESMF_ConfigGetAttribute(diffConfigFile, self%diffu_opt, &
     &                label   = "diffu_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      ! --------------------------------
      ! metType
      !   MERRA2:  
      !   MERRA1:  
      !   FPIT:  
      !   FP:  
      ! --------------------------------

      call ESMF_ConfigGetAttribute(diffConfigFile, self%metType, &
     &                label   = "metType:",&
     &                default = 'MERRA2', rc=STATUS )
      VERIFY_(STATUS)

      ! pbl_mixing_tau : scalar vertical diffusion coefficient (m^2/s)

      call ESMF_ConfigGetAttribute(diffConfigFile, self%pbl_mixing_tau, &
     &                label   = "pbl_mixing_tau:",&
     &                default = 3600.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ! ---------------------------------------------------------------
      ! vert_diffu_coef : scalar vertical diffusion coefficient (m^2/s)
      ! ---------------------------------------------------------------

      call ESMF_ConfigGetAttribute(diffConfigFile, self%vert_diffu_coef, &
     &                label   = "vert_diffu_coef:",&
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------

      if ((self%diffu_opt .LT. 0) .OR. (self%diffu_opt .GT. 2)) then
         PRINT*,TRIM(IAm),': diffu_opt out of range'
         STATUS = 1
         VERIFY_(STATUS)
      end if

      ! Grid box surface area, m^{2}
      ! ----------------------------
      CALL MAPL_GetPointer(impDiff, cellArea, 'AREA', rc=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(self%cellArea(i1:i2,j1:j2), STAT=STATUS)
      VERIFY_(STATUS)
      self%cellArea(i1:i2,j1:j2)=cellArea(i1:i2,j1:j2)

      self%FIRST = .TRUE.

      return

      end subroutine initializeDiffusion
!EOC
!-------------------------------------------------------------------------
!BOP
!
! IROUTINE: runDiffusion
!
! !INTERFACE:
!
      subroutine runDiffusion (self, impDiff, expDiff, nymd, nhms,  &
     &              tdt, rc)
!
! !USES:
!
! !INPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: impDiff ! Import State
      INTEGER, INTENT(IN) :: nymd, nhms          ! time
      REAL,    INTENT(IN) :: tdt                 ! chemical timestep (secs)
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(OUT) ::  rc                ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: expDiff   ! Export State
      type (Diffusion_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Runs the diffusion component.
!
! !LOCAL VARIABLES:
      REAL, POINTER, DIMENSION(:,:)   :: tropp, zpbl
      REAL, POINTER, DIMENSION(:,:,:) :: KH, totalMass, ple, zle

      integer :: numSpecies, ic, STATUS, kR, k
      integer :: i1, i2, j1, j2, km, im, jm
      integer :: k1, k2, ilo, ihi, julo, jhi
      REAL, ALLOCATABLE :: pl(:,:,:)

      REAL(KIND=DBL), ALLOCATABLE :: tropopausePress(:,:)
      REAL(KIND=DBL), ALLOCATABLE :: pctm1(:,:)
      REAL(KIND=DBL), ALLOCATABLE :: pbl(:,:)
      REAL(KIND=DBL), ALLOCATABLE :: kzz(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: mass(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: press3c(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: press3e(:,:,:)
      REAL(KIND=DBL), ALLOCATABLE :: gridBoxHeight(:,:,:)

      REAL(KIND=DBL)              :: cdt

      type (ESMF_Field)               :: FIELD
      type (ESMF_Array)               :: ARRAY
      type (ESMF_FieldBundle)         :: DiffTR
      REAL, POINTER, DIMENSION(:,:,:) :: S

      type (t_GmiArrayBundle), pointer :: concentration(:)
      character(len=ESMF_MAXSTR) :: NAME
      character(len=ESMF_MAXSTR) :: IAm 
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "runDiffusion"

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
      ilo  = i1
      ihi  = i2
      julo = j1
      jhi  = j2

      ! Get the bundles containing the quantities to be diffused, 
      !----------------------------------------------------------

      call ESMF_StateGet(impDiff, 'DiffTR' , DiffTR,     RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(DiffTR, fieldCOUNT=numSpecies, RC=STATUS)
      VERIFY_(STATUS)

      ! Identify the "fixed" tracers
      !-----------------------------

      IF (self%FIRST) THEN
         allocate(self%isFixedConcentration(numSpecies), STAT=STATUS)
         VERIFY_(STATUS)
         self%isFixedConcentration(:) = .FALSE.
      END IF

      ALLOCATE(concentration(numSpecies), STAT=STATUS)
      VERIFY_(STATUS)

      DO ic = 1, numSpecies
         ! Get field and name from tracer bundle
         !--------------------------------------
         call ESMF_FieldBundleGet(DiffTR, ic, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

         !Identify fixed species such as O2, N2, ad.
         IF (self%FIRST) THEN
            IF (ic .EQ. numSpecies) self%FIRST = .FALSE.

            IF (TRIM(NAME) == 'ACET' .OR. TRIM(NAME) == 'N2'   .OR. &
                TRIM(NAME) == 'O2'   .OR. TRIM(NAME) == 'NUMDENS')  THEN
               self%isFixedConcentration(ic) = .TRUE.
            END IF
         END IF

         ! Get pointer to the quantity
         !----------------------------
         call ESMFL_BundleGetPointerToData(DiffTR, NAME, S, RC=STATUS)
         VERIFY_(STATUS)

         ! The quantity must exist; others are optional.
         !----------------------------------------------
         ASSERT_(associated(S ))

         ALLOCATE(concentration(ic)%pArray3d(i1:i2,j1:j2,km), STAT=STATUS)
         VERIFY_(STATUS)

         ! Do not forget to re-order the vertical levels
         concentration(ic)%pArray3d(:,:,km:1:-1) = S(:,:,:)
      END DO

      allocate(pbl(i1:i2, j1:j2))
      CALL MAPL_GetPointer(impDiff, zpbl, 'ZPBL', RC=STATUS)
      VERIFY_(STATUS)
      pbl(:,:) = zpbl(:,:)

      cdt = tdt

      if (self%diffu_opt == 1) then
         CALL MAPL_GetPointer(impDiff, KH,    'KH',    RC=STATUS); VERIFY_(STATUS)
         CALL MAPL_GetPointer(impDiff, ple,   'PLE',   RC=STATUS); VERIFY_(STATUS)
         CALL MAPL_GetPointer(impDiff, tropp, 'TROPP', RC=STATUS); VERIFY_(STATUS)

         allocate(pctm1(ilo:ihi,julo:jhi),          STAT=STATUS); VERIFY_(STATUS)
         allocate(tropopausePress(i1:i2,j1:j2),     STAT=STATUS); VERIFY_(STATUS)
         allocate(kzz(i1:i2,j1:j2,k1:k2),           STAT=STATUS); VERIFY_(STATUS)
         kzz = 0.0d0
         allocate(press3c(ilo:ihi,julo:jhi,k1:k2),  STAT=STATUS); VERIFY_(STATUS)
         allocate(press3e(ilo:ihi,julo:jhi,k1-1:k2),  STAT=STATUS); VERIFY_(STATUS)
         press3e = 0.0d0
         ALLOCATE(pl(i1:i2,j1:j2,1:km),             STAT=STATUS); VERIFY_(STATUS)

         pl(:,:,1:km)                 = (ple(:,:,0:km-1)+ple(:,:,1:km))*0.50
         tropopausePress(:,:)         = tropp(:,:)*Pa2hPa
         kzz(:,:,k1:k2-1)             = KH(:,:,km-1:1:-1)
         pctm1(i1:i2,j1:j2)           = ple(:,:,km)*Pa2hPa
         press3c(i1:i2,j1:j2,:)       = pl (:,:,km:1:-1)*Pa2hPa
         IF ( (TRIM(self%metType) == 'MERRA2') .OR. &
              (TRIM(self%metType) == 'FPIT')   .OR. &
              (TRIM(self%metType) == 'FP')   ) THEN 
            press3e(i1:i2,j1:j2,k1:k2)   = ple(:,:,km:1:-1)*Pa2hPa
         ELSEIF (TRIM(self%metType) == 'MERRA1') THEN 
            press3e(i1:i2,j1:j2,k1-1:k2) = ple(:,:,km:0:-1)*Pa2hPa
         END IF

         call updateDiffusion (cdt, self%vert_diffu_coef, pbl, &
                    tropopausePress, kzz, press3c, press3e, pctm1, &
                    concentration, self%isFixedConcentration, self%metType,     &
                    i1, i2, j1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)

         deallocate(kzz, pl)
         deallocate(pctm1)
         deallocate(press3e)
         deallocate(press3c)
         deallocate(tropopausePress)
      else if (self%diffu_opt == 2) then
         ALLOCATE(mass(i1:i2,j1:j2,1:km),STAT=STATUS)
         VERIFY_(STATUS)
         allocate(gridBoxHeight(i1:i2, j1:j2, k1:k2), STAT=STATUS)
         VERIFY_(STATUS)

         CALL MAPL_GetPointer(impDiff,   totalMass, 'MASS', RC=STATUS)
         VERIFY_(STATUS)
         CALL MAPL_GetPointer(impDiff,       zle,     'ZLE', RC=STATUS)
         VERIFY_(STATUS)

         gridBoxHeight(:,:,1:km) = zle(:,:,km-1:0:-1)-zle(:,:,km:1:-1)    ! m
         mass(:,:,1:km)          = totalMass(:,:,km:1:-1)                 ! kg

         call Update_PBL_Mixing (cdt, self%pbl_mixing_tau, pbl, &
                    gridBoxHeight, mass, concentration, &
                    i1, i2, j1, j2, k1, k2, numSpecies)

         deallocate(mass, gridBoxHeight)
      end if

      ! Pass the tracers back to the ESMF bundle
      !-----------------------------------------
      DO ic = 1, numSpecies
         call ESMF_FieldBundleGet(DiffTR, ic, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

         call ESMFL_BundleGetPointerToData(DiffTR, NAME, S, RC=STATUS)
         VERIFY_(STATUS)

         ! Do not forget to re-order the vertical levels
         S(:,:,:) = concentration(ic)%pArray3d(:,:,km:1:-1) 
      END DO

      CALL CleanArrayPointer(concentration, STATUS)
      VERIFY_(STATUS)

      deallocate(pbl)

      return

      end subroutine RunDiffusion
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FinalizeDiffusion
!
! !INTERFACE:
!
      subroutine FinalizeDiffusion (self)
!
! !INPUT/OUTPUT PARAMETERS:
      type (Diffusion_GridComp), intent(inout) :: self
!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR) :: IAm = "finalizeDiffusion"
!EOP
!-------------------------------------------------------------------------
!BOC
      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Completed Diffusion"
      ENDIF

      return

      end subroutine FinalizeDiffusion
!EOC
!-------------------------------------------------------------------------
  end module GmiDiffusionMethod_mod
