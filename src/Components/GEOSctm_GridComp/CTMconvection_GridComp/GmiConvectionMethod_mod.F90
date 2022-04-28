!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiConvectionMethod_mod 
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
      module GmiConvectionMethod_mod
!
! !USES:
      use ESMF
      use MAPL_Mod
      USE Chem_UtilMod
      use GmiArrayBundlePointer_mod
      USE GmiESMFrcFileReading_mod
      USE convectiveTransport_mod
!
      implicit none
!
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
!
! !PUBLIC MEMBER FUNCTIONS:
      private 
      public  :: initializeGmiConvection
      public  :: runGmiConvection
      public  :: finalizeGmiConvection
!
! !PUBLIC DATA MEMBERS:
      public  :: gmiConvection_GridComp
!
# include "setkin_par.h"
# include "setkin_depos.h"
# include "setkin_mw.h"
# include "setkin_lchem.h"

      TYPE gmiConvection_GridComp
        character(len=ESMF_MAXSTR) :: name = "GMI Convection"

        ! Dimension
        INTEGER :: i1, i2, im, j1, j2, jm, km

        ! Surface area of grid cells
        REAL(KIND=DBL), POINTER :: cellArea(:,:)

        integer :: convec_opt   ! Convection option
        logical :: det_ent      ! flag for doing detrainment then entrainment
        logical :: do_downdraft ! flag for doing downdrafts
        logical :: do_old_ncar  ! flag for using old ncar routine

        logical :: do_wetdep
        logical :: do_drydep
        logical :: pr_wet_depos
        integer :: met_opt
        integer :: chem_opt
        REAL(KIND=DBL), POINTER :: mw(:)
        logical, pointer :: isFixedConcentration(:)
        integer, pointer :: mapTracer(:)
        logical :: FIRST
      END TYPE gmiConvection_GridComp
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
! !IROUTINE: initializeGmiConvection
!
! !INTERFACE:
!
      subroutine initializeGmiConvection (self, impConv, expConv, nymd, nhms, &
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
      type (gmiConvection_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Convection related variables from the resource file.
!
! !LOCAL VARIABLES:
      ! Grid cell area can be set by initialize
      ! ---------------------------------------
      REAL, POINTER, DIMENSION(:,:) :: cellArea
      integer ::  i1, i2, j1, j2, km, im, jm
      integer :: k1, k2, ilo, ihi, julo, jhi
      integer                    :: numSpecies, ic, is
      integer                    :: STATUS
      type (ESMF_Config)         :: convConfigFile
      REAL, ALLOCATABLE :: var2D(:,:)
      character(len=ESMF_MAXSTR) :: varName, fileName
      character(len=ESMF_MAXSTR) :: IAm = "initializeGmiConvection"
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

      ! ------------------------------
      ! convec_opt
      !   0:  no convection
      !   1:  do DAO2 convection
      !   2:  do NCAR convection
      !   3:  do GMAO GEOS4 convection
      !------ ------------------------

      call ESMF_ConfigGetAttribute(convConfigFile, self%convec_opt, &
                      label   = "convec_opt:",&
                      default = 0, __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%det_ent, &
                     label = "det_ent:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%do_downdraft, &
                     label = "do_downdraft:", default=.false., __RC__  )

      call ESMF_ConfigGetAttribute(convConfigFile, self%do_old_ncar, &
                     label = "do_old_ncar:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%met_opt, &
                      label   = "met_opt:",&
                      default = 0, __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%chem_opt, &
                      label   = "chem_opt:",&
                      default = 0, __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%do_wetdep, &
                     label = "do_wetdep:", default=.false., __RC__ )

      call ESMF_ConfigGetAttribute(convConfigFile, self%do_drydep, &
                     label = "do_drydep:", default=.false., __RC__ )

      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT*," -----> det_ent      = ", self%det_ent
         PRINT*," -----> convec_opt   = ", self%convec_opt
         PRINT*," -----> do_downdraft = ", self%do_downdraft
         PRINT*," -----> do_old_ncar  = ", self%do_old_ncar
         PRINT *,"Done Reading the Convection Resource File"
      END IF

      !##############################
      ! End reading the resource file
      !##############################

     ! ---------------------------------------------------------------
     ! Check option ranges.  Note that as new options are added, these
     ! range checks will have to be modified.
     ! ---------------------------------------------------------------

      if ((self%convec_opt .LT. 0) .OR. (self%convec_opt .GT. 3)) then
         PRINT*,TRIM(IAm),': convec_opt out of range'
         STATUS = 1
         VERIFY_(STATUS)
      end if

      ! Grid box surface area, m^{2}
      ! ----------------------------
      CALL MAPL_GetPointer(impConv, cellArea, 'AREA', rc=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(self%cellArea(i1:i2,j1:j2), STAT=STATUS)
      VERIFY_(STATUS)
      self%cellArea(i1:i2,j1:j2)=cellArea(:,:)

      return

      end subroutine initializeGmiConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! IROUTINE: runGmiConvection
!
! !INTERFACE:
!
      subroutine runGmiConvection (self, impConv, expConv, nymd, nhms,  &
                    tdt, rc)
!
! !USES:
!
! !INPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: impConv ! Import State
      INTEGER, INTENT(IN) :: nymd, nhms          ! time
      REAL,    INTENT(IN) :: tdt                 ! chemical timestep (secs)
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(OUT) ::  rc                ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ESMF_State), INTENT(INOUT) :: expConv   ! Export State
      type (gmiConvection_GridComp), intent(inOut) :: self
!
! !DESCRIPTION:
! Runs the convection component.
!
! !LOCAL VARIABLES:
      REAL, POINTER, DIMENSION(:,:)   :: zpbl, lwi
      REAL, POINTER, DIMENSION(:,:,:) :: ple, zle, totalMass
      REAL, POINTER, DIMENSION(:,:,:) :: Q, CNV_MFC, CNV_MFD, T

      integer                    :: STATUS
      integer :: numSpecies, ic, kR, ik, is
      integer ::  i1, i2, j1, j2, km, im, jm
      integer :: k1, k2, ilo, ihi, julo, jhi
      REAL, ALLOCATABLE :: pl(:,:,:)

      integer       , ALLOCATABLE :: lwi_flags(:,:)
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
      REAL(KIND=DBL), allocatable :: humidity   (:, :, :)
      REAL(KIND=DBL), allocatable :: wet_depos  (:, :, :)
      REAL(KIND=DBL) :: cdt

      type (ESMF_Field)                   :: FIELD
      type (ESMF_Array)                   :: ARRAY
      type (ESMF_FieldBundle)             :: ConvTR
      REAL, POINTER, DIMENSION(:,:,:)     :: S

      type (ESMF_Config)         :: convConfigFile
      type (t_GmiArrayBundle), pointer :: concentration(:)
      character(len=ESMF_MAXSTR) :: IAm = "runGmiConvection"
      character(len=ESMF_MAXSTR) :: NAME, speciesName 
      integer :: ilong, ivert
      character (len=20) :: chem_mecha
      character (len=128) :: metdata_name_org
      character (len=128) :: metdata_name_model
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

         ALLOCATE(self%mw(numSpecies), STAT=STATUS)
         VERIFY_(STATUS)
         self%mw(:) = -999.00

         ALLOCATE(self%mapTracer(numSpecies), STAT=STATUS)
         VERIFY_(STATUS)
         self%mapTracer(:) = -999

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

!           IF ( MAPL_AM_I_ROOT() ) THEN
!             PRINT*,'Finding mw for ConvTR field <'//TRIM(NAME)//'>'
!           ENDIF
   
            !Identify fixed species such as O2, N2, ad.
            if (TRIM(NAME) == 'GMICHEM::ACET' .OR. TRIM(NAME) == 'GMICHEM::N2'   .OR. &
                TRIM(NAME) == 'GMICHEM::O2'   .OR. TRIM(NAME) == 'GMICHEM::NUMDENS')  THEN
               self%isFixedConcentration(ic) = .TRUE.
            end if

            IF (TRIM(NAME) .EQ. "GMICHEM::AOADAYS") THEN
               self%mw(ic) = 1.0
               self%mapTracer(ic) = NSP
            END IF

            isLOOP: DO is = 1, NSP
               speciesName = lchemvar(is)
               IF (TRIM(speciesName) ==     "O3") speciesName="OX"
               IF (TRIM(speciesName) ==  "CFCl3") speciesName="CFC11"
               IF (TRIM(speciesName) == "CF2Cl2") speciesName="CFC12"
   
               if (TRIM(NAME) == 'GMICHEM::'//TRIM(speciesName)) then
                  self%mw(ic)        = mw_data(is)
                  self%mapTracer(ic) = is
                  exit isLOOP
               end if
            END DO isLOOP
         END DO

         convConfigFile = ESMF_ConfigCreate(rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigLoadFile(convConfigFile, TRIM(rcfilen), rc=STATUS )
         VERIFY_(STATUS)

         ! This is only done if we are not doing GMI Chemistry
         call rcEsmfReadTable(convConfigFile, self%mw, "mw::", rc=STATUS)

         ASSERT_(COUNT(self%mw(:) < 0) == 0)

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
      CALL MAPL_GetPointer(impConv,       lwi,     'LWI', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,         Q,       'Q', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,         T,       'T', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,      zpbl,    'ZPBL', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,       ple,     'PLE', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv, totalMass,    'MASS', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,       zle,     'ZLE', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,   CNV_MFD, 'CNV_MFD', RC=STATUS); VERIFY_(STATUS)
      CALL MAPL_GetPointer(impConv,   CNV_MFC, 'CNV_MFC', RC=STATUS); VERIFY_(STATUS)

      allocate(lwi_flags(i1:i2,j1:j2),           STAT=STATUS); VERIFY_(STATUS)
      allocate(pbl(i1:i2,j1:j2),                 STAT=STATUS); VERIFY_(STATUS)
      allocate(press3c(ilo:ihi,julo:jhi,k1:k2),  STAT=STATUS); VERIFY_(STATUS)
      allocate(press3e(ilo:ihi,julo:jhi,k1-1:k2),STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(pl(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(mass(i1:i2,j1:j2,k1:k2),          STAT=STATUS); VERIFY_(STATUS)
      allocate(gridBoxHeight(i1:i2,j1:j2,k1:k2), STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(cmf(i1:i2,j1:j2,k1:k2),           STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(dtrain(i1:i2,j1:j2,k1:k2),        STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(eu(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(ed(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(md(i1:i2,j1:j2,k1:k2),            STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(kel(i1:i2,j1:j2,k1:k2),           STAT=STATUS); VERIFY_(STATUS)
      ALLOCATE(humidity(i1:i2,j1:j2,k1:k2),      STAT=STATUS); VERIFY_(STATUS)

      ! Increment land-water-ice flag to GMI
      ! specification: 1=water 2=land 3=ice
      ! ------------------------------------
      lwi_flags(:,:)     = 1+lwi(:,:)

      pbl     (:,:)       = zpbl(:,:)
      humidity(:,:,k1:k2) = Q        (:,:,km:1:-1)*ToGrPerKg
      kel     (:,:,k1:k2) = T        (:,:,km:1:-1)
      cmf     (:,:,k1:k2) = CNV_MFC  (:,:,km-1:0:-1)
      dtrain  (:,:,k1:k2) = CNV_MFD  (:,:,km:1:-1)
      mass    (:,:,k1:k2) = totalMass(:,:,km:1:-1)
      pl      (:,:,k1:k2) = (ple(:,:,0:km-1)+ple(:,:,1:km))*0.50

      press3c(i1:i2,j1:j2,:)       = pl (:,:,km:1:-1)*Pa2hPa
      press3e(i1:i2,j1:j2,k1-1:k2) = ple(:,:,km:0:-1)*Pa2hPa

      !eu(:,:,k2) = -cmf(:,:,k2) ! = 0.0
      !do ik = k1, k2-1
      !   eu(:,:,ik) = cmf(:,:,ik+1) - cmf(:,:,ik) + dtrain(:,:,ik)
      !enddo
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

      !#### wet_deps will have to be part of the Import / Export states ####
      ALLOCATE(wet_depos(i1:i2, j1:j2, numSpecies), STAT=STATUS)
      VERIFY_(STATUS)
      wet_depos(:,:,:) = 0.0
     
      cdt = tdt

      ilong = i2 - i1 + 1
      ivert = k2 - k1 + 1

      chem_mecha         = 'strat_trop'
      metdata_name_org   = 'GMAO'
      metdata_name_model = 'GEOS5'

      call updateGmiConvection (chem_mecha, metdata_name_org, metdata_name_model,&
                 self%det_ent, self%do_downdraft, self%do_wetdep, self%pr_wet_depos, &
                 self%chem_opt, IH2O2, IHNO3, lwi_flags, cdt, self%mw, self%mapTracer, &
                 self%cellArea, pbl, cmf, dtrain, eu, ed, md, gridBoxHeight, mass, &
                 wet_depos, kel, press3e, concentration, self%isFixedConcentration, &
#ifdef MICRO_AEROSOL
                 humidity, press3c, REL_SCAV_EFF_new, &
                 i1, i2, j1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, numSpecies)
#else
                 i1, i2, j1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, numSpecies)
#endif

      deallocate(pbl, lwi_flags)
      deallocate(cmf, dtrain, kel, humidity, eu, ed, md)
      deallocate(press3e, press3c, mass, gridBoxHeight, pl)

      ! Pass back the tracers to the ESMF Bundle
      !-----------------------------------------
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

      end subroutine runGmiConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FinalizeGmiConvection
!
! !INTERFACE:
!
      subroutine FinalizeGmiConvection (self)
!
! !INPUT/OUTPUT PARAMETERS:
      type (gmiConvection_GridComp), intent(inout) :: self
!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR) :: IAm = "finalizeGmiConvection"
!EOP
!-------------------------------------------------------------------------
!BOC
      IF ( MAPL_AM_I_ROOT() ) THEN
         PRINT *," "
         PRINT *,TRIM(IAm)//":"
         PRINT *,"Completed Convection"
      ENDIF

      return

      end subroutine FinalizeGmiConvection
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateGmiConvection
!
! !INTERFACE:
!
      subroutine updateGmiConvection (chem_mecha, metdata_name_org, metdata_name_model, &
                    det_ent, do_downdraft, do_wetdep, pr_wet_depos, chem_opt,      &
                    ih2o2_num, ihno3_num, lwi_flags, tdt, mw, mapTracer, mcor, pbl,&
                    cldmas, dtrn, eu, ed, md, grid_height, mass, wet_depos, kel,   &
                    press3e, concentration, isFixedConcentration, &
#ifdef MICRO_AEROSOL
                    humidity, press3c, REL_SCAV_EFF_new, &
                    i1, i2, ju1, j2, k1, k2, ilo, ihi,   &
                    julo, jhi, ilong, ivert, num_species)
#else
                    i1, i2, ju1, j2, k1, k2, ilo, ihi,   &
                    julo, jhi, ilong, ivert, num_species)
#endif
!
! !USES:
      use convectiveTransport_mod

#ifdef MICRO_AEROSOL
#     include "gmi_micro_aerosol.h"
#elif GOCARTaerosol
#     include "gocart_aerosol.h"
#elif strat_trop_aerosol
#     include "gocart_aerosol.h"
#elif strat_trop
#     include "gocart_aerosol.h"
!#else
!#     include "gmi_aerosol.h"
#endif

!c?   Tight coupling to setkin?
#     include "setkin_par.h"
#     include "setkin_depos.h"

#ifdef MICRO_AEROSOL
#     include "gmi_time_constants.h"
#     include "umaerosol.h"
#endif

!
! !INPUT PARAMETERS:
      character (len=*)   :: chem_mecha
      character (len=*)   :: metdata_name_org
      character (len=*)   :: metdata_name_model
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: ilong, ivert, num_species
      integer, intent(in) :: chem_opt                              ! Chemistry option
      integer, intent(in) :: ih2o2_num                             ! const array index for H2O2
      integer, intent(in) :: ihno3_num                             ! const array index for HNO3
      integer, intent(in) :: mapTracer(num_species)                ! 
      integer, intent(in) :: lwi_flags(i1:i2, ju1:j2)              ! land water ice flags (1-water 2-land 3-ice)
      logical, intent(in) :: isFixedConcentration(:)               ! determine ifa species is fixed or not
      logical, intent(in) :: det_ent                               ! flag for doing detrainment then entrainment
      logical, intent(in) :: do_downdraft                          ! flag for doing downdrafts
      logical, intent(in) :: do_wetdep                             ! flag for doing wet deposition
      logical, intent(in) :: pr_wet_depos                          ! flag for accumulated wet deposition
      real*8 , intent(in) :: tdt                                   ! model time step  (s)
      real*8 , intent(in) :: mw(num_species)                       ! array of species' molecular weights (g/mol)
      real*8 , intent(in) :: mcor       (i1:i2,  ju1:j2)           ! horizontal area of each grid box    (m^2)
      real*8 , intent(in) :: pbl        (i1:i2,  ju1:j2)           ! planetary boundary layer thickness  (m)
      real*8 , intent(in) :: cldmas     (i1:i2,  ju1:j2,  k1:k2)   ! convective mass flux in     updraft (kg/m^2/s)
      real*8 , intent(in) :: dtrn       (i1:i2,  ju1:j2,  k1:k2)   ! detrainment rate (DAO:kg/m^2*s, NCAR:s^-1)
      real*8 , intent(in) :: eu         (i1:i2,  ju1:j2,  k1:k2)   ! entrainment into convective updraft (s^-1)
      real*8 , intent(in) :: ed         (i1:i2,  ju1:j2,  k1:k2)   ! entrainment into convective downdraft (s^-1)
      real*8 , intent(in) :: md         (i1:i2,  ju1:j2,  k1:k2)   ! convective mass flux in downdraft (kg/m^2/s)
      real*8 , intent(in) :: grid_height(i1:i2,  ju1:j2,  k1:k2)   ! grid box height  (m)
      real*8 , intent(in) :: mass       (i1:i2,  ju1:j2,  k1:k2)   ! mass of air in each grid box (kg)
      real*8 , intent(in) :: kel        (ilo:ihi,julo:jhi,k1:k2)   ! temperature      (degK)
      real*8 , intent(in) :: press3e    (ilo:ihi,julo:jhi,k1-1:k2) ! pressure at the edge of each grid box (mb)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 ,                intent(inOut) :: wet_depos(i1:i2,ju1:j2,num_species) ! accumulated wet deposition (kg/m^2)
      type(t_GmiArrayBundle), intent(inout) :: concentration(num_species) ! species concentration, known at zone centers (mixing ratio)

#ifdef MICRO_AEROSOL
      real*8 , intent(in ) :: humidity   (i1:i2,  ju1:j2,  k1:k2)  ! specific humidity (g/kg)
      real*8 , intent(in ) :: press3c    (ilo:ihi,julo:jhi,k1:k2)  !  pressure at the center of each grid box (mb)
      real*8 , intent(out) :: REL_SCAV_EFF_new(i1:i2,ju1:j2,k1:k2,num_species) ! relative scavenging efficiency for aerosols (0-1)
#endif
!
! !DESCRIPTION:
! Interface to Convective Transport.
! It formats variables to satisfy the Convective Transport routine.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: MWTAIR =  28.96d0 ! Molecular weights (g/mol)
      real*8, parameter :: GMI_G  =   9.81d0 ! mean surface gravity accel. (m/s^2)
      real*8, parameter :: PASPMB = 100.00d0 ! pascals  per millibar

!     ----------------------------------------------------------------
!     REL_SCAV_EFF : scavenging efficiency of each aerosol relative to
!                    sulfate; all set to 1.0 for now (unitless)
!     ----------------------------------------------------------------

      integer, parameter :: NUM_AEROSOL = 37
      integer, parameter :: NUM_WAV_AER = 1

      real*8, parameter :: REL_SCAV_EFF(NUM_AEROSOL) =  &
                             (/ 1.0d0, 0.3d0, 1.0d0, 0.3d0, 1.0d0,  &
                                1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
                                1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
                                1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
                                1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  &
                                1.0d0, 0.3d0, 0.3d0, 0.3d0, 0.3d0,  &
                                0.3d0, 0.4d0, 0.4d0, 1.0d0, 1.0d0,  &
                                1.0d0, 0.3d0 /)

      real*8, parameter :: MBSTH = 1.0d-15  ! threshold below which we treat 
                                            ! mass fluxes as zero (mb/s)

#ifdef nonZeroInd
      integer, parameter :: IFSO2_l = 1
      integer, parameter :: INSO2_l = 1
#elif nonZeroInd_tracers
      integer, parameter :: IFSO2_l = 1
      integer, parameter :: INSO2_l = 1
#else
      integer, parameter :: IFSO2_l = IFSO2
      integer, parameter :: INSO2_l = INSO2
#endif

!EOP
!-----------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: iku
      integer :: il, ij, ik, ic
      integer :: il2g                      ! gathered index to operate over
#ifdef MICRO_AEROSOL
      integer :: inum
#endif

      integer :: itdt_conv
      integer :: num_conv_steps

      real*8  :: rnum_conv_steps
      real*8  :: tdt_conv

      real*8  :: xmbsth

      integer :: ideep(ilong)              ! gathering array
      integer :: pbli (ilong)              ! index of pbl height

      real*8  :: col_sum_depos   (i1:i2)   ! sum of column wet_deposition loss
      real*8  :: kloss           (i1:i2)   ! wet loss rate constant at all
                                           ! longitudes (s^-1)
      real*8  :: updraft_velocity(i1:i2)   ! velocity in convective updraft
                                           ! (m/s)

      real*8  :: dpi(ilong,  k1:k2)         ! delta pressure between interfaces
      real*8  :: dui(ilong,  k1:k2)         ! mass detraining from updraft
      real*8  :: eui(ilong,  k1:k2)         ! mass entraining into updraft
      real*8  :: mui(ilong,  k1:k2)         ! mass flux up
      real*8  :: mdi(ilong,  k1:k2)         ! mass flux down

      real*8  ::  &
     &  fracis(i1:i2, k1:k2, num_species)  ! insoluble fraction of tracer
      real*8  ::  &
     &  qq    (i1:i2, k1:k2, num_species)  ! tracer array including moisture

! mk
      real*8 :: nemui(ilong, ivert)     ! non-entraining mass flux up
      real*8 :: eneui(ilong, ivert)    ! mass entraining into non-entraining updraft
      real*8 :: dneui(ilong, ivert)    ! mass detraining from non-entraining updraft
      real*8 :: zero(ilong, ivert)     ! array of zeros for non-ent updraft calc

#ifdef MICRO_AEROSOL
!-micro_aerosol----begin------------------------------------------------------
!     --------------------
!     adding for umaerosol
!     --------------------

      real*8  :: relhume   (i1:i2, ju1:j2, k1:k2)
      real*8  :: h2osat    (i1:i2, ju1:j2, k1:k2)
      real*8  :: h2ogas    (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4mfrac  (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4dens   (i1:i2, ju1:j2, k1:k2)
      real*8  :: wetmas    (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4aer    (i1:i2, ju1:j2, k1:k2, naer)
      real*8  :: so4radv   (i1:i2, ju1:j2, k1:k2, nso4)
#ifndef NONON
      real*8  :: so4non    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: aernon    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: xnonum    (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: xso4non   (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: fso4non   (i1:i2, ju1:j2, k1:k2, nnon)
      real*8  :: rso4min   (i1:i2, ju1:j2, k1:k2)
#endif
      real*8  :: xcldnum   (i1:i2, ju1:j2, k1:k2)
      real*8  :: so4cldcoag(i1:i2, ju1:j2, k1:k2)

!-micro_aerosol----end--------------------------------------------------------
#endif

      ideep(:) = 0
      pbli (:) = 0

      kloss(:) = 0.0d0

      updraft_velocity(:) = 0.0d0

      dpi(:,:) = 0.0d0
      dui(:,:) = 0.0d0
      eui(:,:) = 0.0d0
      mui(:,:) = 0.0d0
      mdi(:,:) = 0.0d0

      fracis(:,:,:) = 0.0d0

      nemui(:,:) = 0.0d0
      eneui(:,:) = 0.0d0
      dneui(:,:) = 0.0d0
      zero(:,:) = 0.0d0


      xmbsth = MBSTH


!cc!? For now just hardwire a few values for the radon/lead problem
!c    (chem_opt = 1).  Setkin will eventually provide this information
!c    for all species when doing a full chemistry calculation.

      if (chem_opt == 1) then         ! Radon/Lead chemistry

        aerosol(1) = 0
        aerosol(2) = 1

        hstar_wet(1)   = 9.3d-3
        hstar_wet(2)   = 0.0d0

        delH_298_over_R_wet(1) = 2600.0d0
        delH_298_over_R_wet(2) =    0.0d0

        retention_eff  (1) =    0.0d0
        retention_eff  (2) =    0.0d0

      else if (chem_opt == 6) then     ! Beryllium chemistry

        aerosol(1) = 15
        aerosol(2) = 15

        hstar_wet(1)   = 0.0d0
        hstar_wet(2)   = 0.0d0

        delH_298_over_R_wet(1) =    0.0d0
        delH_298_over_R_wet(2) =    0.0d0

        retention_eff  (1) =    0.0d0
        retention_eff  (2) =    0.0d0

      else if (chem_opt == 8) then    ! Sulfur chemistry

        PRINT*,'Cannot run chem_opt 8 --  stopping'
        STOP
!       hstar_wet(IFSO2_l) = 600.0d0
!       hstar_wet(INSO2_l) = 600.0d0

      end if

       if (IFSO2 /= 0 .OR. INSO2 /= 0) then
         PRINT*,'Code is not able to run the FSO2,NSO2 case - stopping'
         STOP
!.sds. make like gocart
!      if(IFSO2 .gt. 0) hstar_wet(IFSO2_l) = 600.0d0
!      if(INSO2 .gt. 0) hstar_wet(INSO2_l) = 600.0d0
!.sds. 
       end if


!     -----------------------------------------------------------
!     Calculate any needed sub-cycling of the convection operator
!     by comparing the mass flux over the full time step to the
!     mass within the grid box.
!     -----------------------------------------------------------

      num_conv_steps =  &
     &  Maxval (tdt * cldmas(:,:,:) * Spread (mcor(:,:), 3, k2-k1+1) /  &
     &          mass(:,:,:)) + 1.0d0

      rnum_conv_steps = num_conv_steps
      tdt_conv        = tdt / rnum_conv_steps

#ifdef MICRO_AEROSOL
!-micro_aerosol----begin------------------------------------------------------
!     if(chem_mecha == 'micro_aerosol') then
!     --------------------
!     adding for umaerosol
!     --------------------

!     compute relative humidity (0-1)

      relhume(:,:,:) = 1.0d0 - (373.15d0 / kel(i1:i2,ju1:j2,k1:k2))
      relhume(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * relhume(:,:,:)    -  &
     &                    1.9760d0 * relhume(:,:,:)**2 -  &
     &                    0.6445d0 * relhume(:,:,:)**3 -  &
     &                    0.1299d0 * relhume(:,:,:)**4)

      relhume(:,:,:) =  &
     &  humidity(:,:,:) * MWTAIR / 18.0d0 /  &
     &  GPKG * press3c(i1:i2,ju1:j2,k1:k2) / relhume(:,:,:)

      relhume(:,:,:) = Max (Min (relhume(:,:,:), 0.95d0), 0.0d0)

!c    compute so4 aerosol size (m)

      so4aer(:,:,:,1) = concentration(ISO4M1)%pArray3D(:,:,:)
      so4aer(:,:,:,2) = concentration(ISO4N1)%pArray3D(:,:,:)
      so4aer(:,:,:,3) = concentration(ISO4M2)%pArray3D(:,:,:)
      so4aer(:,:,:,4) = concentration(ISO4N2)%pArray3D(:,:,:)

      do il = i1, i2
        do ij = ju1, j2
          do ik = k1, k2
            h2osat(il,ij,ik) = h2osat_f(kel(il,ij,ik))
            h2ogas(il,ij,ik) = relhume(il,ij,ik)*h2osat(il,ij,ik)

            so4mfrac(il,ij,ik) = so4mfrac_f(kel(il,ij,ik),  &
     &                                      h2osat(il,ij,ik),  &
     &                                      h2ogas(il,ij,ik))
            so4dens(il,ij,ik) = so4dens_f(kel(il,ij,ik),  &
     &                                    so4mfrac(il,ij,ik))
          end do
        end do
      end do

      do iso4 = 1, nso4

        iso4n = iso4 * nmomso4
        iso4m = iso4n - 1

        wetmas(:,:,:) = max(r2so4min, so4aer(:,:,:,iso4m)  &
     &                / max(epsilo,   so4aer(:,:,:,iso4n)))  &
     &                / so4mfrac(:,:,:)

        so4radv(:,:,:,iso4) = (r3q * wetmas(:,:,:)  &
     &                      / (GMI_PI * so4dens(:,:,:)))**r1td
        so4radv(:,:,:,iso4) = max(so4radvmin,  &
     &                        min(so4radv(:,:,:,iso4),so4radvmax))

      end do
#ifndef NONON
      do ic = 1, nnon
         ! from ISO4NOC to ISO4S4
         so4non(:,:,:,ic) = concentration(ISO4NOC+ic-1)%pArray3D(:,:,:)
         ! from INOC to ISSLT4
         aernon(:,:,:,ic) = concentration(INOC   +ic-1)%pArray3D(:,:,:)
      end do
!      so4non(:,:,:,1:nnon) = concentration(ISO4NOC:ISO4S4)%pArray3D(:,:,:)
!      aernon(:,:,:,1:nnon) = concentration(INOC:ISSLT4)%pArray3D(:,:,:)

!c    non-so4 aerosol radius (m) and mass (kg/particle)
      do inon = 1, nnon

        radvolm(:) = radgnon(:,inon)  &
     &             * exp(r3h*log(siggnon(:,inon))**2)
        radvnon(inon) = sum(fracnon(:,inon)*radvolm(:)**3)
        pmsnon(inon) = r4td*GMI_PI*rhonon(inon)  &
     &               * radvnon(inon)
        radvnon(inon) = radvnon(inon)**r1td

      end do

!c    non-so4 aerosol number concentration (#particles/kg air)
      do inon = 1, nnon
        xnonum(:,:,:,inon) = aernon(:,:,:,inon) / pmsnon(inon)
      end do

!c    so4 molecule # on non-so4 aerosol surface (#molec/particle)
      xso4non(:,:,:,1:nnon) = max(1.0d-30,so4non(:,:,:,1:nnon)) /  &
     &  max(r1,xnonum(:,:,:,1:nnon)) / so4min

!c    radius of each so4 molecule (m)
      rso4min(:,:,:) =  &
     &  (so4min * r3q /  &
     &  (GMI_PI*so4mfrac(:,:,:)*so4dens(:,:,:)))**r1td

!c    area fraction of so4 molecules on non-so4 aerosol surface
      do inon = 1, nnon

        fso4non(:,:,:,inon) =  &
     &    rso4min(:,:,:)**2.0d0 * xso4non(:,:,:,inon)  &
     &    / (4.0d0*(radvnon(inon)+rso4min(:,:,:))**2.0d0)

      end do
#endif

!c    compute cloud droplet number (cm^-3) to be the sum of
!c    mode 2 so4 number and non-so4 aerosol

      xcldnum(:,:,:) = so4aer(:,:,:,4)

#ifndef NONON
      do inon = 1, nnon

        if (inon <= 9) then                  ! 5 oc/bc, dust 1-4
          xcldnum(:,:,:) = xcldnum(:,:,:) + xnonum(:,:,:,inon)  &
     &                   * min(fso4non(:,:,:,inon)/fso4crit, r1)
        else                                 ! ss 1-4
          xcldnum(:,:,:) = xcldnum(:,:,:) + xnonum(:,:,:,inon)
        end if

      end do
#endif

!c    change unit from kg^-1 to cm^-3
      xcldnum(:,:,:) = min(r3000,max(r10,xcldnum(:,:,:) *  &
     &  (press3c(i1:i2,ju1:j2,k1:k2) * MWTAIR * MB2CGS) /  &
     &  (GPKG * kel(i1:i2,ju1:j2,k1:k2) * BOLTZMN_E * AVOGAD)))

!c    compute Brownian coagulation coefficient K12 (cm3 s-1)
!c    for the 1st so4 mode with cloud droplets
!c    fit from Table 12.3 and Fig.12.5, in Seinfeld and Pandis (1997)
!c    as log10(K12) = -1.8436 * log10(D(um)) - 9.21

      so4radv(:,:,:,1) = min(so4radv(:,:,:,1),1.0d-7)

      so4cldcoag(:,:,:) =  &
     &  10 ** (a11 * log10(so4radv(:,:,:,1)*2.0d+6) + a12)

!     ==============================
      ICLOOP: do ic = 1, num_species
!     ==============================

!c      calculate scavenging efficiency of aerosols

        if (isFixedConcentration(ic) .or.  &
     &     (aerosol(mapTracer(ic)) <= 0)  .or.  &
     &     (aerosol(mapTracer(ic)) >  NUM_AEROSOL)) then
!         ============
          cycle ICLOOP
!         ============
        end if

!c      SO4 mode 1 (assume cloud lifetime 4 hr)
        if (aerosol(mapTracer(ic)) == 1) then

          REL_SCAV_EFF_new(:,:,:,mapTracer(ic)) = min(r1,  &
     &      so4cldcoag(:,:,:) * 4.0d0 * SECPHR * xcldnum(:,:,:))

        else
#ifdef NONON
          REL_SCAV_EFF_new(:,:,:,mapTracer(ic)) = REL_SCAV_EFF(aerosol(mapTracer(ic)))
#else
          inon = aerosol(mapTracer(ic)) - nso4

          if (inon>0 .and. inon<=9) then         ! 5 oc/bc, dust1-4

            REL_SCAV_EFF_new(:,:,:,mapTracer(ic)) =  &
     &        min(fso4non(:,:,:,inon)/fso4crit, r1)

          else                             ! mode 2 so4, ss 1-4

            REL_SCAV_EFF_new(:,:,:,mapTracer(ic)) = REL_SCAV_EFF(aerosol(mapTracer(ic)))

          end if
#endif
        end if

!     =============
      end do ICLOOP
!     =============
!-micro_aerosol----end--------------------------------------------------------
!     endif     ! chem_mecha == 'micro_aerosol'
#endif

!     =======
      ijloop: do ij = ju1, j2
!     =======

        il2g = 0

!       =======
        illoop: do il = i1, i2
!       =======

!         ----------------------------------------------------
!         Verify that there is convection in the current cell.
!         ----------------------------------------------------

          if (Maxval (cldmas(il,ij,:)) >= 0.0d0) then

            il2g = il2g + 1

            ideep(il2g) = il

            dpi(il2g,:) =  &
     &        (press3e(il,ij,k2-1:k1-1:-1) -  &
     &         press3e(il,ij,k2:k1:-1)) *  &
     &        PASPMB


            mui(il2g,:) = cldmas(il,ij,k2:k1:-1) * GMI_G

            if (metdata_name_org(1:3) == 'DAO' .or. &
     &          (metdata_name_org(1:4) == 'GMAO' .and.  &
     &           metdata_name_model(1:5) == 'GEOS5')  ) then

              eui(il2g,k1+1:k2-1) = Max (0.0d0,  &
     &               cldmas(il,ij,k2-1:k1+1:-1) -  &
     &               cldmas(il,ij,k2-2:k1  :-1) +  &
     &               dtrn  (il,ij,k2-1:k1+1:-1))

              eui(il2g,:) = eui(il2g,:) * GMI_G

              if (det_ent) dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * GMI_G

              if (metdata_name_model(1:2) /= 'GS') then

                if(do_downdraft) mdi(il2g,:) = md(il,ij,k2:k1:-1)*GMI_G

              end if

            else

              eui(il2g,:) = eu    (il,ij,k2:k1:-1) * dpi(il2g,:)
!

              if (det_ent)  &
     &          dui(il2g,:) = dtrn(il,ij,k2:k1:-1) * dpi(il2g,:)

              if(do_downdraft) then

!               --------------------------------------------------------------
!               Special algorithm for calculating downdraft mass flux from the
!               GISS downdraft entrainment. Limited to grid boxes 1-7 and
!               limited to be less than or equal to updraft mass flux.
!               --------------------------------------------------------------

                if (metdata_name_org(1:4) == 'GISS') then
!mk
                  mdi(il2g,1:k2-7) = 0.0
!mk jan 21
!                 mdi(il2g,k2-6:k2) = ed(il,ij,1:7) * dpi(il2g,k2-6:k2)
!mk
                  mdi(il2g,k2-6:k2) = ed(il,ij,7:1:-1) * dpi(il2g,k2-6:k2)

                  do ik = k2-6, k2-1

                    mdi(il2g,ik) = mdi(il2g,ik-1) - 2.0d0 * mdi(il2g,ik)

                  end do

                  mdi(il2g,k2) = 0.0d0

                  where (abs(mdi(il2g,:)) > mui(il2g,:))  &
     &              mdi(il2g,:) = -mui(il2g,:)

                else

                  mdi(il2g,:) = md(il,ij,k2:k1:-1) * GMI_G

                end if

              end if

            end if

          end if

!       =============
        end do illoop
!       =============

!mk
!... parameterize deep convection in GISS fields
!
!... need to parameterize "deep" (non-entraining) convection from
!...  the existing GISS fields. 'eui' has been converted to top-down
!...  so need to work from level 5 from the bottom upward
!...  The method was suggested by UCI (Prather)
!
        if(metdata_name_org(1:4) == 'GISS' .and. il2g /= 0) then
!... start with 12% of value at level 5 of met field (empirical)
          nemui(1:il2g,k2-4) = 0.12d0*mui(1:il2g,k2-4)
!... build upward without entrainment
          do ik=k2-5,k1,-1
            nemui(1:il2g,ik) = Min(nemui(1:il2g,ik+1),mui(1:il2g,ik))
          enddo
!... build downward without detrainment
          do ik=k2-3,k2
            nemui(1:il2g,ik) = Min(nemui(1:il2g,ik-1),mui(1:il2g,ik))
          enddo

!... now remove this flux from the original total mass flux up field
!...  => mui now becomes the "shallow" mass flux up
          mui(1:il2g,:) = mui(1:il2g,:) - nemui(1:il2g,:)

!... calculate the entrainment into the "deep" updraft
          eneui(1:il2g,:) =  Max (0.0d0  &
     &      ,(nemui(1:il2g,1:k2-1) - nemui(1:il2g,2:k2)) )

!... redo "shallow" entrainment so total maintained
!...  => eui now becomes the "shallow" entrainment into updraft
          eui(1:il2g,:) =  eui(1:il2g,:) - eneui(1:il2g,:)

!... calculate the detrainment from the deep updraft
!mk
!         if (det_ent) then
!           dneui(1:il2g,:) =  Max (0.0d0
!    &       ,(nemui(1:il2g,2:k2) - nemui(1:il2g,1:k2-1)) )
!... redo "shallow" detrainment so total maintained
!...  => dui now becomes the "shallow" detrainment from updraft
!           dui(1:il2g,:) =  dui(1:il2g,:) - dneui(1:il2g,:)
!         endif

        endif


!       =======
        il2gif: if (il2g /= 0) then
!       =======

          if (do_wetdep) then

!           ----------------------------------------------------------------
!           Calculate the insoluble fraction at each level for each species.
!           ----------------------------------------------------------------

            where ((lwi_flags(i1:i2,ij) == 1) .or.  &
     &             (lwi_flags(i1:i2,ij) == 3))
              updraft_velocity(:) =  5.0d0
            elsewhere
              updraft_velocity(:) = 10.0d0
            end where

            do ik = k1, k2

              iku = k2 - ik + 1

              do ic = 1, num_species

#ifdef MICRO_AEROSOL
                if (mapTracer(ic) == ihno3_num .or. mapTracer(ic) == ISO4G) then
#else
                if (mapTracer(ic) == ihno3_num) then
#endif

                  kloss(:) = 5.0d-3

                else if (aerosol(mapTracer(ic)) >= 1) then


#ifdef MICRO_AEROSOL
                  kloss(:) = 5.0d-3  &
     &                     * REL_SCAV_EFF_new(:,ij,ik,mapTracer(ic))
#else
                  kloss(:) = 5.0d-3 * REL_SCAV_EFF(aerosol(mapTracer(ic)))
#endif

                else if (hstar_wet(mapTracer(ic)) > 0.0d0) then

!                 =======================
                  call Calc_Wet_Loss_Rate  &
!                 =======================
     &              (.true., mapTracer(ic), ih2o2_num, delH_298_over_R_wet(mapTracer(ic)),  &
     &               hstar_wet(mapTracer(ic)), retention_eff(mapTracer(ic)), press3e(:,ij,ik),  &
     &               kel(:,ij,ik), kloss(:), i1, i2, ilo, ihi)


                  kloss(:) = 5.0d-3 * kloss(:)

                else

                  kloss(:) = 0.0d0

                end if

                fracis(:,iku,ic) =  &
     &            Exp (-kloss(:) * grid_height(:,ij,ik) /  &
     &                 updraft_velocity(:))

              end do

            end do

          else

            fracis(:,:,:) = 1.0d0

          end if

!         ----------------------------------------------------------------
!         Find the index of the top of the planetary boundary layer (pbl).
!         Convection will assume well mixed tracers below that level.
!         ----------------------------------------------------------------

          do il = 1, il2g

            pbli(il) = 0

            ikloop: do ik = k1, k2

              if (pbl(ideep(il),ij) < Sum (grid_height(ideep(il),ij,k1:ik))) then
                pbli(il) = ik
                exit ikloop
              end if

            end do ikloop

            if (pbli(il) == 0) then

              err_msg = 'Could not find pbl in Do_Convec_Ncar.'
              PRINT*, err_msg
              stop
              !call GmiPrintError (err_msg, .true., 2, ideep(il), ij, 1, pbl(ideep(il),ij), 0.0d0)

            end if

            pbli(il) = k2 - pbli(il)

          end do



!            write(*,*)'hi mk : num_conv_steps = ',num_conv_steps
!         ==========
          itdtcloop: do itdt_conv = 1, num_conv_steps
!         ==========

            do ic = 1, num_species
            qq(:,k2:k1:-1,ic) = concentration(ic)%pArray3D(:,ij,k1:k2)
            end do

!mk
!... call the "deep" convection first
            if(metdata_name_org(1:4) == 'GISS' ) then
!              ==============
               call convectiveTransport  &
!              ==============
     &          (il2g, tdt_conv, xmbsth, ideep, pbli, dneui, eneui,  &
     &           nemui, zero, dpi, fracis, qq, isFixedConcentration, &
     &           i1, i2, k1, k2, ilong, num_species)
            endif

!... call the "shallow" convection for GISS or total for other met fields
!           ==============
            call convectiveTransport  &
!           ==============
     &        (il2g, tdt_conv, xmbsth, ideep, pbli, dui, eui,  &
     &         mui, mdi, dpi, fracis, qq, isFixedConcentration, &
     &           i1, i2, k1, k2, ilong, num_species)


            if (pr_wet_depos) then

!             ---------------------------------------------------
!             Calculate the wet deposition from the difference in
!             mixing ratio.
!             ---------------------------------------------------

              do ic = 1, num_species

                if ((hstar_wet(mapTracer(ic)) > 0.0d0) .or.  &
     &              (aerosol(mapTracer(ic)) >= 1)  .or.  &
     &              (mapTracer(ic) == ihno3_num)) then

                  col_sum_depos(i1:i2) = 0.0d0

                  do ik = k1, k2

                    iku = k2 - ik + 1

                    col_sum_depos(:) =  &
     &                col_sum_depos(:) +  &
     &                (concentration(ic)%pArray3D(:,ij,ik) - qq(:,iku,ic)) *  &
     &                mass(:,ij,ik)

                  end do


                  where (col_sum_depos(i1:i2) > 0.0d0)  &
     &              wet_depos(i1:i2,ij,ic) =  &
     &                wet_depos(i1:i2,ij,ic) +  &
     &                col_sum_depos(i1:i2) * mw(ic) /  &
     &                MWTAIR / mcor(i1:i2,ij)

                end if

              end do

            end if

            do ic = 1, num_species
              concentration(ic)%pArray3D(:,ij,k1:k2) = qq(:,k2:k1:-1,ic)
            end do

!         ================
          end do itdtcloop
!         ================

!       =============
        end if il2gif
!       =============

!     =============
      end do ijloop
!     =============

      return

      end subroutine updateGmiConvection
 end module GmiConvectionMethod_mod
