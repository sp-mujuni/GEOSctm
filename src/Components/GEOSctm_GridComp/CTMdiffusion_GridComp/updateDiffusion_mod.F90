!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: updateDiffusion_mod
!
! !INTERFACE:
!
#include "MAPL_Generic.h"
!
module updateDiffusion_mod
!
! !USES:
      use ESMF
      use MAPL_Mod
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      !use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: updateDiffusion
      public  :: Update_PBL_Mixing
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
! !IROUTINE: updateDiffusion
!
      subroutine updateDiffusion (tdt, vert_diffu_coef, pbl, tropp, kzz,&
                       press3c, press3e, pctm1, concentration, &
                       isFixedConcentration, metType, &
                       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)
!
! !INPUT PARAMETERS:

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: numSpecies
      real*8 , intent(in) :: tdt                               ! model time step (s)
      real*8 , intent(in) :: vert_diffu_coef                   ! scalar vertical diffusion coefficient (m^2/s)
      real*8 , intent(in) :: pbl(i1:i2, ju1:j2)                ! planetary boundary layer depth (m
      real*8 , intent(in) :: tropp(i1:i2, ju1:j2)              ! tropopause pressure (mb)
      real*8 , intent(in) :: press3c(ilo:ihi,julo:jhi,k1:k2)   ! atmospheric pressure at the center 
                                                               ! of each grid box (mb)
      real*8 , intent(in) :: press3e(ilo:ihi,julo:jhi,k1-1:k2) ! atmospheric pressure at the edge
                                                               ! of each grid box (mb)
      real*8 , intent(in) :: pctm1(ilo:ihi, julo:jhi)          ! surface pressure field at t1, 
                                                               ! known at zone centers (mb)
      logical, intent(in) :: isFixedConcentration(numSpecies)
      character(len=ESMF_MAXSTR), intent(in) :: metType
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: kzz(i1:i2, ju1:j2, k1:k2)         ! array of vertical diffusion coefficients (m^2/s)
                             ! species concentration (mixing ratio)
      type (t_GmiArrayBundle), intent(inout) :: concentration(numSpecies)
!
! !DESCRIPTION:
! Updates diffusion.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      call Adjust_Kzz (press3c, tropp, kzz, i1, i2, ju1, j2, k1, k2, &
                       ilo, ihi, julo, jhi)

      call Do_Vert_Diffu(tdt, vert_diffu_coef, pbl, kzz, press3c, &
                         press3e, pctm1, concentration, &
                         isFixedConcentration, metType, i1, i2, ju1, j2,&
                         k1, k2, ilo, ihi, julo, jhi, numSpecies)

      return

      end subroutine updateDiffusion 
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_PBL_Mixing
!
! !INTERFACE: 
!
      subroutine Update_PBL_Mixing (tdt, pbl_mixing_tau, pbl, grid_height, &
                                    mass, concentration, i1, i2, ju1, j2,  &
                                    k1, k2, numSpecies)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: numSpecies
      real*8 , intent(in) :: tdt                             ! model time step (s)
      real*8 , intent(in) :: pbl_mixing_tau                  ! time scale for pbl mixing (s)
      real*8 , intent(in) :: pbl        (i1:i2,ju1:j2)       ! planetary boundary layer depth (m)
      real*8 , intent(in) :: grid_height(i1:i2,ju1:j2,k1:k2) ! height of each grid box (m)
      real*8 , intent(in) :: mass       (i1:i2,ju1:j2,k1:k2) ! air mass in each grid box (kg)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inout) :: concentration(numSpecies)
!
! !DESCRIPTION:
!   Updates diffusion through a boundary layer mixing model.
!
! !LOCAL VARIABLES:
      integer :: ic, il, ij, ik, kpbl
      real*8  :: frac_pbl, height, mass_pbl, mixing_factor, well_mixed

!EOP
!------------------------------------------------------------------------------
!BOC
      mixing_factor = Exp(-tdt / pbl_mixing_tau)

      !---------------------------------------
      ! Loop over all longitudes and latitudes.
      !---------------------------------------

      do il = i1, i2
         do ij = ju1, j2
            mass_pbl = 0.0d0
            height = 0.0d0

            !-----------------------------------
            ! Find the total mass within the pbl.
            !-----------------------------------

            kloop: do ik = k1, k2
               height = height + grid_height(il,ij,ik)
               if (pbl(il, ij) >= height) then
                  mass_pbl = mass_pbl + mass(il, ij, ik)
               else
                  !-------------------------------------------------------
                  ! This grid box is fractionally in the PBL. Save part of
                  ! its mass also, save the k index, and exit the loop.
                  !-------------------------------------------------------

                  frac_pbl = (pbl(il,ij) -  &
                             (height - grid_height(il,ij,ik))) / grid_height(il,ij,ik)

                  mass_pbl = mass_pbl + mass(il,ij,ik) * frac_pbl

                  kpbl = ik

                  exit kloop
               end if
            end do kloop

            !--------------------------------------------------------
            ! Now if a real PBL was found (kpbl > k1) loop over the
            ! species, calculate a well mixed
            ! mixing ratio, and apply it to the grids within the pbl
            ! and then to the grid box
            ! which is fractionally within the pbl.
            !---------------------------------------------------------

            if (kpbl > k1) then
               do ic = 1, numSpecies
                  well_mixed = &
                     (Sum(concentration(ic)%pArray3D(il,ij,k1:kpbl-1) *  &
                          mass(il,ij,k1:kpbl-1)) +  &
                     frac_pbl * concentration(ic)%pArray3D(il,ij,kpbl) * &
                     mass(il,ij,kpbl)) / mass_pbl

                  concentration(ic)%pArray3D(il,ij,1:kpbl-1) =  &
                        concentration(ic)%pArray3D(il,ij,k1:kpbl-1) *  &
                        mixing_factor + well_mixed * (1.0d0 - mixing_factor)

                  concentration(ic)%pArray3D(il,ij,kpbl) = &
                        concentration(ic)%pArray3D(il,ij,kpbl) + frac_pbl *  &
                       (concentration(ic)%pArray3D(il,ij,kpbl)*(mixing_factor - 1.0d0) +  &
                        well_mixed * (1.0d0 - mixing_factor))
               end do
            end if
         end do
      end do

      return

      end subroutine Update_PBL_Mixing
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust_Kzz
!
! !INTERFACE:
!
      subroutine Adjust_Kzz(press3c, tropp, kzz, i1, i2, ju1, j2, &
                            k1, k2, ilo, ihi, julo, jhi)
!
! !INPUT PARAMETERS
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in) :: press3c(ilo:ihi,julo:jhi,k1:k2) ! atmospheric pressure at the center 
                                                             ! of each grid box (mb)
      real*8 , intent(in) :: tropp(i1:i2,ju1:j2)             ! tropopause pressure (mb)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inout) :: kzz(i1:i2,ju1:j2,k1:k2)      ! vertical diffusion coefficients (m^2/s)

!
! !DESCRIPTION:
!   Makes the kzz (vertical diffusion) zero above the tropopause.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (Maxval (kzz) > 0.0d0) then
         where (press3c(i1:i2,ju1:j2,:) < (Spread (tropp(:,:), 3, k2) + 100.0d0))
               kzz(i1:i2,ju1:j2,:) = 0.0d0
         end where
      end if

      return

      end subroutine Adjust_Kzz
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !IROUTINE: Do_Vert_Diffu
!
! !INTERFACE:
!
      subroutine Do_Vert_Diffu(tdt, vert_diffu_coef, pbl, kzz, press3c,&
                               press3e, psf, concentration, &
                               isFixedConcentration, metType, &
                               i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: numSpecies
      real*8 , intent(in) :: tdt                                 ! model time step (s)
      real*8 , intent(in) :: vert_diffu_coef                     ! scalar vertical diffusion coefficient (m^2/s)
      real*8 , intent(in) :: pbl(i1:i2, ju1:j2)                  ! planetary boundary layer depth (m)
      real*8 , intent(in) :: kzz(i1:i2, ju1:j2, k1:k2)           ! vertical diffusion coefficients (m^2/s)
      real*8 , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)   ! atmospheric pressure at the center of each grid box (mb)
      real*8 , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2) ! atmospheric pressure at the edge   of each grid box (mb)
      real*8 , intent(in) :: psf(ilo:ihi, julo:jhi)              ! surface pressure field known at zone centers (mb)
      logical, intent(in) :: isFixedConcentration(numSpecies)
      character(len=ESMF_MAXSTR), intent(in) :: metType
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inout) :: concentration(numSpecies)
!
! !DECLARED VARIABLES:
      real*8, parameter :: HZERO  = 8000.0d0  ! atmospheric scale height (m)
      real*8, parameter :: HZERO2 = HZERO * HZERO  ! (m^2)
!
! !DESCRIPTION:
!   Calculates vertical diffusion using an implicit technique.
!   Note: Constituent mass is conserved during mixing.
!          Constant field should remain constant.
!
! !LOCAL VARIABLES:
      integer :: ik, ic
      real*8  :: avdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: bvdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: cvdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: vdiff (i1:i2, ju1:j2, k1:k2-1)
      real*8  :: bb    (i1:i2, ju1:j2, k1:k2-1)
      real*8  :: cc    (i1:i2, ju1:j2, k1:k2-1)
      real*8  :: rr    (i1:i2, ju1:j2, k1:k2-1)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      avdiff(:,:,:) = 0.0d0; bvdiff(:,:,:) = 0.0d0
      cvdiff(:,:,:) = 0.0d0; vdiff (:,:,:) = 0.0d0

      bb(:,:,:) = 0.0d0; cc(:,:,:) = 0.0d0
      rr(:,:,:) = 0.0d0

      ! --------------------------------------------------------
      ! If there was no kzz in the met file, it was set to < 0.
      ! In this case calculate a vertical diffusion coefficient,
      ! otherwise use kzz.
      ! --------------------------------------------------------

      if (Maxval (kzz) < 0.0d0) then

         ! ------------------------------------------------------
         ! Set up the vertical diffusion coefficients.
         ! First put height (m) into vdiff, then put simple
         ! parameterization for diffusion coefficient into vdiff,
         ! then put it into units of per mb
         ! ------------------------------------------------------

         vdiff = 0.0d0

         do ik = 1, 10

            vdiff(:,:,ik) = -HZERO *  &
                            Log (press3e(i1:i2,ju1:j2,ik) / psf(i1:i2,ju1:j2))

            where (vdiff(:,:,ik) < pbl(:,:))
               vdiff(:,:,ik) = vert_diffu_coef * (1.0d0 -  &
                               2.0d0 * (Abs (vdiff(:,:,ik) - (pbl(:,:) * 0.5d0)) / pbl(:,:)))
            elsewhere
               vdiff(:,:,ik) = 0.0d0
            end where

            vdiff(:,:,ik) = vdiff(:,:,ik) / HZERO2 *  &
                            (press3e(i1:i2,ju1:j2,ik) * press3e(i1:i2,ju1:j2,ik)) /  &
                            (press3c(i1:i2,ju1:j2,ik) - press3c(i1:i2,ju1:j2,ik+1))
         end do
      else
         vdiff(:,:,k1:k2-1) = kzz(:,:,k1:k2-1) / HZERO2 *  &
                             (press3e(i1:i2,ju1:j2,k1:k2-1) *  &
                              press3e(i1:i2,ju1:j2,k1:k2-1)) /  &
                             (press3c(i1:i2,ju1:j2,k1:k2-1) -  &
                              press3c(i1:i2,ju1:j2,k1+1:k2))
      end if

      ! -----------------------------------------------------------
      ! Construct the tri-diagonal matrix which will be used in the
      ! implicit vertical diffusion solution.
      !   avdiff is the lower band
      !   bvdiff is the main diagonal
      !   cvdiff is the upper band
      ! -----------------------------------------------------------

      ! -----------------------
      ! First do the top layer.
      ! -----------------------

      avdiff(:,:,k2-1) = 0.0d0
      bvdiff(:,:,k2-1) =  vdiff(:,:,k2-2)
      cvdiff(:,:,k2-1) = -vdiff(:,:,k2-2)


      ! ----------------------
      ! Now do all the layers.
      ! ----------------------

      avdiff(:,:,k1+1:k2-2) = -vdiff(:,:,k1+1:k2-2)
      bvdiff(:,:,k1+1:k2-2) = (vdiff(:,:,k1+1:k2-2) + vdiff(:,:,k1  :k2-3))
      cvdiff(:,:,k1+1:k2-2) = -vdiff(:,:,k1  :k2-3)

      ! ------------------------------------
      ! Now do the layer next to the ground.
      ! ------------------------------------

      avdiff(:,:,k1) = -vdiff(:,:,k1)
      bvdiff(:,:,k1) =  vdiff(:,:,k1)
      cvdiff(:,:,k1) = 0.0d0

      ! ---------------------------------------------------
      ! Now convert them to the proper dimensionless units.
      ! ---------------------------------------------------

      IF ( (TRIM(metType) == 'MERRA2') .OR.  &
           (TRIM(metType) == 'FPIT')   .OR.  &
           (TRIM(metType) == 'FP')    ) THEN
         avdiff(i1:i2,ju1:j2,k1:k2-1) = avdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
               (press3e(i1:i2,ju1:j2,k1:k2-1) - press3e(i1:i2,ju1:j2,k1+1:k2))

         bvdiff(i1:i2,ju1:j2,k1:k2-1) = bvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
               (press3e(i1:i2,ju1:j2,k1:k2-1) - press3e(i1:i2,ju1:j2,k1+1:k2)) + 1.0d0

         cvdiff(i1:i2,ju1:j2,k1:k2-1) = cvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt  /  &
               (press3e(i1:i2,ju1:j2,k1:k2-1) - press3e(i1:i2,ju1:j2,k1+1:k2))
      ELSEIF (TRIM(metType) == 'MERRA1') THEN
         avdiff(i1:i2,ju1:j2,k1:k2-1) = avdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
               (press3e(i1:i2,ju1:j2,k1-1:k2-2) - press3e(i1:i2,ju1:j2,k1:k2-1))

         bvdiff(i1:i2,ju1:j2,k1:k2-1) = bvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
               (press3e(i1:i2,ju1:j2,k1-1:k2-2) - press3e(i1:i2,ju1:j2,k1:k2-1)) + 1.0d0

         cvdiff(i1:i2,ju1:j2,k1:k2-1) = cvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt  /  &
               (press3e(i1:i2,ju1:j2,k1-1:k2-2) - press3e(i1:i2,ju1:j2,k1:k2-1))
      END IF

      ! ---------------------
      ! Now solve the system.
      ! ---------------------

      icloop: do ic = 1, numSpecies

         !   ======================================
         if (isFixedConcentration(ic)) cycle icloop
         !   ======================================

         ! -----------------------------------------
         ! Set up the right hand side of the system.
         ! -----------------------------------------

         bb(:,:,k2-1) = bvdiff(:,:,k2-1)
         rr(:,:,k2-1) = concentration(ic)%pArray3D(:,:,k2-1)

         ! ---------------------------------
         ! Eliminate the lower diagonal (a).
         ! ---------------------------------

         do ik = k2 - 1, k1 + 1, -1
            cc(:,:,ik) = cvdiff(:,:,ik) / bb(:,:,ik)
            rr(:,:,ik) =     rr(:,:,ik) / bb(:,:,ik)

            bb(:,:,ik-1) = bvdiff(:,:,ik-1) - avdiff(:,:,ik-1) * cc(:,:,ik)
            rr(:,:,ik-1) = concentration(ic)%pArray3D(:,:,ik-1) -  &
                           avdiff(:,:,ik-1) * rr(:,:,ik)
         end do

         ! -------------------------------
         ! Solve for the new mixing ratio.
         ! -------------------------------

         concentration(ic)%pArray3D(:,:,k1) = rr(:,:,k1) / bb(:,:,k1)

         do ik = k1 + 1, k2 - 1
            concentration(ic)%pArray3D(:,:,ik) = rr(:,:,ik) -  &
                      concentration(ic)%pArray3D(:,:,ik-1) * cc(:,:,ik)
         end do
      end do icloop

      return

      end subroutine Do_Vert_Diffu
!EOC
!------------------------------------------------------------------------------

end module updateDiffusion_mod
