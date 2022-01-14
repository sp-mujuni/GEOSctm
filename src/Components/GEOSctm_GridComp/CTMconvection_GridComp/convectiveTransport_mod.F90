!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: convectiveTransport_mod 
!
! !INTERFACE:
!
module convectiveTransport_mod

      USE GmiArrayBundlePointer_mod

      implicit none

      private 
      public  :: doConvectiveTransport
!
! !DESCRIPTION:
! Module for the convective transport.
!
!EOP
!-----------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doConvectiveTransport
!
! !INTERFACE:
!
      subroutine doConvectiveTransport (det_ent, do_downdraft, pbl, cldmas, &
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
                  err_msg = 'Could not find pbl in doConvectiveTransport.'
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

      end subroutine doConvectiveTransport
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convectiveTransport
!
! !INTERFACE:
!
      subroutine convectiveTransport (il2g, delt, xmbsth, ideep, pbli, &
                           dui, eui, mui, mdi, dpi, fracis, qq,        &
                           isFixedConcentration, i1, i2, k1, k2, ilong, numSpecies)
!
! !INPUT PARAMETERS
      logical, intent(in) :: isFixedConcentration(:)
      integer, intent(in) :: i1, i2, k1, k2, ilong, numSpecies
      integer, intent(in) :: il2g    ! gathered max lon indices over which to operate
      real*8,  intent(in) :: delt    ! convection time step (s)
      real*8,  intent(in) :: xmbsth  ! threshold below which we treat mass fluxes as zero (mb/s)
      integer, intent(in) :: ideep (ilong) ! gathering array
      integer, intent(in) :: pbli  (ilong) ! index of planetary boundary layer
      real*8,  intent(in) :: dui   (ilong, k1:k2) ! mass detraining from updraft
      real*8,  intent(in) :: eui   (ilong, k1:k2) ! mass entraining into updraft
      real*8,  intent(in) :: mui   (ilong, k1:k2) ! mass flux up
      real*8,  intent(in) :: mdi   (ilong, k1:k2) ! mass flux down
      real*8,  intent(in) :: dpi   (ilong, k1:k2) ! delta pressure between interfaces
      real*8,  intent(in) :: fracis(i1:i2, k1:k2, numSpecies) ! insoluble fraction of tracer
!
! !INPUT/OUTPUT PARAMETERS:
                   ! tracer array including moisture (mixing ratio)
      real*8, intent(inOut) :: qq (i1:i2, k1:k2, numSpecies)
!
! !DESCRIPTION:
! This routine performs convective transport of trace species.
! Note that we assume that the tracers are in a moist mixing ratio.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: SMALL = 1.0d-36
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: il, ik, ic
      integer :: km1, kp1
      real*8  :: avg_pbl    ! average mixing ratio in pbl
      real*8  :: cabv       ! mix ratio of constituent above
      real*8  :: cbel       ! mix ratio of constituent below
      real*8  :: cdifr      ! normalized diff. between cabv and cbel
      real*8  :: fluxin     ! flux coming into  each box at a k level
      real*8  :: fluxout    ! flux going out of each box at a k level
      real*8  :: maxc
      real*8  :: minc
      real*8  :: scav       ! mixing ratio of a scavenged tracer
      real*8  :: sqrt_fisg
      real*8  :: temp_conu  !.sds - new temp variable for updraft concen calc

      real*8  :: chat  (il2g, k1:k2) ! mix ratio in env.      at interfaces
      real*8  :: conu  (il2g, k1:k2) ! mix ratio in updraft   at interfaces
      real*8  :: cond  (il2g, k1:k2) ! mix ratio in downdraft at interfaces
      real*8  :: dcondt(il2g, k1:k2) ! gathered tend.  array
      real*8  :: edi   (il2g, k1:k2) ! gathered downdraft entrainment
      real*8  :: fisg  (il2g, k1:k2) ! gathered insoluble frac. of tracer
      real*8  :: xconst(il2g, k1:k2) ! gathered tracer array

      ICLOOP: do ic = 1, numSpecies
         if (isFixedConcentration(ic)) cycle icloop

         ! ------------------------------------------------
         ! Gather up the constituent and set tend. to zero.
         ! ------------------------------------------------
         do ik = k1, k2
            do il = 1, il2g
               fisg  (il,ik) = fracis(ideep(il),ik,ic)
               xconst(il,ik) = qq    (ideep(il),ik,ic)
            end do
         end do

         ! -----------------------------------------
         ! From now on work only with gathered data.
         ! -----------------------------------------

         ! ----------------------------------------------------
         ! Interpolate environment tracer values to interfaces.
         ! ----------------------------------------------------

         IKLOOP1: do ik = k1, k2
            km1 = Max (k1, ik-1)
            kp1 = Min (k2, ik+1)

            edi(:,ik) = mdi(:,ik) - mdi(:,kp1)

            ILLOOP1: do il = 1, il2g
               minc = Min (xconst(il,km1), xconst(il,ik))
               maxc = Max (xconst(il,km1), xconst(il,ik))

               if (minc < 0.0d0) then
                  cdifr = 0.0d0
               else
                  cdifr = Abs (xconst(il,ik) - xconst(il,km1)) /  &
                          Max (maxc, SMALL)
               end if

               if (cdifr > 1.0d-6) then
                  ! -------------------------------------------------------
                  ! The two layers differ significantly, so use a geometric
                  ! averaging procedure.
                  ! -------------------------------------------------------
                  cabv = Max (xconst(il,km1), (maxc * 1.0d-12))
                  cbel = Max (xconst(il,ik),  (maxc * 1.0d-12))
                  chat(il,ik) = Log (cabv/cbel) / (cabv - cbel) * cabv * cbel
               else
                  ! ---------------------------------------------------
                  ! The two layers have only a small difference, so use
                  ! arithmetic mean.
                  ! ---------------------------------------------------
                  chat(il,ik) = 0.5d0 * (xconst(il,ik) + xconst(il,km1))
               end if

               ! -------------------------------------
               ! Provisional up and down draft values.
               ! -------------------------------------
               conu(il,ik) = chat(il,ik)
               cond(il,ik) = chat(il,km1)

               ! ------------------
               ! Provisional tends.
               ! ------------------
               dcondt(il,ik) = 0.0d0
            end do ILLOOP1
         end do IKLOOP1

         where (edi(:,:) < 0.0d0) edi(:,:) = 0.0d0

         ILLOOP2: do il = 1, il2g
            ! ----------------------------------------------------------------------
            ! Find the mixing ratio in the downdraft, top of atmosphere down.
            ! NOTE: mass flux   in downdraft (mdi) will be zero or negative.
            !       entrainment in downdraft (edi) will be zero or positive.
            ! ----------------------------------------------------------------------
            IKLOOPDD: do ik = 2, k2
               if ((mdi(il,ik-1) < -xmbsth) .or. (edi(il,ik) > xmbsth)) then
                  cond(il,ik) = (cond(il,ik-1)*mdi(il,ik-1)  &
                          - xconst(il,ik)*edi(il,ik))  &
                         / (mdi(il,ik-1) - edi(il,ik))
               end if
            end do IKLOOPDD

            ! ------------------------------------------------------
            ! Calculate updrafts with scavenging from bottom to top.
            ! Include the downdrafts.
            ! ------------------------------------------------------

            ! -------------------------------------
            ! Do the bottom most levels in the pbl.
            ! -------------------------------------

            if (Sum (dpi(il,k2:pbli(il):-1)) == 0.0d0)  then
               err_msg = 'Problem in convectiveTransport.'
               PRINT*, err_msg
               STOP
            endif

            avg_pbl = Sum (xconst(il,k2:pbli(il):-1) * dpi(il,k2:pbli(il):-1)) /  &
                      Sum (dpi(il,k2:pbli(il):-1))

            sqrt_fisg = Sqrt (fisg(il,pbli(il)))

            scav    = avg_pbl * (1.0d0 -  sqrt_fisg)

            conu(il,pbli(il)) = avg_pbl * sqrt_fisg

            fluxin  = (mui(il,pbli(il)) + mdi(il,pbli(il))) *  &
                      Min (chat(il,pbli(il)), xconst(il,pbli(il)-1))  &
                        - (mdi(il,pbli(il)) * cond(il,pbli(il)))

            fluxout = mui(il,pbli(il)) * conu(il,pbli(il)) +  &
                          mui(il,pbli(il)) * scav

            dcondt(il,k2) = (fluxin - fluxout) /  &
                          Sum (dpi(il,k2:pbli(il):-1))

            if (avg_pbl > 0.0d0) then
               do ik = k2, pbli(il), -1
                  qq(ideep(il),ik,ic) = (avg_pbl + (dcondt(il,k2)*delt)) /  &
                                         avg_pbl * xconst(il,ik)
               end do
            end if

            ! ---------------------------
            ! Loop over all other levels.
            ! ---------------------------
            IKLOOP2: do ik = pbli(il)-1, 2, -1
               kp1 = ik + 1
               km1 = ik - 1

               ! -----------------------------
               ! Find mixing ratio in updraft.
               ! -----------------------------
               if ( (mui(il,kp1) - dui(il,ik) + eui(il,ik)) > xmbsth) then
                  sqrt_fisg   = Sqrt (fisg(il,ik))

                  if ( abs(mui(il,kp1) + eui(il,ik)) .ge. xmbsth) then
                     !... mix updraft and entrainment, then detrain this concentration
                     temp_conu  = ( mui(il,kp1) * conu(il,kp1) +  &
                                   eui(il,ik)  * xconst(il,ik) ) /  &
                                  (mui(il,kp1) + eui(il,ik))
                     !... if updraft and entrainment flux sum to 0.0, then take average concen
                  else
                     temp_conu  = (xconst(il,ik) + conu(il,kp1)) / 2.0d0
                  endif
!
                  scav = ((mui(il,kp1) * conu(il,kp1) -  &
                           dui(il,ik) * temp_conu ) *  &
                          (1.0d0 - fisg(il,ik)) +  &
                          eui(il,ik)  * xconst(il,ik) *  &
                          (1.0d0 - sqrt_fisg)) /  &
                         (mui(il,kp1) - dui(il,ik) + eui(il,ik))

                  conu(il,ik) = ((mui(il,kp1) * conu(il,kp1) -  &
                           dui(il,ik) * temp_conu ) *  &
                          fisg(il,ik) +  &
                          eui(il,ik)  * xconst(il,ik) *  &
                          sqrt_fisg) /  &
                         (mui(il,kp1) - dui(il,ik) + eui(il,ik))
               else
                  conu(il,ik) = xconst(il,ik)
                  scav        = 0.0d0
               end if

               ! -------------------------------------------------------
               ! Calculate fluxes into and out of box.  With scavenging
               ! included the net flux for the whole column is no longer
               ! guaranteed to be zero.
               ! Include the downdrafts.
               ! -------------------------------------------------------
               fluxin  = mui(il,kp1) * conu(il,kp1) +  &
                        (mui(il,ik)+mdi(il,ik))  &
                        * Min (chat(il,ik), xconst(il,km1))  &
                        - (mdi(il,ik) * cond(il,ik))

               fluxout = mui(il,ik)  * conu(il,ik) +  &
                        (mui(il,kp1) + mdi(il,kp1))  &
                        * Min (chat(il,kp1), xconst(il,ik)) +  &
                        (mui(il,ik) + mui(il,kp1)) * 0.5d0 * scav  &
                        - (mdi(il,kp1) * cond(il,kp1))

               dcondt(il,ik) = (fluxin - fluxout) / dpi(il,ik)

               ! --------------------------------------------
               ! Update and scatter data back to full arrays.
               ! --------------------------------------------
               qq(ideep(il),ik,ic) = xconst(il,ik) + (dcondt(il,ik) * delt)

            end do IKLOOP2
         end do ILLOOP2
      end do ICLOOP

      return

      end subroutine convectiveTransport
!EOC
!------------------------------------------------------------------------------
end module convectiveTransport_mod
