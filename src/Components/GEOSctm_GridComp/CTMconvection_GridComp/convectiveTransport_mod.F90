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

      implicit none

      private 
      public  :: convectiveTransport
!
! !DESCRIPTION:
! Module for the convective transport.
!
!EOP
!-----------------------------------------------------------------------------
CONTAINS
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
