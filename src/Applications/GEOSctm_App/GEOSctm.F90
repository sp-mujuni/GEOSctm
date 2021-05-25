!-------------------------------------------------------------------------
!         NASA/GSFC, Software Systems Support Office, Code 610.3         !
!-------------------------------------------------------------------------
! *********************************************************************
! *****                      Main Program                          ****
! *****                     GEOS-5 as a CTM                        ****
! *********************************************************************

#define I_AM_MAIN
#include "MAPL_Generic.h"
!----------------------------------------------------------------------
!BOP

Program GEOSctm
   
! !USES:
   use MAPL

   use GEOS_ctmGridCompMod, only: ROOT_SetServices    => SetServices
   implicit NONE

   character(len=*), parameter :: Iam = 'GEOSctm'

   type (MAPL_FlapCLI) :: cli
   type (MAPL_CapOptions) :: cap_options

   integer :: status

!EOP
!----------------------------------------------------------------------
!BOC

   cli = MAPL_FlapCLI(description = 'GEOS CTM', & 
                      authors     = 'GMAO')
   cap_options = MAPL_CapOptions(cli)
   cap = MAPL_Cap('CTM', ROOT_SetServices, cap_options = cap_options)
   call cap%run(_RC)

end program GEOSctm
!EOC
!----------------------------------------------------------------------
