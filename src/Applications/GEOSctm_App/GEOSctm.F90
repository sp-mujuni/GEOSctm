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
   use MAPL_Mod

   use GEOS_ctmGridCompMod, only: ROOT_SetServices    => SetServices

   implicit NONE

   integer :: STATUS
   logical :: am_I_root

   character(len=*), parameter :: Iam = 'GEOSctm'

!EOP
!----------------------------------------------------------------------
!BOC
   
    call MAPL_CAP(ROOT_SetServices, AmIRoot = am_I_root, rc=STATUS)

    call exit(STATUS)

end program GEOSctm
!EOC
!----------------------------------------------------------------------
