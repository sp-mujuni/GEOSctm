!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: CTM_pTracersGridCompMod - The pTracers Grid Component
!
! !INTERFACE:
!
   Module CTM_pTracersGridCompMod
!
! !USES:
!

   use ESMF

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!   This is a stub for the {\tt pTracers} gridded component.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for pTracers Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  28Jan2013  Kouatchou  First crack.
!
!EOP
!-------------------------------------------------------------------------

  if ( present(rc) ) rc = 0
  
  end subroutine SetServices

end Module CTM_pTracersGridCompMod
