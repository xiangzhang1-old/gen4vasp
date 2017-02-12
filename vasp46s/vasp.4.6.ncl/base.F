#include "symbol.inc"
!************************************************************************
! RCS:  $Id: base.F,v 1.2 2001/02/20 14:44:56 kresse Exp $
!
! this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
      END MODULE

      MODULE BASE
      USE prec
      INCLUDE "base.inc"
      CONTAINS
!
! small subroutine which tries to give good dimensions for 1 dimension
!
      SUBROUTINE MAKE_STRIDE (N)
      INTEGER N,NEW
      INTEGER, PARAMETER :: NGOOD=16

      NEW=(N+NGOOD)/NGOOD
      NEW=NEW*NGOOD+1
      N=NEW

      END SUBROUTINE

      END MODULE








