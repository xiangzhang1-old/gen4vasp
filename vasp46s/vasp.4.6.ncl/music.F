

#include "symbol.inc"
      SUBROUTINE MUSIC(YES,IWHICH)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! You want some hymnus at the end? --- Voila!!
      LOGICAL YES
      INTEGER SYSTEM
      CHARACTER*80 HYMNUS
! Where to find it / which command(s) to use ...
      IF (IWICH==0) HYMNUS='/tph21/u/furth/C/bell -s h'
      IF (IWICH==1) HYMNUS='/tph21/u/furth/C/bell -s o'
! You do not want??
      IF (.NOT.YES) RETURN

! Sorry, currently UNIX-support only ..., maybe commment
! out next line if function system not provided/accesible:
      L=LENGTH(HYMNUS)
      IERR=SYSTEM(HYMNUS(1:L))

      RETURN
      END
