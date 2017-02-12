#include "symbol.inc"
      SUBROUTINE FTFLOP(L,IFLOP)
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
!***********************************************************************
! RCS:  $Id: main.F,v 1.1 1997/07/23 13:58:44 kresse
!                                                                      *
! Estimates the number of floating point operations (IFLOP) needed for *
! a FFT of length L according to the factorization from routine FTFACT *
!                                                                      *
!***********************************************************************

      DIMENSION IFAC(19)

! Single factor transforms according to structure of SPASSM:
      IF (L<=5) THEN
         IF (L==2) IFLOP=4
         IF (L==3) IFLOP=16
         IF (L==4) IFLOP=16
         IF (L==5) IFLOP=36
         RETURN
      END IF
! Get the factorization:
      CALL FTFACT(L,IFAC)
      IFLOP=0
! Count operations according to the structure of FPASSM/IPASSM:
      LA=1
      DO 30 I=2,IFAC(1)+1
         IFACT=IFAC(I)
         M=L/IFACT
         DO 10 J=1,LA
            IF (IFACT==2) IFLOP=IFLOP+4
            IF (IFACT==3) IFLOP=IFLOP+16
            IF (IFACT==4) IFLOP=IFLOP+16
            IF (IFACT==5) IFLOP=IFLOP+36
   10    CONTINUE
         IF (LA/=M) THEN
            DO 20 J=LA+1,M,LA
             DO 20 JJ=1,LA
              IF (IFACT==2) IFLOP=IFLOP+10
              IF (IFACT==3) IFLOP=IFLOP+28
              IF (IFACT==4) IFLOP=IFLOP+34
              IF (IFACT==5) IFLOP=IFLOP+52
   20       CONTINUE
         ENDIF
         LA=LA*IFACT
   30 CONTINUE
! That is all ...
      RETURN
      END

      FUNCTION  FTFTOT(NGPTAR)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NGPTAR(3)

      FTFTOT=0
      call ftflop(ngptar(1),iflop)
      FTFTOT=FTFTOT+iflop*ngptar(2)*ngptar(3)
      call ftflop(ngptar(2),iflop)
      FTFTOT=FTFTOT+iflop*ngptar(1)*ngptar(3)
      call ftflop(ngptar(3),iflop)
      FTFTOT=FTFTOT+iflop*ngptar(1)*ngptar(2)
      RETURN
      END



      SUBROUTINE FTFACT(L,IFAC)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
!***********************************************************************
!                                                                      *
! Factorization routine that first extracts all factors of 4, then all *
! factors of 2, all factors of 3 and all factors of 5 (allowed maximum *
! are 18 factors). On output IFAC(1) contains the number of factors of *
! L (or -99 in case of error), IFAC(2)...IFAC(IFAC(1)+1) the factors.  *
! One nice feature has been included in this routine which deviates a  *
! little from the conventions used in SCILIB routine CFTFAX: If L<0 we *
! factorize ABS(L) and take some other order of the factors: Instead   *
! of taking first all factors 4, then all factors 2, 3, 5 we rearrange *
! the factors so that the last factor is the largest of all and that   *
! the first factor is the second largest one! This has a quite simple  *
! reason: The routines CFFTM etc. use for the inner passes which swap  *
! between portions of the WORK array the highest efficient layout for  *
! the data arrangement which garuantees highest performances. But the  *
! first and the last pass are not necessarily as efficient as the all  *
! other passes (depending on the increments ...), hence we must assume *
! that they are 'low performance passes'. Our choice for the factors   *
! tries to minimize the number of operations to be used for these two  *
! passes and hence to minimize the 'slow-down' due to these passes.    *
!                                                                      *
!***********************************************************************

! The following parameter decides whether to use factors 4 or not!
! The number of necessary floating point operations is always minimal
! if we use only factors 2, but on some machines factors 4 might have
! the advantage that more work is done in one block/loop what could
! reduce overhead for starting up loops/pipelines (vector processors!)
      PARAMETER(NOFOUR=0)
! FTFACT   m u s t   be called by all routines doing FFTs. So it is an
! appropriate place to decide some other things ... . When do you want
! to norm your output? For ISIGN=+1/-1? Not at all? Set INORM!
      PARAMETER(INORM=0)

      DIMENSION IFAC(19)

      COMMON /FTNRM/ IFTNRM
!
      IFTNRM=INORM
      IORD=0
! Negative L means reordering of the factors (makes only sense if more
! than one factor is contained in L, single factor is trivial) ... :
      IF (L<-5) IORD=1
! L must be larger than 1!
      IF (ABS(L)>1) GOTO 10
      IFAC(1)=-99
      RETURN
   10 NN=ABS(L)
      K=1
! Test for factors of 4:
      IF (NOFOUR/=0) GOTO 30
   20 IF (MOD(NN,4)/=0) GOTO 30
      K=K+1
      IF (K>19) GOTO 90
      IFAC(K)=4
      NN=NN/4
      IF (NN==1) GOTO 80
      GOTO 20
! Test for (extra) factor(s) of 2
   30 IF (MOD(NN,2)/=0) GOTO 40
      K=K+1
      IF (K>19) GOTO 90
! The smallest possible amount of floating point operations can be
! achieved if we use for the factors the order 2,4,3,5 (increasing
! operation count per elementary step) - this has been found in an
! empirical 'FLOP-counting' experiment for all possible orders ...:
      IF (NOFOUR==0) THEN
         DO 35 IFOUR=K-1,2,-1
   35    IFAC(IFOUR+1)=IFAC(IFOUR)
         IFAC(2)=2
      ELSE
         IFAC(K)=2
      ENDIF
      NN=NN/2
      IF (NN==1) GOTO 80
      IF (NOFOUR/=0) GOTO 30
! Test for factors of 3
   40 IF (MOD(NN,3)/=0) GOTO 50
      K=K+1
      IF (K>19) GOTO 90
      IFAC(K)=3
      NN=NN/3
      IF (NN==1) GOTO 80
      GOTO 40
! Test for factors of 5
   50 IF (MOD(NN,5)/=0) GOTO 90
      K=K+1
      IF (K>19) GOTO 90
      IFAC(K)=5
      NN=NN/5
      IF (NN==1) GOTO 80
      GOTO 50
! IFAC(1) now contains number of factors
   80 IFAC(1)=K-1
      IF (IORD/=0) GOTO 100
      RETURN
! Error: other factors than 2,3,4 or 5 or more than 18 factors!
   90 IFAC(1)=-99
      RETURN
  100 CONTINUE
      L1FAC=1
      L2FAC=1
! Largest factor:
      DO 110 I=2,K
         IF (L1FAC<IFAC(I)) THEN
            L1FAC=IFAC(I)
            INDEX1=I
         END IF
  110 CONTINUE
! Second largest factor:
      DO 120 I=2,K
         IF (I==INDEX1) GOTO 120
         IF (L2FAC<IFAC(I)) THEN
            L2FAC=IFAC(I)
            INDEX2=I
         END IF
  120 CONTINUE
! Reorder:
      IF (INDEX2/=2) THEN
         ISAVE=IFAC(2)
         IFAC(2)=IFAC(INDEX2)
         IFAC(INDEX2)=ISAVE
      END IF
      IF (INDEX1==2) INDEX1=INDEX2
      IF (INDEX1/=K) THEN
         ISAVE=IFAC(K)
         IFAC(K)=IFAC(INDEX1)
         IFAC(INDEX1)=ISAVE
      END IF
      RETURN
      END
