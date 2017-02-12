#include "symbol.inc"
      MODULE hamil
      USE prec
      CONTAINS
!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions
!  FFT and non-local projections of wavefunctions must be supplied
!***********************************************************************

      SUBROUTINE ECCP(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV, CE)
      USE prec
      USE wave
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)


      TYPE (wavefun1) :: W1,W2
      TYPE (wavedes1) :: WDES1
      TYPE (grid_3d)  :: GRID

      INTEGER NGVECTOR, ISPINOR
      GDEF      CNL
      RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!=======================================================================
! calculate the local contribution
!=======================================================================
#if defined (use_cray_ptr)
      INTEGER            :: m, mm1, mm2, mm3, ispinor1, ispinor2, &
                            mgvector, mrspinors, mmp, mmplwv

      REAL (q)           :: w_datake(*)
      POINTER               (w_datake_p, w_datake)
#ifdef gammareal
      REAL (q)           :: w1_cp(*), w2_cp(*)
      COMPLEX (q)        :: w1_cr(*), w2_cr(*), w1_cw(*), w2_cw(*)
#else
      COMPLEX (q)        :: w1_cp(*), w1_cr(*), w1_cw(*), &
                            w2_cp(*), w2_cr(*), w2_cw(*)
#endif
      POINTER               (w1_cp_p, w1_cp), (w2_cp_p, w2_cp), &
                            (w1_cr_p, w1_cr), (w2_cr_p, w2_cr), &
                            (w1_cw_p, w1_cw), (w2_cw_p, w2_cw)
#endif
      CLOCAL=0
      NGVECTOR=WDES1%NGVECTOR

#if defined (use_cray_ptr)
      w1_cp_p = LOC (W1%CPROJ)
      w1_cw_p = LOC (W1%CW)
      w1_cr_p = LOC (W1%CR)

      w2_cp_p = LOC (W2%CPROJ)
      w2_cw_p = LOC (W2%CW)
      w2_cr_p = LOC (W2%CR)

      mgvector  = wdes1%ngvector
      mrspinors = wdes1%nrspinors

      mmp    = grid%rl%np
      mmplwv = grid%mplwv

      DO ispinor2 = 0, mrspinors - 1
         mm2 = ispinor2 * mmplwv
         DO ispinor1 = 0, mrspinors - 1
            mm1 = ispinor1 * mmplwv
            mm3 = 1 + ispinor1 + 2 * ispinor2
            DO m = 1, mmp
               clocal = clocal + sv(m,mm3) * w1_cr(m+mm1) * CONJG(w2_cr(m+mm2))
            END DO
         END DO
      END DO
#else
      DO ISPINOR =0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1
         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *GRID%MPLWV
            MM_=M+ISPINOR_*GRID%MPLWV
            CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
         ENDDO
      ENDDO
      ENDDO
#endif

      CLOCAL=CLOCAL/GRID%NPLWV
!=======================================================================
! kinetic energy contribution
!=======================================================================
      CKIN=0
#if defined (use_cray_ptr)
      DO ispinor2 = 0, mrspinors - 1
         mm2 = ispinor2 * mgvector
         w_datake_p = LOC (WDES1%DATAKE(1,ispinor2+1))
         DO m = 1, mgvector
            ckin = ckin + w1_cw(m+mm2) * CONJG (w2_cw(m+mm2)) * w_datake(m)
         END DO
      END DO
#else

      DO ISPINOR=0,WDES1%NRSPINORS-1
         DO M=1,NGVECTOR
            MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!           CKIN=CKIN+W1%CW(MM)*CONJG(W2%CW(MM))*WDES1%DATAKE(M)
            CKIN=CKIN+W1%CW(MM)*CONJG(W2%CW(MM))*WDES1%DATAKE(M,ISPINOR+1)
!-MM- end of alterations
         ENDDO
      ENDDO
#endif
!=======================================================================
! non local contribution
!=======================================================================
      CNL =0
      NPRO=0

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1

      NPRO =ISPINOR *WDES1%NPRO/2
      NPRO_=ISPINOR_*WDES1%NPRO/2

      NIS =1
      DO NT=1,WDES1%NTYP
        LMMAXC=WDES1%LMMAX(NT)
        IF (LMMAXC==0) GOTO 310
        DO NI=NIS,WDES1%NITYP(NT)+NIS-1
         CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
         NPRO = LMMAXC+NPRO
         NPRO_= LMMAXC+NPRO_
        ENDDO
  310   NIS = NIS+WDES1%NITYP(NT)
      ENDDO
      ENDDO
      ENDDO spinor

      CE=GREAL(CLOCAL+CKIN+CNL)
      CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

      RETURN
      END SUBROUTINE


      END MODULE

!************************* SUBROUTINE ECCP_NL   ************************
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions; non local part only
! for one ion only
! I have put this in a seperate routine because optimization
! is than easier
!***********************************************************************

      SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      GDEF      CNL
      OVERLAP CDIJ(LMDIM,LMDIM)
      GDEF CPROJ1(LMMAXC),CPROJ2(LMMAXC)

      DO L=1,LMMAXC
      DO LP=1,LMMAXC
        CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*GCONJG(CPROJ2(L))
      ENDDO; ENDDO
      END SUBROUTINE

!************************* SUBROUTINE VHAMI  ***************************
!
! this subroutine calculates the product of the local potential
! with a wavefunction stored in CR and returns the result in CVR
!
!***********************************************************************

      SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
      USE prec
      USE mgrid
      USE wave
      IMPLICIT NONE

      TYPE (grid_3d)     GRID
      TYPE (wavedes1)    WDES1

      RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential

      COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
! local variables
      REAL(q) RINPLW
      INTEGER ISPINOR,ISPINOR_,M,MM,MM_

      RINPLW=1._q/GRID%NPLWV
!
! calculate the local contribution (store result in CWORK1)
! 
      IF (WDES1%NRSPINORS==1) THEN
         DO M=1,GRID%RL%NP
            CVR(M)= SV(M,1) *CR(M)*RINPLW
         ENDDO
      ELSE
         CVR(1:GRID%MPLWV*2)=0
         DO ISPINOR =0,1
         DO ISPINOR_=0,1
            DO M=1,GRID%RL%NP
               MM =M+ISPINOR *GRID%MPLWV
               MM_=M+ISPINOR_*GRID%MPLWV
               CVR(MM)= CVR(MM)+ SV(M,1+ISPINOR_+2*ISPINOR) *CR(MM_)*RINPLW
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE VHAMIL


!************************* SUBROUTINE HAMILT ***************************
!
! this subroutine calculates the H acting onto a wavefuntion
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!   H- EVALUE Q |phi>
! where Q is the nondiagonal part of the overlap matrix S
!***********************************************************************

      SUBROUTINE HAMILT( &
     &    WDES1,W1,NONLR_S,NONL_S,GRID,  LREAL,EVALUE, &
     &    LMDIM,CDIJ,CQIJ, SV,CH)
      USE prec

      USE mpimy
      USE mgrid
      USE wave
      USE nonl
      USE nonlr
      IMPLICIT NONE


      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun1)    W1
      TYPE (wavedes1)    WDES1
      INTEGER  LMDIM, NGVECTOR, ISPINOR, ISPINOR_, MM, MM_

      RGRID   SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
              CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CH(WDES1%NPL)
      REAL(q)      EVALUE
      LOGICAL   LREAL
! local variables
      REAL(q) RINPLW; INTEGER M
      COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS)

      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR
!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
      CALL VHAMIL(WDES1,GRID,SV(1,1),W1%CR(1),CWORK1(1)) 
!=======================================================================
! non-local contribution in real-space
!=======================================================================
      IF (LREAL) THEN
         CALL RACC(NONLR_S,WDES1,W1, LMDIM,CDIJ,CQIJ,EVALUE, CWORK1)

         DO ISPINOR=0,WDES1%NRSPINORS-1
            CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%MPLWV),CH(1+ISPINOR*NGVECTOR),GRID,.FALSE.)
            DO M=1,NGVECTOR
               MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!              CH(MM)=CH(MM)+W1%CW(MM)* WDES1%DATAKE(M)
               CH(MM)=CH(MM)+W1%CW(MM)* WDES1%DATAKE(M,ISPINOR+1)
!-MM- end of alterations
            ENDDO
         ENDDO
      ELSE
!=======================================================================
! calculate the non local contribution in reciprocal space
!=======================================================================
         CALL VNLACC(NONL_S,WDES1,W1, LMDIM,CDIJ,CQIJ,EVALUE, CH)
         DO ISPINOR=0,WDES1%NRSPINORS-1
            CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*GRID%MPLWV),CH(1+ISPINOR*NGVECTOR),GRID,.TRUE.)
            DO M=1,NGVECTOR
               MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!              CH(MM)=W1%CW(MM)* WDES1%DATAKE(M)+CH(MM)
               CH(MM)=W1%CW(MM)* WDES1%DATAKE(M,ISPINOR+1)+CH(MM)
!-MM- end of alterations
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE

!************************* SUBROUTINE HAMILTMU *************************
!
! this subroutine calculates the H acting onto a set of wavefuntions
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!   H- EVALUE Q |phi>
! where Q is the nondiagonal part of the overlap matrix S
!***********************************************************************

      SUBROUTINE HAMILTMU( &
     &    WDES1,W1,NONLR_S,NONL_S,GRID,  LREAL,EVALUE, &
     &    LMDIM,CDIJ,CQIJ, SV,CH,LD, NSIM, LDO)
      USE prec

      USE mpimy
      USE mgrid
      USE wave
      USE nonl
      USE nonlr
      IMPLICIT NONE


      INTEGER NSIM,NP,LD
      INTEGER LMDIM, NGVECTOR, ISPINOR, ISPINOR_, MM, MM_ 
      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun1)    W1(NSIM)
      TYPE (wavedes1)    WDES1

      RGRID      SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      OVERLAP    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CH(LD,NSIM)
      REAL(q)    EVALUE(NSIM)
      LOGICAL LREAL
      LOGICAL LDO(NSIM)
! local variables
      REAL(q) RINPLW; INTEGER M

#if defined (use_malloc)
      COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)
      POINTER               (cwork1_p, cwork1)
#else
! I would prefer to allocate the array CWORK1 from the stack
! this however fails for large systems with no obvious reason
! hence allocate it from the heap
!      COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)
#endif

#if defined (use_cray_ptr)
      INTEGER            :: mgvector, mrspinors, ispinor1, &
                            ispinor2, mm1, mm2, mm3, mmp, mmplwv

      REAL (q)           :: w_datake(*)
      POINTER               (w_datake_p, w_datake)
      COMPLEX (q)        :: w1_cr(*), w1_cw(*)
      POINTER               (w1_cr_p, w1_cr), (w1_cw_p, w1_cw)
#endif

#if defined (use_malloc)
      cwork1_p = malloc(2*q*grid%mplwv*wdes1%nrspinors*nsim)
#else
      ALLOCATE(CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM))
#endif

      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR

#if defined (use_cray_ptr)
      mgvector  = wdes1%ngvector
      mrspinors = wdes1%nrspinors

      mmp    = grid%rl%np
      mmplwv = grid%mplwv
#endif
!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
#if defined (use_cray_ptr)
      DO np = 1, nsim
         IF (ldo(np)) THEN
            w1_cr_p = LOC (W1(NP)%CR)

            cwork1(:,np) = (0.0_q, 0.0_q)

            DO ispinor2 = 0, mrspinors - 1
               mm2 = ispinor2 * mmplwv
               DO ispinor1 = 0, mrspinors - 1
                  mm1 = ispinor1 * mmplwv
                  mm3 = 1 + ispinor1 + 2 * ispinor2
                  DO m = 1, mmp
                     cwork1(m+mm2,np) = cwork1(m+mm2,np) + sv(m,mm3) *  &
                                        w1_cr(m+mm1) * rinplw
                  END DO
               END DO
            END DO
         END IF
      END DO
#else
      DO NP=1,NSIM
        IF ( LDO(NP) ) THEN
        CWORK1(:,NP)  = 0
        DO ISPINOR =0,WDES1%NRSPINORS-1
        DO ISPINOR_=0,WDES1%NRSPINORS-1
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              MM_=M+ISPINOR_*GRID%MPLWV
              CWORK1(MM,NP)=  CWORK1(MM,NP)+SV(M,1+ISPINOR_+2*ISPINOR) *W1(NP)%CR(MM_)*RINPLW
           ENDDO
        ENDDO
        ENDDO
        ENDIF
      ENDDO
#endif
!=======================================================================
! non-local contribution in real-space
!=======================================================================
      IF (LREAL) THEN
         CALL RACCMU(NONLR_S,WDES1,W1, LMDIM,CDIJ,CQIJ,EVALUE,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)

         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN
#if defined (use_cray_ptr)
               w1_cw_p = LOC (W1(NP)%CW)
               DO ispinor2 = 0, mrspinors - 1
                  mm2 = ispinor2 * mgvector
                  CALL fftext (mgvector, wdes1%nindpw(1), &
                               cwork1(1+ispinor2*wdes1%mplwv,np), &
                               ch(1+ispinor2*mgvector,np), grid, .FALSE.)
                  w_datake_p = LOC (wdes1%datake(1,ispinor2+1))
                  DO m = 1, mgvector
                     ch(m+mm2,np) = ch(m+mm2,np) + w1_cw(m+mm2) * w_datake(m)
                  END DO
               END DO
#else

            DO ISPINOR=0,WDES1%NRSPINORS-1
               CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.FALSE.)
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!                 CH(MM,NP)=CH(MM,NP)+W1(NP)%CW(MM)*WDES1%DATAKE(M)
                  CH(MM,NP)=CH(MM,NP)+W1(NP)%CW(MM)*WDES1%DATAKE(M,ISPINOR+1)
!-MM- end of alterations
               ENDDO
            ENDDO
#endif
            ENDIF
         ENDDO
!=======================================================================
! calculate the non local contribution in reciprocal space
!=======================================================================
      ELSE

         DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            CALL VNLACC(NONL_S,WDES1,W1(NP), LMDIM,CDIJ,CQIJ,EVALUE(NP),CH(1,NP))
#if defined (use_cray_ptr)
            w1_cw_p = LOC (W1(NP)%CW)
            DO ispinor2 = 0, mrspinors - 1
               mm2 = ispinor2 * mgvector
               CALL fftext (mgvector, wdes1%nindpw(1), & 
                            cwork1(1+ispinor2*wdes1%mplwv,np), &
                            ch(1+ispinor2*mgvector,np), grid, .TRUE.)
               w_datake_p = LOC (wdes1%datake(1,ispinor2+1))
               DO m = 1, mgvector
                  ch(m+mm2,np) = w1_cw(m+mm2) * w_datake(m) + ch(m+mm2,np)
               END DO
            END DO
#else
            DO ISPINOR=0,WDES1%NRSPINORS-1
               CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!                 CH(MM,NP)=W1(NP)%CW(MM)* WDES1%DATAKE(M)+CH(MM,NP)
                  CH(MM,NP)=W1(NP)%CW(MM)* WDES1%DATAKE(M,ISPINOR+1)+CH(MM,NP)
!-MM- end of alterations
               ENDDO
            ENDDO
#endif
        ENDIF
        ENDDO
      ENDIF

#if defined (use_malloc)
      CALL free (cwork1_p)
#else
      DEALLOCATE(CWORK1)
#endif

      RETURN
      END SUBROUTINE
