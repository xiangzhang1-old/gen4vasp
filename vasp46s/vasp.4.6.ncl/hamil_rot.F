#include "symbol.inc"
!======================================================================
! RCS:  $Id: hamil_rot.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  MODUL contains some routines only used for all bands simultaneous
!  update scheme implemented in rot.f
!  dont worry about this routines if you do not use rot.f
!======================================================================

!************************* SUBROUTINE HAMNNP ***************************
!
! this subroutine calculates the Hamiltonian enclosed between
! one wavefunction and a set of wavefunctions
!   H_n' n = <C_n | H | C_n'>
!***********************************************************************

      SUBROUTINE HAMNNP(GRID,NONLR_S,NONL_S,W,WDES,CR, LREAL,SV,N,NK, &
         LMDIM,CDIJ,CHAM)
      USE prec

      USE nonl
      USE nonlr
      USE mpimy
      USE mgrid
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun)     W
      TYPE (wavedes)     WDES

      LOGICAL LREAL

      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS)
      COMPLEX(q)   CR(GRID%RL%NP)
      RGRID     SV(GRID%RL%NP)
      GDEF      CHAM(WDES%NBANDS,WDES%NBANDS,NK)
! local work arrays
      GDEF      DWORK(WDES%NBANDS)
      COMPLEX(q),ALLOCATABLE :: CWORK1(:),CWORK2(:)
      GDEF   ,ALLOCATABLE :: CPROF(:)

      ALLOCATE(CWORK1(GRID%MPLWV),CWORK2(GRID%MPLWV),CPROF(WDES%NPROD))

      NPL=WDES%NPLWKP(NK)
!-----------------------------------------------------------------------
! calculate the local contribution
!-----------------------------------------------------------------------
      RINPLW=1._q/GRID%NPLWV
      DO M=1,GRID%RL%NP
        CWORK1(M)= REAL( SV(M) ,KIND=q) *CR(M)*RINPLW
      ENDDO
!-----------------------------------------------------------------------
! transforme back to reciprocal space and add kinetic energy
!-----------------------------------------------------------------------
      CALL FFTEXT(NPL,WDES%NINDPW(1,NK),CWORK1,CWORK2,GRID,.FALSE.)
      DO M=1,NPL
       CWORK2(M)=GCONJG(CWORK2(M)+W%CPTWFP(M,N,NK)* WDES%DATAKE(M,NK))
      ENDDO
!=======================================================================
! calculate the non local part
!=======================================================================
      CPROF=0

      NIS=1
      NPRO=0

      typ: DO NT=1,NONL_S%NTYP
      LMMAXC=NONL_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 100
      ion: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1

      DO L =1,LMMAXC
      DO LP=1,LMMAXC
        CPROF(L+NPRO)=CPROF(L+NPRO) + &
                      GCONJG(CDIJ(LP,L,NI)*W%CPROJ(LP+NPRO,N,NK))
      ENDDO
      ENDDO

      NPRO= LMMAXC+NPRO
      ENDDO ion
  100 NIS=NIS+NONL_S%NITYP(NT)
      ENDDO typ
!=======================================================================
! calculate the Hamilton-matrix defined as
! CHAM(n2,n1) = < C(n1) | H | C(n2) >
!=======================================================================
      CALL GGEMV( 'T',m_ NPL , WDES%NBANDS-N+1 ,one ,W%CPTWFP(1,N,NK), &
     &             m_ WDES%NRPLWV, CWORK2, 1 , zero ,  DWORK, 1)

      CALL GGEMV( 'T',  WDES%NPRO, WDES%NBANDS-N+1 ,one ,W%CPROJ(1,N,NK), &
     &             WDES%NPROD, CPROF, 1 , one ,  DWORK, 1)

      N1=N
      DO N2=N1,WDES%NBANDS
        CHAM(N2,N1,NK)= GREAL(DWORK(N2-N1+1))
        CHAM(N1,N2,NK)= GREAL(GCONJG(DWORK(N2-N1+1)))
      ENDDO

      DEALLOCATE(CWORK1,CWORK2,CPROF)
      RETURN
      END SUBROUTINE



!************************* SUBROUTINE HAMIL0 ***************************
!
! this subroutine calculates the gradient vector of the band N, if the
! eigenmatrix is not diagonal   i.e.
!   H | C_n >- H_n n' S | C_n' >
! especially  due to the overlap-operator this is a somewhat complicated
! operation
!
! CHAM contains the (possibly not diagonal) sub-space Hamilton matrix
! the  wavefunction must be given in reciprocal space C(..,N) and real
! space CR and the projection operatores must be stored in CPROJ(..,N)
! CACC contains the result in reciprocal space
!***********************************************************************

      SUBROUTINE HAMIL0(GRID,NONLR_S,NONL_S,W,WDES,CR, LREAL,SV,N,NK, &
         LMDIM,CQIJ,CDIJ,CHAM,  CACC)
      USE prec

      USE nonl
      USE nonlr
      USE mpimy
      USE mgrid
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (wavefun)     W

      LOGICAL LREAL

      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS),CQIJ(LMDIM,LMDIM,WDES%NIONS)
      COMPLEX(q)   CR(GRID%RL%NP)
      RGRID     SV(GRID%RL%NP)
      GDEF      CHAM(WDES%NBANDS,WDES%NBANDS)
      COMPLEX(q)   CACC(WDES%NRPLWV)
! work arrays
      GDEF      CWORK(WDES%NBANDS)

      COMPLEX(q),ALLOCATABLE :: CPTACC(:)
      GDEF   ,ALLOCATABLE :: CWORK2(:),CWORK3(:)

      ALLOCATE(CPTACC(GRID%MPLWV),CWORK2(WDES%NPROD),CWORK3(WDES%NPROD))

      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)
!=======================================================================
! calculate sum_i n' Q_ij j H_nn' C_in
!=======================================================================

      NPL=WDES%NPLWKP(NK)

!---First store the elements of the Hamiltonmatrix in the correct order
      DO NP=1,WDES%NBANDS
        CWORK(NP)=CHAM(N,NP)
      ENDDO
!-----------------------------------------------------------------------
!  now sum_n  -H_nn'  C_i n'     and   sum_n -H_nn' C_Gn'
!-----------------------------------------------------------------------
      CALL GGEMM( 'N', 'N' ,m_ NPL , 1 , WDES%NBANDS , -one, &
     &             W%CW(1,1,NK),m_ WDES%NRPLWV , CWORK , WDES%NBANDS , &
     &             zero , CACC ,m_ WDES%NRPLWV )

      CALL GGEMM( 'N', 'N' , WDES%NPRO , 1 , WDES%NBANDS  , -one, &
     &             W%CPROJ(1,1,NK),  WDES%NPROD , CWORK , WDES%NBANDS , &
     &             zero , CWORK3 , WDES%NPROD  )

!-----------------------------------------------------------------------
!  now sum_ij Q_ij  (sum_np -H_n np  C_i np)   + D_ij C_in
!-----------------------------------------------------------------------
      CWORK2=0

      NIS =1
      NPRO=0

      typ: DO NT=1,NONL_S%NTYP
      LMMAXC=NONL_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 100
      ion: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
           CWORK2(L+NPRO)=CWORK2(L+NPRO)+CDIJ(LP,L,NI)*W%CPROJ (LP+NPRO,N,NK) &
                                     +CQIJ(LP,L,NI)*CWORK3(LP+NPRO)
        ENDDO
        ENDDO
      NPRO= LMMAXC+NPRO
      ENDDO ion
  100 NIS = NIS+NONL_S%NITYP(NT)
      ENDDO typ
!=======================================================================
! calculate the local contribution (result in CPTACC)
!=======================================================================
      RINPLW=1._q/GRID%NPLWV
      DO M=1,GRID%RL%NP
        CPTACC(M)= REAL( SV(M) ,KIND=q) *CR(M)*RINPLW
      ENDDO
!=======================================================================
! non-local contribution using real-space projection add to CPTACC
!=======================================================================
      IF (LREAL) THEN
      CALL RACC0(NONLR_S,WDES1,CWORK2,CPTACC)
      ENDIF
!=======================================================================
! calculate the non local contribution using reciprocal space
!=======================================================================
      IF (.NOT.LREAL) THEN
        CALL VNLAC0(NONL_S,WDES1,CWORK2,CACC)
!=======================================================================
! calculate the total Hamiltonian
!=======================================================================
        CALL FFTEXT(NPL,WDES%NINDPW(1,NK),CPTACC,CACC,GRID,.TRUE.)
        DO M=1,NPL
          CACC(M)=CACC(M)+W%CW(M,N,NK)* WDES%DATAKE(M,NK)
        ENDDO
      ELSE
        CALL FFTEXT(NPL,WDES%NINDPW(1,NK),CPTACC,CACC,GRID,.TRUE.)
        DO M=1,NPL
          CACC(M)=CACC(M)+W%CW(M,N,NK)* WDES%DATAKE(M,NK)
        ENDDO
      ENDIF
      DEALLOCATE(CPTACC,CWORK2,CWORK3)

      RETURN
      END SUBROUTINE
