!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part         --------------------------------
!
!   charge density: full grid mode
!
!
!   charge density complex
!
!
!   wavefunctions: full grid mode
!
!
!   wavefunctions complex
!
!
!   common definitions
!





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
      COMPLEX(q)      CNL
      COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!=======================================================================
! calculate the local contribution
!=======================================================================
      CLOCAL=0
      NGVECTOR=WDES1%NGVECTOR

      DO ISPINOR =0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1
         DO M=1,GRID%RL%NP
            MM =M+ISPINOR *GRID%MPLWV
            MM_=M+ISPINOR_*GRID%MPLWV
            CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
         ENDDO
      ENDDO
      ENDDO

      CLOCAL=CLOCAL/GRID%NPLWV
!=======================================================================
! kinetic energy contribution
!=======================================================================
      CKIN=0

      DO ISPINOR=0,WDES1%NRSPINORS-1
         DO M=1,NGVECTOR
            MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!           CKIN=CKIN+W1%CPTWFP(MM)*CONJG(W2%CPTWFP(MM))*WDES1%DATAKE(M)
            CKIN=CKIN+W1%CPTWFP(MM)*CONJG(W2%CPTWFP(MM))*WDES1%DATAKE(M,ISPINOR+1)
!-MM- end of alterations
         ENDDO
      ENDDO
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

      CE=(CLOCAL+CKIN+CNL)
      

      RETURN
      END SUBROUTINE


      END MODULE

!************************* SUBROUTINE ECCP_NL   ************************
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions; non local part only
! for (1._q,0._q) ion only
! I have put this in a seperate routine because optimization
! is than easier
!***********************************************************************

      SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      COMPLEX(q)      CNL
      COMPLEX(q) CDIJ(LMDIM,LMDIM)
      COMPLEX(q) CPROJ1(LMMAXC),CPROJ2(LMMAXC)

      DO L=1,LMMAXC
      DO LP=1,LMMAXC
        CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*CONJG(CPROJ2(L))
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

      COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential

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

      COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
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
!              CH(MM)=CH(MM)+W1%CPTWFP(MM)* WDES1%DATAKE(M)
               CH(MM)=CH(MM)+W1%CPTWFP(MM)* WDES1%DATAKE(M,ISPINOR+1)
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
!              CH(MM)=W1%CPTWFP(MM)* WDES1%DATAKE(M)+CH(MM)
               CH(MM)=W1%CPTWFP(MM)* WDES1%DATAKE(M,ISPINOR+1)+CH(MM)
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

      COMPLEX(q)      SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      COMPLEX(q)    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CH(LD,NSIM)
      REAL(q)    EVALUE(NSIM)
      LOGICAL LREAL
      LOGICAL LDO(NSIM)
! local variables
      REAL(q) RINPLW; INTEGER M

! I would prefer to allocate the array CWORK1 from the stack
! this however fails for large systems with no obvious reason
! hence allocate it from the heap
!      COMPLEX(q) :: CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)


      ALLOCATE(CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM))

      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR

!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
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
!=======================================================================
! non-local contribution in real-space
!=======================================================================
      IF (LREAL) THEN
         CALL RACCMU(NONLR_S,WDES1,W1, LMDIM,CDIJ,CQIJ,EVALUE,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)

         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN

            DO ISPINOR=0,WDES1%NRSPINORS-1
               CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.FALSE.)
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!                 CH(MM,NP)=CH(MM,NP)+W1(NP)%CPTWFP(MM)*WDES1%DATAKE(M)
                  CH(MM,NP)=CH(MM,NP)+W1(NP)%CPTWFP(MM)*WDES1%DATAKE(M,ISPINOR+1)
!-MM- end of alterations
               ENDDO
            ENDDO
            ENDIF
         ENDDO
!=======================================================================
! calculate the non local contribution in reciprocal space
!=======================================================================
      ELSE

         DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            CALL VNLACC(NONL_S,WDES1,W1(NP), LMDIM,CDIJ,CQIJ,EVALUE(NP),CH(1,NP))
            DO ISPINOR=0,WDES1%NRSPINORS-1
               CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!                 CH(MM,NP)=W1(NP)%CPTWFP(MM)* WDES1%DATAKE(M)+CH(MM,NP)
                  CH(MM,NP)=W1(NP)%CPTWFP(MM)* WDES1%DATAKE(M,ISPINOR+1)+CH(MM,NP)
!-MM- end of alterations
               ENDDO
            ENDDO
        ENDIF
        ENDDO
      ENDIF

      DEALLOCATE(CWORK1)

      RETURN
      END SUBROUTINE
