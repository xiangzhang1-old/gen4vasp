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





      MODULE dfast
      USE prec
      CONTAINS
!***********************************************************************
! RCS:  $Id: dfast.F,v 1.2 2001/02/20 14:44:57 kresse Exp $
!
! This modul contains the parts of the program which are important for
! the performance (hot spots).
! routines in this file call BLAS subroutines
! compiled with an auto-double compiler
! following BLAS routines are used:
!      ZGEMM, ZGEMV, ZDOTC, ZAXPY, ZDSCAL
! MIND: there are also some hot spots in
!   nonl.F and nonlr.F (at least as important)
!   but I decided to keep them where they are now
!***********************************************************************

!************************** SUBROUTINE PROALL***************************
!
! this subroutine
! calculates the scalar product of the current wavefunctions
! stored in W with the projectors
!
!  C_lme ion,n = < b_lme ion | psi_n >
!
! and stores the result
!
!***********************************************************************

      SUBROUTINE PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,WDES,W,LOVERL,LREAL,LMDIM)
      USE prec

      USE mpimy
      USE mgrid
      USE nonl
      USE nonlr
      USE poscar
      USE wave
      USE lattice
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W

      LOGICAL LREAL,LOVERL
!=======================================================================
! calculate the projections
! loop over all special points
!=======================================================================
      DO NK=1,WDES%NKPTS
        IF (LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES,0.0_q,0.0_q,0.0_q)
        CALL RPRO(NONLR_S,WDES,W,GRID,NK)

        ELSE
        CALL PHASE(WDES,NONL_S,NK)
        CALL PROJ(NONL_S,WDES,W,NK)
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE

!************************* SUBROUTINE ORTHON ***************************
!
! orthogonalize a wavefunction W1 to all other bands
! including the current band
! the subroutine uses BLAS 3 calls,
!***********************************************************************

      SUBROUTINE ORTHON(WDES,NK,W,W1, LOVERL,LMDIM,CQIJ,ISP)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)


      TYPE (wavedes)    WDES
      TYPE (wavespin)   W
      TYPE (wavefun1)   W1

      LOGICAL LOVERL
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NRSPINORS*WDES%NRSPINORS)
! work arrays
      COMPLEX(q),ALLOCATABLE:: CPRO(:),CWORK(:)

      ALLOCATE(CPRO(WDES%NBANDS),CWORK(WDES%NPRO))
!=======================================================================
! caclulate P |w>
!=======================================================================
      IF (LOVERL) THEN
         CWORK=0

         NPRO=0
         spinor: DO ISPINOR=0,WDES%NRSPINORS-1
         DO ISPINOR_=0,WDES%NRSPINORS-1

         NPRO =ISPINOR *WDES%NPRO/2
         NPRO_=ISPINOR_*WDES%NPRO/2

         NIS =1
         DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 230
         DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
            DO LP=1,LMMAXC
               CWORK(L+NPRO)=CWORK(L+NPRO)+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*W1%CPROJ(LP+NPRO_)
            ENDDO
            ENDDO

            NPRO = LMMAXC+NPRO
            NPRO_= LMMAXC+NPRO_
         ENDDO

 230     NIS = NIS+WDES%NITYP(NT)
         ENDDO
         ENDDO
         ENDDO spinor
      ENDIF
!=======================================================================
! performe orthogonalisation
!=======================================================================

     CALL ZGEMV( 'C' ,  WDES%NPLWKP(NK) , WDES%NBANDS ,(1._q,0._q) , W%CPTWFP(1,1,NK,ISP), &
     &              WDES%NRPLWV, W1%CPTWFP(1) , 1 , (0._q,0._q) ,  CPRO(1), 1)

      IF (LOVERL ) THEN
       IF (WDES%NPRO /= 0) &
       CALL ZGEMV( 'C' ,  WDES%NPRO , WDES%NBANDS ,(1._q,0._q) , W%CPROJ(1,1,NK,ISP) , &
                    WDES%NPROD, CWORK(1), 1 , (1._q,0._q) ,  CPRO(1), 1)
      ENDIF
      

      CALL ZGEMM( 'N', 'N' ,  WDES%NPLWKP(NK) , 1 , WDES%NBANDS , -(1._q,0._q) , &
                   W%CPTWFP(1,1,NK,ISP),  WDES%NRPLWV , CPRO(1) , WDES%NBANDS , &
                   (1._q,0._q) , W1%CPTWFP(1) ,  WDES%NRPLWV )

      IF (WDES%NPRO /= 0) &
      CALL ZGEMM( 'N', 'N' ,  WDES%NPRO , 1 , WDES%NBANDS  , -(1._q,0._q) , &
                   W%CPROJ(1,1,NK,ISP) ,  WDES%NPROD , CPRO(1) , WDES%NBANDS , &
                   (1._q,0._q) , W1%CPROJ(1) ,  WDES%NPROD  )

      DEALLOCATE(CPRO,CWORK)

      RETURN
      END SUBROUTINE

!************************* SUBROUTINE CNORMN  **************************
!
! this subroutine normalises a wavefunction
! subroutine is not important for performance
!***********************************************************************

      SUBROUTINE CNORMN(WDES1,W, LOVERL,LMDIM,CQIJ,WSCAL)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W

      COMPLEX(q) ZDOTC
      LOGICAL LOVERL
      COMPLEX(q)      CP
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)

      WFMAG=ZDOTC(WDES1%NPL,W%CPTWFP(1),1,W%CPTWFP(1),1)
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
      IF (LOVERL) THEN
         CP  =0
         NPRO=0
         spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
         DO ISPINOR_=0,WDES1%NRSPINORS-1

         NPRO =ISPINOR *WDES1%NPRO/2
         NPRO_=ISPINOR_*WDES1%NPRO/2

         NIS =1
         DO NT=1,WDES1%NTYP
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC==0) GOTO 230

         DO NI=NIS,WDES1%NITYP(NT)+NIS-1
            DO L=1,LMMAXC
            DO LP=1,LMMAXC
               CP=CP+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*W%CPROJ(LP+NPRO_)*CONJG(W%CPROJ(L+NPRO))
            ENDDO
            ENDDO

            NPRO = LMMAXC+NPRO
            NPRO_= LMMAXC+NPRO_
         ENDDO
 230     NIS = NIS+WDES1%NITYP(NT)
         ENDDO
         ENDDO
         ENDDO spinor

         WFMAG=WFMAG+CP
      ENDIF

      

!-----check that it is non-(0._q,0._q)
      IF(WFMAG<=0) THEN
!=======================================================================
! if it is (0._q,0._q) halt the execution
!=======================================================================
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
        WSCAL=WFMAG
      ELSE
        WSCAL= 1._q/SQRT(WFMAG)
        CALL ZDSCAL( WDES1%NPL ,WSCAL,W%CPTWFP(1),1)
        CALL ZDSCAL( WDES1%NPRO,WSCAL,W%CPROJ(1),1)
      ENDIF
      END SUBROUTINE

!************************* SUBROUTINE CNORMA  **************************
!
! this subroutine calculates the norm of a wavefunction
!
!***********************************************************************

      SUBROUTINE CNORMA(WDES1,W, LOVERL,LMDIM,CQIJ,WSCAL)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W

      COMPLEX(q) ZDOTC
      LOGICAL LOVERL
      COMPLEX(q)      CP
      COMPLEX(q)   CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)

      WFMAG=ZDOTC(WDES1%NPL,W%CPTWFP(1),1,W%CPTWFP(1),1)
      CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
      IF (LOVERL) THEN
         NPRO=0

         spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
         DO ISPINOR_=0,WDES1%NRSPINORS-1

         NPRO =ISPINOR *WDES1%NPRO/2
         NPRO_=ISPINOR_*WDES1%NPRO/2

         NIS =1
         DO NT=1,WDES1%NTYP
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC==0) GOTO 230

         DO NI=NIS,WDES1%NITYP(NT)+NIS-1
            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CP=CP+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*W%CPROJ(LP+NPRO_)*CONJG(W%CPROJ(L+NPRO))
               ENDDO
            ENDDO

            NPRO = LMMAXC+NPRO
            NPRO_= LMMAXC+NPRO_
         ENDDO
 230     NIS = NIS+WDES1%NITYP(NT)
         ENDDO
         ENDDO
         ENDDO spinor

         WFMAG=WFMAG+CP
      ENDIF

      
!-----check that it is non-(0._q,0._q)
      IF(WFMAG<=0) THEN
!=======================================================================
! if it is (0._q,0._q) halt the execution
!=======================================================================
        WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
        WSCAL=WFMAG
      ELSE
        WSCAL= 1._q/SQRT(WFMAG)
      ENDIF
      END SUBROUTINE


!************************* SUBROUTINE PROJCN ***************************
!
! this subroutine projects out from (1._q,0._q) wavefunction
! CF another wavefunction  CPRO
! subroutine is not important for performance
!
!***********************************************************************

      SUBROUTINE PROJCN(WDES1,W1,W2, LOVERL,LMDIM,CQIJ, CSCPD)
      USE prec

      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1,W2
      COMPLEX(q)    CADD

      COMPLEX(q) ZDOTC
      LOGICAL LOVERL
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)

      CSCPD= (ZDOTC(WDES1%NPL,W2%CPTWFP(1),1,W1%CPTWFP(1),1))
!=======================================================================
! if necessary caclulate <p| P |w>
!=======================================================================
      IF (LOVERL) THEN
         CADD=0
         NPRO=0
         spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
         DO ISPINOR_=0,WDES1%NRSPINORS-1

         NPRO =ISPINOR *WDES1%NPRO/2
         NPRO_=ISPINOR_*WDES1%NPRO/2

         NIS =1
         DO  NT=1,WDES1%NTYP
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC==0) GOTO 230

         DO NI=NIS,WDES1%NITYP(NT)+NIS-1
            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CADD=CADD+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*W1%CPROJ(LP+NPRO_)*CONJG(W2%CPROJ(L+NPRO))
               ENDDO
            ENDDO

            NPRO = LMMAXC+NPRO
            NPRO_= LMMAXC+NPRO_
         ENDDO
 230     NIS = NIS+WDES1%NITYP(NT)
         ENDDO
         ENDDO
         ENDDO spinor

         CSCPD=(CSCPD+CADD)
      ENDIF
!=======================================================================
! performe orthogonalisations
!=======================================================================
      

      CALL ZAXPY(WDES1%NPL ,-CSCPD,W2%CPTWFP(1)   ,1,W1%CPTWFP(1)   ,1)
      CALL ZAXPY(WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE ORSP   ****************************
!
! this subroutine perfomes a gram-schmidt orthogonalistion of a set
! of vectors (all elements on local node)
! the subroutine uses BLAS 3 calls,
!***********************************************************************

      SUBROUTINE ORSP(NBANDS,NPL,NRPLWV,CPTWFP)
      USE prec
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) CPTWFP(NRPLWV,NBANDS)
      COMPLEX(q) CPRO(NBANDS)
      COMPLEX(q) ZDOTC
      CPRO=0
!=======================================================================
! loop over all vectors
!=======================================================================
      DO N=1,NBANDS
!=======================================================================
! normalise the vector
!=======================================================================
       WFMAG=ZDOTC(NPL,CPTWFP(1,N),1,CPTWFP(1,N),1)
       CALL ZDSCAL(NPL,1/SQRT(WFMAG),CPTWFP(1,N),1)
!=======================================================================
! now orthogonalise all higher vectors to the
! present vector
!=======================================================================
       IF (NBANDS/=N ) THEN
       CALL ZGEMV( 'C', NPL , NBANDS-N ,(1._q,0._q) , CPTWFP(1,N+1), &
     &             NRPLWV, CPTWFP(1,N), 1 , (0._q,0._q) ,  CPRO, 1)

       DO I=1,NBANDS
         CPRO(I)=CONJG(CPRO(I))
       ENDDO

       CALL ZGEMM( 'N', 'T' , NPL , NBANDS-N , 1 , -(1._q,0._q) , &
     &             CPTWFP(1,N), NRPLWV , CPRO , NBANDS , &
     &             (1._q,0._q) , CPTWFP(1,N+1) , NRPLWV )
      ENDIF
      ENDDO
      RETURN
      END SUBROUTINE

      END MODULE

!************************ SUBROUTINE LINCOM ****************************
!
! build linear combinations of wavefunctions according to matrix CTRANS
! this subroutine performes implicitly a MATRIX x MATRIX multiplication,
! it is needed for the unitary transformation of the wavefunctions or
! for orthogonalisation routines and uses a blocked algorithm
! to save storage
!
!  COUT_n k =  CH_ kp k CIN_n kp
! where n  =  1 ... NPL    (leading dimension of array = NPLDIM)
!       k  =  1 ... NOUT
!       kp =  1 ... NIN    (NIN must be greater or equal NOUT)
! Important: on exit the input array will be overwritten by output data!
!
!  MODE determines the mode for the transformation
!  "U"   CTRANS   upper triangle set
!  "L"   CTRANS   lower triangle set
!  "F"   CTRANS   all components set
!
!  "A"   used only for Davidson
!  "B"   used only for Davidson
!***********************************************************************


      SUBROUTINE LINCOM(MODE,CF,CPROF,CTRANS,NIN,NOUT,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN,NBLK,CFA,CPROFA)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER*1 MODE
      DIMENSION   CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      COMPLEX(q)        CPROF(NPROD,NIN),CPROFA(NPROD,NIN)
      COMPLEX(q)        CTRANS(LDTRAN,NIN)

! work array
      COMPLEX(q),ALLOCATABLE ::   CBLOCK(:,:)
      ALLOCATE(CBLOCK(NBLK,LDTRAN))

      CALL LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT, NPL, &
     &             NPLDIM,LDTRAN, NBLK,CFA)
      IF (NPRO/=0) THEN
      CALL LINBAS(MODE,CPROF,CBLOCK,CTRANS,NIN,NOUT, NPRO, &
     &             NPROD,LDTRAN, NBLK,CPROFA)
      ENDIF
      DEALLOCATE(CBLOCK)

      RETURN
      END SUBROUTINE


!************************* SUBROUTINE OVERL ***************************
!
! calculate the result of the overlap-operator acting onto a set of
! wavefunctions
!**********************************************************************

      SUBROUTINE OVERL(WDES, LOVERL,LMDIM,CQIJ, CPROF,CRESUL)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes1) WDES
      COMPLEX(q) CRESUL(WDES%NPROD,WDES%NBANDS),CPROF(WDES%NPROD,WDES%NBANDS)

      LOGICAL LOVERL
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NRSPINORS*WDES%NRSPINORS)

      IF (LOVERL) THEN

       BANDS: DO NB=1,WDES%NBANDS

       DO NP=1,WDES%NPRO
        CRESUL(NP,NB)=0
       ENDDO

       spinor: DO ISPINOR=0,WDES%NRSPINORS-1
       DO ISPINOR_=0,WDES%NRSPINORS-1

       NPRO =ISPINOR *WDES%NPRO/2
       NPRO_=ISPINOR_*WDES%NPRO/2

       NIS =1
       DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 100
         DO NI=NIS,WDES%NITYP(NT)+NIS-1

           DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
           DO LP=1,LMMAXC
             CRESUL(L+NPRO,NB)=CRESUL(L+NPRO,NB)+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO_,NB)
           ENDDO
           ENDDO
           NPRO = LMMAXC+NPRO
           NPRO_= LMMAXC+NPRO_
         ENDDO
 100     NIS = NIS+WDES%NITYP(NT)
       ENDDO
       ENDDO
       ENDDO spinor

       ENDDO BANDS
      ENDIF

      RETURN
      END SUBROUTINE


!************************* SUBROUTINE OVERL1 **************************
!
! calculate the result of the overlap-operator acting onto (1._q,0._q)
! wavefunction
!
!**********************************************************************

      SUBROUTINE OVERL1(WDES1, LMDIM,CDIJ,CQIJ, EVALUE, CPROF,CRESUL)
      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes1) WDES1
      COMPLEX(q) CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
              CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)

      NPRO=0
      CRESUL=0

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES1%NRSPINORS-1

      NPRO =ISPINOR *WDES1%NPRO/2
      NPRO_=ISPINOR_*WDES1%NPRO/2

      NIS =1
      DO NT=1,WDES1%NTYP
        LMMAXC=WDES1%LMMAX(NT)
        IF (LMMAXC==0) GOTO 100
        DO NI=NIS,WDES1%NITYP(NT)+NIS-1

          DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
          DO LP=1,LMMAXC
            CRESUL(L+NPRO)=CRESUL(L+NPRO)+(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                                    EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
          ENDDO; ENDDO
          NPRO = LMMAXC+NPRO
          NPRO_= LMMAXC+NPRO_
        ENDDO
 100    NIS = NIS+WDES1%NITYP(NT)
      ENDDO
      ENDDO
      ENDDO spinor

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE ORTH1 *****************************
! calculates (1._q,0._q) strip of the overlap between two set of vectors
!   O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
! only lower part of O is calculated
! (i.e. resulting matrix should be hermitian)
! use ORTH2 if this is not the case
!***********************************************************************

      SUBROUTINE ORTH1(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NBLK, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)
      CHARACTER*(*) CSEL

      IF (NSTRIP+NPOS-1 > NBANDS) THEN
        WRITE(*,*)'internal error in ORTH1: dim=',NSTRIP+NPOS,NBANDS
        STOP
      ENDIF
!
! update of lower triangular part
!
    IF (CSEL(1:1) == 'L' .OR. CSEL(1:1) == 'l') THEN
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF
!
! update of upper triangular part
!
    ELSE IF (CSEL(1:1) == 'U' .OR. CSEL(1:1) == 'u') THEN
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

    ELSE
      WRITE(*,*)'internal error in ORTH1: CSEL=',CSEL
    ENDIF

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE ORTH1 *****************************
! calculates (1._q,0._q) strip of the overlap between two set of vectors
!   O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
! only lower part of O is calculated
! (i.e. resulting matrix should be hermitian)
! use ORTH2 if this is not the case
!***********************************************************************

      SUBROUTINE ORTH1_old(CPTWFP,CFW,CPROJ,CPROW,NBANDS,NBLK, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)

! Try to get best load balance, maximum block size < NBLK ...
      IF (NPL/=0) THEN
      NBLOCK= NPL
      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
     &         CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
     &         CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE ORTH2 *****************************
! calculates (1._q,0._q) strip of the overlap between two set of vectors
!   O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
! for general case (i.e. resulting matrix must not be Hermitian)
!***********************************************************************

      SUBROUTINE ORTH2(CPTWFP,CFW,CPROJ,CPROW,NBANDS,NBLK, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)

! Try to get best load balance, maximum block size < NBLK ...
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
     &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN
      NBLOCK=NPRO

      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
     &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      RETURN
      END SUBROUTINE



!************************ SUBROUTINE LINBAS ****************************
!
! build linear combinations of set of vectors according to matrix CTRANS
! this subroutine performes implicitly a MATRIX x MATRIX multiplication,
! it is needed for the unitary transformation of the wavefunctions or
! for orthogonalisation routines and uses a blocked algorithm
! to save storage
! LINBAS is only called from LINCOM
!***********************************************************************

      SUBROUTINE LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT,NPL, &
     &           NPLDIM,LDTRAN,NBLK,CFA)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

! Might be that ZGEMM performs much much faster than ZTRMM ... ??
      PARAMETER(IUSETR=1)

      CHARACTER*1 MODE
      LOGICAL     LTRI,LADD,LBOTH
      COMPLEX(q)     CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      COMPLEX(q)     CBLOCK(NBLK,LDTRAN)
      COMPLEX(q)     CTRANS(LDTRAN,NIN)

      IF (NOUT>NIN) THEN
         WRITE(*,1)
 1       FORMAT('internal error in routine LINBAS: wrong arguments, NOUT>NIN')
         STOP
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
      LADD =(MODE=='A').OR.(MODE=='a')
      LBOTH=(MODE=='B').OR.(MODE=='b')
      IF (LTRI.AND.(IUSETR==0)) THEN
         DO 4 N2=1,NIN
            IF ((MODE=='L').OR.(MODE=='l')) THEN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2-1
               CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ELSE
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+1,NIN
               CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDIF
    4    ENDDO
      ENDIF

! Try to get best load balance, maximum block size < NBLK ...
      NBLOCK=NBLK

      DO 70 IBLOCK=0,NPL-1,NBLOCK
         ILENPL=MIN(NBLOCK,NPL-IBLOCK)
         IADDPL=MIN(IBLOCK,NPL-1)
         ILENPL=MAX(ILENPL,0)

         IF (LTRI.AND.(IUSETR/=0)) THEN
! 'Triangular update':
            CALL ZTRMM &
     &                ('R',MODE,'N','N', ILENPL,NOUT,(1._q,0._q), &
     &                 CTRANS,LDTRAN,CF(1+IADDPL,1), NPLDIM)
         ELSE
! 'Full update':
            IF (LBOTH.OR.(.NOT.LADD)) THEN
               DO 30 N1=1,NIN
                  DO 10 M=1,ILENPL
                     CBLOCK(M,N1)=CF(M+IADDPL,N1)
   10             CONTINUE
   30          CONTINUE
               CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, (0._q,0._q), CF(IADDPL+1,1), NPLDIM)
            ENDIF
            IF (LBOTH.OR.LADD) THEN
               IADDT=0
               IF (LBOTH) IADDT=NIN
               DO 60 N1=1,NIN
                  DO 40 M=1,ILENPL
                     CBLOCK(M,N1)=CFA(M+IADDPL,N1)
   40             CONTINUE
   60          CONTINUE
               CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &             CBLOCK(1,1), NBLK, CTRANS(1+IADDT,1), &
     &                   LDTRAN, (1._q,0._q), CF(IADDPL+1,1), NPLDIM)
            ENDIF
         ENDIF

   70 CONTINUE

      RETURN
      END SUBROUTINE
