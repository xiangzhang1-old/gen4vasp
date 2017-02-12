!#define debug
#include "symbol.inc"
      MODULE nonl
      
      USE prec
      INCLUDE "nonl.inc"
      CONTAINS

!****************** subroutine NONL_ALLOC  *****************************
! RCS:  $Id: nonl.F,v 1.2 2002/08/14 13:59:41 kresse Exp $
!
! allocate required arrays
! base on T_INFO and P structure
! i.e.
! number of ions and types is taken from T_INFO
! LMDIM  is taken from P structure
!***********************************************************************

      SUBROUTINE  NONL_ALLOC(NONL_S,T_INFO,P,WDES, LREAL)
      USE prec
      USE pseudo
      USE poscar
      USE wave
      IMPLICIT NONE


      TYPE (nonl_struct) NONL_S
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (wavedes)     WDES
      LOGICAL LREAL

! local var
      INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT

      NIONS =T_INFO%NIONS
      NTYPD =T_INFO%NTYPD
      LMDIM =P(1)%LMDIM
      NRPLWV=WDES%NGDIM
      NKPTS =WDES%NKPTS

      NONL_S%LRECIP=.NOT. LREAL
      NONL_S%NTYP  =T_INFO%NTYP
      NONL_S%NK    =0
      NONL_S%NIONS =T_INFO%NIONS
      NONL_S%NITYP =>T_INFO%NITYP
      NONL_S%POSION=>T_INFO%POSION
!-MM- spin spiral stuff
      NONL_S%LSPIRAL=WDES%LSPIRAL
!-MM- end of addition

      ALLOCATE(NONL_S%LMMAX(NONL_S%NTYP))
      DO NT=1,T_INFO%NTYP
        NONL_S%LMMAX(NT)=P(NT)%LMMAX
      ENDDO
!-MM- Original allocation statement
!     IF (.NOT. LREAL) &
!     ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
!              NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS), &
!              NONL_S%CQFAK(LMDIM,NTYPD))
      IF ((.NOT.LREAL).AND.(.NOT.NONL_S%LSPIRAL)) &
      ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
               NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
               NONL_S%CQFAK(LMDIM,NTYPD))
      IF ((.NOT.LREAL).AND. NONL_S%LSPIRAL) &
      ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
               NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
               NONL_S%CQFAK(LMDIM,NTYPD))
!-MM- end of alteration
      RETURN
      END SUBROUTINE

!****************** subroutine NONL_ALLOC_SPHPRO ***********************
!
! allocate required arrays to describe projectors for
! one ion for one k-point
!
!***********************************************************************

      SUBROUTINE  NONL_ALLOC_SPHPRO(NONL_S,P,WDES)
      USE prec
      USE pseudo
      USE poscar
      USE wave
      IMPLICIT NONE


      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P
      TYPE (wavedes)     WDES
! local var
      INTEGER NIONS,NTYPD,LMDIM,NRPLWV,NKPTS,NT

      NIONS =1
      NTYPD =1
      LMDIM =P%LMDIM
      NRPLWV=WDES%NGDIM
      NKPTS =WDES%NKPTS

      NONL_S%LRECIP=.TRUE.
      NONL_S%NTYP  =NTYPD
      NONL_S%NK    =0
      NONL_S%NIONS =NIONS
!-MM- spin spiral stuff
      NONL_S%LSPIRAL=WDES%LSPIRAL
!-MM- end of addition
      ALLOCATE(NONL_S%NITYP(NTYPD)); NONL_S%NITYP(1)=1
      ALLOCATE(NONL_S%LMMAX(NTYPD)); NONL_S%LMMAX(1)=P%LMMAX
!-MM- Original allocation statement
!     ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
!              NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS), &
!              NONL_S%CQFAK(LMDIM,NTYPD))
      IF (.NOT.NONL_S%LSPIRAL) &
      ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
               NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,1), &
               NONL_S%CQFAK(LMDIM,NTYPD))
      IF (NONL_S%LSPIRAL) &
      ALLOCATE(NONL_S%CREXP(NRPLWV,NIONS), &
               NONL_S%QPROJ(NRPLWV,LMDIM,NTYPD,NKPTS,2), &
               NONL_S%CQFAK(LMDIM,NTYPD))
!-MM- end of alteration
      RETURN
      END SUBROUTINE

!****************** subroutine NONL_DEALLOC_SPHPRO *********************
!
! allocate required arrays to describe projectors for
! one ion for one k-point
!
!***********************************************************************

      SUBROUTINE  NONL_DEALLOC_SPHPRO(NONL_S)
      USE prec
      IMPLICIT NONE

      TYPE (nonl_struct) NONL_S
      DEALLOCATE(NONL_S%CREXP,NONL_S%QPROJ,NONL_S%CQFAK)
      DEALLOCATE(NONL_S%NITYP)
      DEALLOCATE(NONL_S%LMMAX)

      RETURN
      END SUBROUTINE


!****************** subroutine SPHER  ********************************
!  subroutine SPHER calculates the sperical harmonics multiplied
!  by the pseudopotential QPROJ on the  compressed reciprocal
!  lattice grid
!  the full non local Pseudopotential is given by
!  <G|V(K)|GP> = SUM(L,LP,site,R)
!                     L                                  + LP
!    QPROJ(G,L site) i * DIJ(L,LP site) QPROJ(GP,LP site) i *
!                        EXP(-i(G-GP) R)
!  i^L is stored in CQFAK
!  the arrangement of the routine is, optimized for a vector-
!  processing-facility, on a scalar cpu different arrangements,
!  might result in higher efficency
!*********************************************************************

      SUBROUTINE SPHER(GRID,NONL_S,P,WDES,LATT_CUR,  IZERO,BI)
      USE prec

      USE pseudo
      USE wave
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      USE asa

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(NONL_S%NTYP)
      TYPE (wavedes)     WDES
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR

      DIMENSION BI(3,3)
! work arrays
      REAL(q),ALLOCATABLE :: GLEN(:),XS(:),YS(:),ZS(:),VPS(:),FAKTX(:), &
                             YLM(:,:)

      LYDIM=MAXL(NONL_S%NTYP,P)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      NA=WDES%NGDIM

      ALLOCATE(GLEN(NA),XS(NA),YS(NA),ZS(NA),VPS(NA),FAKTX(NA),YLM(NA,LMYDIM))

!-MM- spin spiral stuff
!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply zero when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
      QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
      QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
      QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!     QX=0;QY=0;QZ=0
!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
      IF (WDES%LSPIRAL) THEN 
         NSPINORS=2
      ELSE
         NSPINORS=1
      ENDIF
      
      spinor: DO ISPINOR=1,NSPINORS
!-MM- end of addition

!=======================================================================
! main loop over all special points
!=======================================================================
      kpoint: DO NK=1,WDES%NKPTS
!=======================================================================
! now calculate the necessary tables:
! containing the length, phasefactor exp(i m phi) and
! sin(theta),cos(theta)
!=======================================================================
      IND=1
      col: DO NC=1,GRID%RC%NCOL
      N2=GRID%RC%I2(NC) ; G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
      N3=GRID%RC%I3(NC) ; G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
      row: DO N1=1,GRID%RC%NROW

      G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
#ifdef wNGXhalf
      IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)<0) GOTO 210
      IF (GRID%LPCTX(N1)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTZ(N3)<0) GOTO 210
#endif
#ifdef wNGZhalf
      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) GOTO 210
      IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) GOTO 210
#endif
      GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
      GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
      GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI
!=======================================================================
! check to see if the kinetic energy of the plane wave is less than
! enmax in which case the plane wave is included in the set of basis
! states for the k point
!=======================================================================
      ENERGI=HSQDTM*((GIX*GIX)+(GIY*GIY)+(GIZ*GIZ))
      WSETFACT1
      WSETFACT(SQRT(2._q))

      IF(ENERGI<WDES%ENMAX) THEN
!-MM- changes to accomodate spin spirals
! original statements
!       GX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)) *TPI
!       GY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)) *TPI
!       GZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)) *TPI
        GX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)-QX) *TPI
        GY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)-QY) *TPI
        GZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)-QZ) *TPI
!-MM- end of alteration
        GLEN(IND)=MAX(SQRT(GX*GX+GY*GY+GZ*GZ),1E-10_q)
        FAKTX(IND)=WMULFACT 1
        XS(IND)  =GX/GLEN(IND)
        YS(IND)  =GY/GLEN(IND)
        ZS(IND)  =GZ/GLEN(IND)
        IND=IND+1
      ENDIF
  210 CONTINUE ! jump if point is not included
      ENDDO row
      ENDDO col
      INDMAX=IND-1
!=======================================================================
!  compare INDMAX with
!=======================================================================
      IF (INDMAX/=(WDES%NGVECTOR(NK))) THEN
        WRITE(*,*)'internal ERROR: SPHER:  INDMAX != NPLWKP'
        WRITE(*,*)'possibly the WAVECAR file can not be read'
        STOP
      ENDIF
!=======================================================================
! now calculate the tables containing the spherical harmonics
! multiply by the pseudopotential and 1/(OMEGA)^(1/2)
!=======================================================================
      CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)

      typ: DO NT=1,NONL_S%NTYP

      LMIND=1
      l_loop: DO L=1,P(NT)%LMAX
!=======================================================================
! first interpolate the non-local pseudopotentials
! and multiply by 1/(OMEGA)^(1/2)
!=======================================================================
      FAKT= 1/SQRT(LATT_CUR%OMEGA)
      ARGSC=NPSNL/P(NT)%PSMAXN

!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        ARG=(GLEN(IND)*ARGSC)+1
        NADDR=INT(ARG)

        IF (NADDR<NPSNL-2) THEN
          REM=MOD(ARG,1.0_q)
          V1=P(NT)%PSPNL(NADDR-1,L)
          V2=P(NT)%PSPNL(NADDR,L  )
          V3=P(NT)%PSPNL(NADDR+1,L)
          V4=P(NT)%PSPNL(NADDR+2,L)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          VPS(IND)=(T0+REM*(T1+REM*(T2+REM*T3)))*FAKT*FAKTX(IND)
        ELSE
          VPS(IND)=0
        ENDIF
      ENDDO
!=======================================================================
! initialize to 0
!=======================================================================
      LL=P(NT)%LPS(L)
      MMAX=2*LL

      IF (IZERO==1 .OR.IZERO==-1) THEN
      DO LM=0,MMAX
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
!-MM- changes to accomodate spin spirals
! original statement
!        NONL_S%QPROJ(IND,LMIND+LM,NT,NK)=0
         NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)=0
!-MM- end of alteration
      ENDDO; ENDDO
      ENDIF
      IF (IZERO==-1) THEN
      DO IND=1,INDMAX
         VPS(IND)=-VPS(IND)
      ENDDO
      ENDIF

!=======================================================================
! now multiply with the spherical harmonics
! and set "phase factor" CSET
!=======================================================================
      IF (LL==0) THEN
	CSET=1.0_q
      ELSE IF (LL==1) THEN
        CSET=(0.0_q,1.0_q)
      ELSE IF (LL==2) THEN
        CSET=-1.0_q
      ELSE IF (LL==3) THEN
        CSET=(0.0_q,-1.0_q)
      ENDIF

      DO LM=0,MMAX
        NONL_S%CQFAK(LMIND+LM,NT)=CSET
      ENDDO
      LMBASE=LL**2+1

      DO LM=0,MMAX
      DO IND=1,INDMAX
!-MM- alteration to accomodate spin spirals
! original statement
!       NONL_S%QPROJ(IND,LMIND+LM,NT,NK)= NONL_S%QPROJ(IND,LMIND+LM,NT,NK)+ &
!                                         VPS(IND)*YLM(IND,LM+LMBASE)
        NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)= NONL_S%QPROJ(IND,LMIND+LM,NT,NK,ISPINOR)+ &
                                          VPS(IND)*YLM(IND,LM+LMBASE)
!-MM- end of alteration
      ENDDO
      ENDDO

      LMIND=LMIND+MMAX+1
!=======================================================================
! end of loop over L
!=======================================================================

      ENDDO l_loop
      IF (LMIND-1/=P(NT)%LMMAX) THEN
        WRITE(*,*)'internal ERROR: SPHER:  LMMAX is wrong',P(NT)%LMMAX,LMIND-1
        STOP
      ENDIF
      ENDDO typ

!=======================================================================
! and of loop over special-points
!=======================================================================
      ENDDO kpoint

!-MM- spin spiral stuff 
! conjugate phase alteration for spin down: -q/2 -> q/2
      QX=-QX
      QY=-QY
      QZ=-QZ
      ENDDO spinor
!-MM- end of addition      

      DEALLOCATE(GLEN,XS,YS,ZS,VPS,FAKTX,YLM)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE PHASE *****************************
!
! this subroutine calculates the phasefactor CREXP (exp(ig.r))
! for one k-point, on the compressed grid of  lattice vectors
! CREXP is only calculated if NK changes, if ions positions change
! PHASE must be called with NK=0 to force the routine to recalculate
! the phase-factor in the next call
!***********************************************************************

      SUBROUTINE PHASE(WDES,NONL_S,NK)
      USE prec
      USE constant
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (nonl_struct) NONL_S

!=======================================================================
! check if special point changed
!=======================================================================
      IF (NK==0 .OR.NK==NONL_S%NK) THEN
        NONL_S%NK=NK
        RETURN
      ENDIF
      NONL_S%NK=NK
!=======================================================================
! number of G-Vectors on the k-Point
!=======================================================================
      NPL= WDES%NGVECTOR(NK)
!=======================================================================
! set the phase-factor
!=======================================================================
      ion:  DO NI=1,NONL_S%NIONS
      GXDX=NONL_S%POSION(1,NI)
      GYDY=NONL_S%POSION(2,NI)
      GZDZ=NONL_S%POSION(3,NI)
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
          CGDR=CITPI*(WDES%IGX(M,NK)*GXDX+WDES%IGY(M,NK)*GYDY+WDES%IGZ(M,NK)*GZDZ)
          NONL_S%CREXP(M,NI)=EXP(CGDR)
      ENDDO
      ENDDO ion

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE VNLACC  ***************************
!
!  subroutine for calculating the non local contribution of
!  the Hamiltonian, using reciprocal space projection scheme
!  the result of the wavefunction projected on the projection operatores
!  must be given in CPROJ
!
!***********************************************************************

      SUBROUTINE VNLACC(NONL_S,WDES1,W1, &
     &     LMDIM,CDIJ,CQIJ,EVALUE,  CACC)
      USE prec

      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1
      OVERLAP CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
              CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CACC(WDES1%NPL)
! work arrays
      GDEF,ALLOCATABLE    :: CRESUL(:)
      ALLOCATE(CRESUL(WDES1%NPRO))

      CALL OVERL1(WDES1, LMDIM,CDIJ,CQIJ, EVALUE, W1%CPROJ(1),CRESUL(1))

      CACC=0
      CALL VNLAC0(NONL_S,WDES1,CRESUL(1),CACC(1))

      DEALLOCATE(CRESUL)
      RETURN
      END SUBROUTINE

!************************ SUBROUTINE VNLAC0  ***************************
!
! this subroutine calculates a linear combination of
! projection operatores in reciprocal space
! the result is added to  CACC
!               -----
!***********************************************************************

      SUBROUTINE VNLAC0(NONL_S,WDES1,CPROJ_LOC,CACC)
      USE prec
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CACC(WDES1%NRPLWV)
      GDEF    CPROJ_LOC(WDES1%NPRO)

! work arrays
#if defined (use_malloc)
      REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
      POINTER    (work_p, work), (tmp_p, tmp), (cproj_p, cproj)
#else
#ifdef avoidalloc
      REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
      REAL(q),ALLOCATABLE :: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE    :: CPROJ(:)
#endif
#endif

#if defined (use_cray_ptr)
      COMPLEX (q) :: ns_cqfak(*), ns_crexp(*)
      POINTER        (ns_cqfak_p, ns_cqfak), (ns_crexp_p, ns_crexp)
#endif

#if defined (use_malloc)
      work_p  = malloc(q*wdes1%ngvector*2)
      tmp_p   = malloc(q*101*2)
      cproj_p = malloc(2*q*wdes1%npro_tot)
#else
#ifndef avoidalloc
      ALLOCATE(WORK(WDES1%NGVECTOR*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif
#endif
      NPL=WDES1%NGVECTOR

! merge projected wavefunctions from all nodes
      CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))
!=======================================================================
! performe loops over ions
!=======================================================================
      LMBASE= 0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NIS=1
      typ: DO NT=1,NONL_S%NTYP
      LMMAXC=NONL_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

#if defined (use_cray_ptr)
      ns_cqfak_p = LOC (NONL_S%CQFAK(1,NT))
#endif
      ion: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
!=======================================================================
! set TMP
!=======================================================================
#if defined (use_cray_ptr)
      ns_crexp_p = LOC (nonl_s%crexp(1,ni))
#endif
       DO L=1,LMMAXC
#if defined (use_cray_ptr)
       CTMP= CPROJ(LMBASE+L)*CONJG(ns_cqfak(L))
#else
       CTMP= CPROJ(LMBASE+L)*CONJG(NONL_S%CQFAK(L,NT))
#endif
       TMP(L,1)= REAL( CTMP ,KIND=q)
       TMP(L,2)= AIMAG(CTMP)
      ENDDO
!=======================================================================
! initialise accelerations to 0 (real imaginary part)
!=======================================================================
#ifdef do_loops
!DIR$ IVDEP
!OCL NOVREC
      DO K=1,NPL
        WORK(K)=0
        WORK(K+NPL)=0
      ENDDO
!=======================================================================
! call BLAS to calculate SUM(LM=1,LMMAX) QPROJ(K,LM) * TMP(LM)
!=======================================================================
      DO LM=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
      DO K=1,NPL
!-MM- changes to accomodate spin spirals
! original statements
!       WORK(K)    =WORK(K)    +TMP(LM,1)*NONL_S%QPROJ(K,LM,NT,WDES1%NK)
!       WORK(K+NPL)=WORK(K+NPL)+TMP(LM,2)*NONL_S%QPROJ(K,LM,NT,WDES1%NK)
        WORK(K)    =WORK(K)    +TMP(LM,1)*NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
        WORK(K+NPL)=WORK(K+NPL)+TMP(LM,2)*NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
!-MM- ens of alterations
      ENDDO; ENDDO
#else
!-MM- changes to accomodate spin spirals
! original statements
!     CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK), &
!                  WDES1%NGDIM, TMP(1,1) , 1 , 0._q , WORK(1), 1)
!     CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK), &
!                  WDES1%NGDIM, TMP(1,2) , 1 , 0._q , WORK(1+NPL), 1)
      CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, TMP(1,1) , 1 , 0._q , WORK(1), 1)
      CALL DGEMV( 'N' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, TMP(1,2) , 1 , 0._q , WORK(1+NPL), 1)
!-MM- end of alterations
#endif
!=======================================================================
! add acceleration from this ion
!=======================================================================
!DIR$ IVDEP
!OCL NOVREC
      DO K=1,NPL
         KK=K+NPL*ISPINOR
#if defined (use_cray_ptr)
         CACC(KK)=CACC(KK)+ CMPLX( WORK(K) , WORK(K+NPL) ,KIND=q) *CONJG(ns_crexp(K))
#else
         CACC(KK)=CACC(KK)+ CMPLX( WORK(K) , WORK(K+NPL) ,KIND=q) *CONJG(NONL_S%CREXP(K,NI))
#endif
      ENDDO

      LMBASE= LMMAXC+LMBASE
      ENDDO ion
  600 NIS = NIS+NONL_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONL_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition
      ENDDO spinor

#if defined (use_malloc)
      CALL free (work_p)
      CALL free (tmp_p)
      CALL free (cproj_p)
#else
#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
#endif
      RETURN
      END SUBROUTINE

!************************ SUBROUTINE PROJ    ***************************
!
! this subroutine calculates the projection of all bands of one
! specific k-point onto the reciprocal space projection operator
!
!***********************************************************************

      SUBROUTINE PROJ(NONL_S,WDES,W,NK)
      USE prec
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1


      CALL SETWDES(WDES,WDES1,NK)

      DO ISP=1,WDES%ISPIN
      DO N=1,WDES%NBANDS
        CALL SETWAV_(W,W1,N,NK,ISP)
        CALL PROJ1(NONL_S,WDES1,W1)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE PROJ1   ***************************
!
! this subroutine calculates the scalar product of one wavefunction with
! all projector functions in reciprocal space
! thesis gK Equ. (10.34)
!
!***********************************************************************

      SUBROUTINE PROJ1(NONL_S,WDES1,W1)
      USE prec
      USE wave
      
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (wavedes1)    WDES1
      TYPE (wavefun1)    W1

! work array
#if defined (use_malloc)
      REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
      POINTER    (work_p, work), (tmp_p, tmp), (cproj_p, CPROJ)
#else
#ifdef avoidalloc
      REAL(q) :: WORK(WDES1%NGVECTOR*2),TMP(101,2)
      GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
      REAL(q),ALLOCATABLE :: WORK(:),TMP(:,:)
      GDEF,ALLOCATABLE    :: CPROJ(:)
#endif
#endif

#if defined (use_cray_ptr)
      COMPLEX (q) :: ns_cqfak(*), ns_crexp(*), w1_cptwfp(*)
      POINTER        (ns_cqfak_p, ns_cqfak), (ns_crexp_p, ns_crexp), &
                     (w1_cptwfp_p, w1_cptwfp)
#endif
 

#if defined (use_malloc)
      work_p  = malloc(q*wdes1%ngvector*2)
      tmp_p   = malloc(q*101*2)
      cproj_p = malloc(2*q*wdes1%npro_tot)
#else
#ifndef avoidalloc
      ALLOCATE(WORK(WDES1%NGVECTOR*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif
#endif
      NPL=WDES1%NGVECTOR

!=======================================================================
! performe loops over ions
!=======================================================================
      LMBASE= 0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

#if defined (use_cray_ptr)
      w1_cptwfp_p = LOC (w1%cptwfp)
#endif

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NIS=1

#if defined (use_cray_ptr)
      k1  = npl * ispinor
#endif
      typ: DO NT=1,NONL_S%NTYP
      LMMAXC=NONL_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

#if defined (use_cray_ptr)
      ns_cqfak_p = LOC (NONL_S%CQFAK(1,NT))
#endif
      ion: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
#if defined (use_cray_ptr)
      ns_crexp_p = LOC (nonl_s%crexp(1,ni))
#endif
!DIR$ IVDEP
!OCL NOVREC
      DO K=1,NPL                 
#if defined (use_cray_ptr)
        CTMP= ns_crexp(K) * w1_cptwfp(k+k1)
#else
        KK=K+NPL*ISPINOR
        CTMP=    NONL_S%CREXP(K,NI) *W1%CPTWFP(KK)
#endif
        WORK(K)    = REAL( CTMP ,KIND=q)
        WORK(K+NPL)= AIMAG(CTMP)
      ENDDO
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
      l_loop: DO LM=1,LMMAXC
      SUMR=0
      SUMI=0
!DIR$ IVDEP
!OCL NOVREC
      DO K=1,NPL
!-MM- changes to accomodate spin spirals
! original statements
!       SUMR = SUMR + WORK(K)     * NONL_S%QPROJ(K,LM,NT,WDES1%NK)
!       SUMI = SUMI + WORK(K+NPL) * NONL_S%QPROJ(K,LM,NT,WDES1%NK)
        SUMR = SUMR + WORK(K)     * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
        SUMI = SUMI + WORK(K+NPL) * NONL_S%QPROJ(K,LM,NT,WDES1%NK,ISPIRAL)
!-MM- end of alterations
      ENDDO
      CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *NONL_S%CQFAK(LM,NT))
      ENDDO l_loop
#else
!-MM- changes to accomodate spin spirals
! original statements
!     CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK), &
!                  WDES1%NGDIM, WORK(1) , 1 , 0._q ,  TMP(1,1), 1)
!     CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK), &
!                  WDES1%NGDIM, WORK(1+NPL) , 1 , 0._q ,  TMP(1,2), 1)
      CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, WORK(1) , 1 , 0._q ,  TMP(1,1), 1)
      CALL DGEMV( 'T' , NPL, LMMAXC, 1._q , NONL_S%QPROJ(1,1,NT,WDES1%NK,ISPIRAL), &
                   WDES1%NGDIM, WORK(1+NPL) , 1 , 0._q ,  TMP(1,2), 1)
!-MM- end of alterations
      l_loop: DO LM=1,LMMAXC
        SUMR=TMP(LM,1)
        SUMI=TMP(LM,2)
#if defined (use_cray_ptr)
        CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) * ns_cqfak(lm))
#else
        CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *NONL_S%CQFAK(LM,NT))
#endif
      ENDDO l_loop
#endif

      LMBASE=LMBASE+LMMAXC
      ENDDO ion

  600 NIS = NIS+NONL_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONL_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition       
      ENDDO spinor

! distribute the projected wavefunctions to nodes
      CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#if defined (use_malloc)
      CALL free (work_p)
      CALL free (tmp_p)
      CALL free (cproj_p)
#else
#ifndef avoidalloc
      DEALLOCATE(WORK,TMP,CPROJ)
#endif
#endif
      RETURN
      END SUBROUTINE

      END MODULE


!************************ SUBROUTINE FORNL  ****************************
!
! this subroutine calculates the forces related to the non local
! pseudopotential and one plane-wave basis-state C
! the projection of the wavefunction ontot the projection operatores
! must be stored in CPROJ
! Algorithm:
! ACC(G') = SUM(L,L',R)
!         SUM(G)   QPROJ(G ,L)  EXP( iG  R) C(G)  iG *
!         D(L,L') + Evalue  Q(L,L')
!         SUM(GP)  QPROJ(GP,LP) EXP(-iGP R) C(GP)
!
! CHEERS to CRAY this routine does not compile in the MODULE
!***********************************************************************

      SUBROUTINE FORNL(NONL_S,WDES,W,LATT_CUR, LMDIM,NIOND,CDIJ,CQIJ,EINL)
      USE prec
      USE nonl
      USE wave
      USE lattice
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (latt)        LATT_CUR
      TYPE (wavedes1)    WDES1         ! descriptor for one k-point

      DIMENSION EINL(3,NONL_S%NIONS)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! work arrays
      REAL(q) ::  ENL(NONL_S%NIONS)
      GDEF :: CM
      GDEF,ALLOCATABLE :: CX(:) ,CY(:) ,CZ(:)
      GDEF,ALLOCATABLE :: CXL(:),CYL(:),CZL(:)

      N =WDES%NPRO_TOT
      NL=WDES%NPRO

      ALLOCATE(CX(N),CY(N),CZ(N),CXL(NL),CYL(NL),CZL(NL))

      EINL=0._q
      ENL =0._q
!=======================================================================
! loop over special points, and bands
!=======================================================================

      kpoint: DO NK=1,WDES%NKPTS
      NPL= WDES%NGVECTOR(NK)
      CALL SETWDES(WDES,WDES1,NK)

      CALL PHASE(WDES,NONL_S,NK)

      spin:   DO ISP=1,WDES%ISPIN
      band: DO N=1,WDES%NBANDS

      EVALUE=W%CELEN(N,NK,ISP)
      WEIGHT=W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)*WDES%RSPIN
!=======================================================================
! first build up CX, CY, CZ for tables
!=======================================================================
      CX=0
      CY=0
      CZ=0

      LMBASE= 0
      
!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition

      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

      NIS=1
      typ:  DO NT=1,NONL_S%NTYP
      LMMAXC=NONL_S%LMMAX(NT)
      IF (LMMAXC==0) GOTO 100

      ions: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
      l_loop: DO LM=1,LMMAXC
      CMUL=NONL_S%CQFAK(LM,NT)*CITPI
!DIR$ IVDEP
!OCL NOVREC

      DO K=1,NPL                 
         KK=K+NPL*ISPINOR
!-MM- changes to accomodate spin spirals
! original statement
!        CVAL =  NONL_S%QPROJ(K,LM,NT,NK)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP)*CMUL
         CVAL =  NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W%CPTWFP(KK,N,NK,ISP)*CMUL
!-MM- end of changes
         CX(LM+LMBASE)=CX(LM+LMBASE)+WDES%IGX(K,NK)*CVAL
         CY(LM+LMBASE)=CY(LM+LMBASE)+WDES%IGY(K,NK)*CVAL
         CZ(LM+LMBASE)=CZ(LM+LMBASE)+WDES%IGZ(K,NK)*CVAL
      ENDDO

      ENDDO l_loop
      LMBASE= LMMAXC+LMBASE
      ENDDO ions

  100 NIS = NIS+NONL_S%NITYP(NT)
      ENDDO typ
!-MM- spin spiral stuff
      IF (NONL_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition
      ENDDO spinor

      CALL DIS_PROJ(WDES1,CX,CXL)
      CALL DIS_PROJ(WDES1,CY,CYL)
      CALL DIS_PROJ(WDES1,CZ,CZL)
!=======================================================================
! sum up local contributions
! calculate SUM_LP  D(LP,L)-E Q(LP,L) * C(LP)
!=======================================================================
      spinor2 : DO ISPINOR=0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *WDES%NPRO/2
      LMBASE_=ISPINOR_*WDES%NPRO/2

      NIS   =1
      typ2:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 600

      ions2: DO NI=NIS,WDES%NITYP(NT)+NIS-1
      NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index
      l_loop2: DO LM=1,LMMAXC
         CM=0
         DO LMP=1,LMMAXC

            CM=     (CDIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
            GCONJG(W%CPROJ(LM+LMBASE,N,NK,ISP))
            ENL(NIP)=ENL(NIP) + (W%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
            EINL(1,NIP)=EINL(1,NIP)-(2*WEIGHT)*(CXL(LMP+LMBASE_)*CM) 
            EINL(2,NIP)=EINL(2,NIP)-(2*WEIGHT)*(CYL(LMP+LMBASE_)*CM)
            EINL(3,NIP)=EINL(3,NIP)-(2*WEIGHT)*(CZL(LMP+LMBASE_)*CM)

         ENDDO
      ENDDO l_loop2
      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_
      ENDDO ions2
  600 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ2
      ENDDO
      ENDDO spinor2
!=======================================================================
      ENDDO band
      ENDDO spin
      ENDDO kpoint
!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
      CALLMPI( M_sum_d(WDES%COMM, EINL(1,1),NONL_S%NIONS*3))
      CALLMPI( M_sum_d(WDES%COMM, ENL (1)  ,NONL_S%NIONS))

      CALL  DIRKAR(NONL_S%NIONS,EINL,LATT_CUR%B)

      DEALLOCATE(CX,CY,CZ,CXL,CYL,CZL)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE STRENL  ****************************
!
! calculate non-local contributions to stress,
! easiest to implement and definitly quite fast, is an approach based
! on finite differences, the implementation took only 3 hours (and this is
! -I think- the most compelling feature of the routine)
! CHEERS to CRAY this routine does not compile in the MODULE
!
!***********************************************************************

      SUBROUTINE STRENL(GRID,NONL_S,P,W,WDES,LATT_CUR,  BI, &
          LMDIM,NIOND,CDIJ,CQIJ, ISIF,FNLSIF)
      USE prec
      USE nonl
      USE pseudo
      USE wave
      USE mpimy
      USE mgrid
      USE lattice
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(NONL_S%NTYP)
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (wavespin)    W
      TYPE (wavefun1)    W1
      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_FIN

      DIMENSION FNLSIF(3,3)
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      DIMENSION BI(3,3)
! work arrys
      GDEF,ALLOCATABLE,TARGET ::  CWORK(:)

      DIS=1E-5_q
!TEST which precission should be used in the DIS statment above
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
      ALLOCATE(CWORK(WDES%NPROD))

      FNLSIF=0

      DO IDIR=1,3
      DO JDIR=1,3
!=======================================================================
! use central differences to calculate the stress
! set up QPROJ so that central differences can be  calculated
!=======================================================================
      LATT_FIN=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3; DO J=1,3
          LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
        ENDDO; ENDDO
      ELSE
!  all directions
        DO I=1,3
          LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
        ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN)

      IZERO=-1
      CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO,BI)

      LATT_FIN=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3; DO J=1,3
          LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
        ENDDO; ENDDO
      ELSE 
!  all directions
        DO I=1,3
          LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
        ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN)

      IZERO=0
      CALL SPHER(GRID,NONL_S,P,WDES,LATT_FIN,  IZERO,BI)

!=======================================================================
! loop over all k-points spin and bands
!=======================================================================
  kpoint: DO NK=1,WDES%NKPTS

      CALL PHASE(WDES,NONL_S,NK)
      CALL SETWDES(WDES,WDES1,NK)

    spin: DO ISP=1,WDES%ISPIN
    band: DO N=1,WDES%NBANDS
      CALL SETWAV_(W,W1,N,NK,ISP); W1%CPROJ => CWORK
      CALL PROJ1(NONL_S,WDES1,W1)  ! calculate W1%CPROJ (linked to CWORK)

      WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*W1%FERWE
      EVALUE=W1%CELEN

      LBASE= 0
      spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LBASE =ISPINOR *WDES%NPRO/2
      LBASE_=ISPINOR_*WDES%NPRO/2

      NIS=1

      DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 270

      DO NI=NIS,WDES%NITYP(NT)+NIS-1
      DO L=1 ,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
      DO LP=1,LMMAXC
        FNLSIF(IDIR,JDIR) = FNLSIF(IDIR,JDIR)+(2*WEIGHT)*CWORK(LBASE_+LP)* &
                            GCONJG(W%CPROJ(LBASE+L,N,NK,ISP))* &
                           (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
      ENDDO
      ENDDO

      LBASE=  LMMAXC+LBASE
      LBASE_= LMMAXC+LBASE_
      ENDDO
  270 NIS = NIS+WDES%NITYP(NT)
      ENDDO
      ENDDO
      ENDDO spinor
      ENDDO band
      ENDDO spin
      ENDDO kpoint
!
!  only isotrop pressure finish now
!
      IF (ISIF==1) THEN
        FNLSIF(2,2)= FNLSIF(1,1)
        FNLSIF(3,3)= FNLSIF(1,1)
        GOTO 310  ! terminate (not very clean but who cares)
      ENDIF
!=======================================================================
! next direction
!=======================================================================
      ENDDO
      ENDDO

  310 CONTINUE
      CALLMPI( M_sum_d(WDES%COMM, FNLSIF, 9))

      FNLSIF=FNLSIF/DIS/2

#ifdef debug
      WRITE(*,'(E10.3,3E14.7)')DIS,((FNLSIF(I,J),I=1,3),J=1,3)
!      IF (DIS>1E-10) GOTO 1000
!      GOTO 1000
#endif
!=======================================================================
! recalculate the projection operators
! (the array was used as a workspace)
!=======================================================================
      IZERO=1
      CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR,  IZERO,BI)

      DEALLOCATE(CWORK)
      RETURN
      END SUBROUTINE
