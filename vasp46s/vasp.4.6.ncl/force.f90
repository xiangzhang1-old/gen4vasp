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





      MODULE force
      USE prec
      CONTAINS
      SUBROUTINE mpi_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
      END SUBROUTINE
      END MODULE
!************************ SUBROUTINE STRELO ****************************
! RCS:  $Id: force.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
! this calculates the stress on the unit cell
! which is related to the change in local pseudopotential on changing
! the unit vectors (EISIF + PSCSIF)
! and the stress related to the change of the hartree energy (DSIF)
! see formulas (10.30) (10.51) and (10.52) in thesis gK
!
!***********************************************************************

      SUBROUTINE STRELO(GRIDC,P,T_INFO,LATT_CUR, &
           CHTOT,CSTRF, NELECT,DSIF,EISIF,PSCSIF)
      USE prec

      USE pot
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      COMPLEX(q) CHTOT(GRIDC%RC%NP)
      REAL(q)    EISIF(3,3),DSIF(3,3),PSCSIF(3,3),NELECT

! work arrays
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK2(:)

      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
!=======================================================================
! first caclulate the local potential and its derivative
!=======================================================================
      CALL POTION(GRIDC,P,LATT_CUR,T_INFO,CWORK1,CWORK2,CSTRF,PSCENC)
!=======================================================================
! there are two contributions to the force on the unit cell,
! (1._q,0._q) due to the tendency to put the charge density  at the largest
! values of the pseudopotential and the other due to the
! 1/(cell volume) dependence
! see formulas (10.51) and (10.52) in thesis gK
!=======================================================================
      EISIF =0    ! stress due to local potential
      DSIF  =0    ! stress due to Hartree potential
      PSCSIF=0    ! stress due to q-> 0 behaviour (10.30) in thesis gK

      PSCV=0      ! 1/ volume part for local PP
      DENC=0      ! 1/ volume part for Hartree pot

      NG=1
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW

        FACTM=1
        

        GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
        GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
        GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

        GSQU=GX**2+GY**2+GZ**2
        GMOD=SQRT(GSQU)
        PSCV=PSCV+ REAL( CONJG(CHTOT(NG))*CWORK1(NG) ,KIND=q)
!     G=0 part is handled somewhere else
        IF (GMOD>1E-7_q) THEN ! avoid division by 0

         PSCD=  REAL( CONJG(CHTOT(NG))*CWORK2(NG) ,KIND=q)
         EISIF(1,1)=EISIF(1,1)+PSCD*GX*GX/GMOD
         EISIF(2,2)=EISIF(2,2)+PSCD*GY*GY/GMOD
         EISIF(3,3)=EISIF(3,3)+PSCD*GZ*GZ/GMOD
         EISIF(1,2)=EISIF(1,2)+PSCD*GX*GY/GMOD
         EISIF(2,3)=EISIF(2,3)+PSCD*GY*GZ/GMOD
         EISIF(3,1)=EISIF(3,1)+PSCD*GZ*GX/GMOD

         DUM= CHTOT(NG)*CONJG(CHTOT(NG))/GSQU
         DENC=DENC+DUM
         DSIF(1,1)=DSIF(1,1)-DUM*GX*GX/GSQU
         DSIF(2,2)=DSIF(2,2)-DUM*GY*GY/GSQU
         DSIF(3,3)=DSIF(3,3)-DUM*GZ*GZ/GSQU
         DSIF(1,2)=DSIF(1,2)-DUM*GX*GY/GSQU
         DSIF(2,3)=DSIF(2,3)-DUM*GY*GZ/GSQU
         DSIF(3,1)=DSIF(3,1)-DUM*GZ*GX/GSQU

        ENDIF
        NG=NG+1
      ENDDO row
      ENDDO col
!=======================================================================
! scale the forces on the unit cell 2 Pi (and e^2/Omea/(2 pi)^2
! and add the  contribution to the total force from the
! 1/(unit cell volume)  dependence of the electron-ion energy
!=======================================================================
      EISIF(2,1)=EISIF(1,2)
      EISIF(3,2)=EISIF(2,3)
      EISIF(1,3)=EISIF(3,1)

      DSIF(2,1)=DSIF(1,2)
      DSIF(3,2)=DSIF(2,3)
      DSIF(1,3)=DSIF(3,1)

      SCALE= EDEPS/LATT_CUR%OMEGA/TPI**2
      DENC = DENC*SCALE/2

      DO K=1,3
        DO M=1,3
          EISIF(M,K)=EISIF(M,K)*TPI
          DSIF(M,K) =DSIF(M,K) *SCALE
        ENDDO
        PSCSIF(K,K)=PSCSIF(K,K)+PSCENC ! q-> 0 term
        DSIF  (K,K)=DSIF  (K,K)+DENC   ! add 1/ omega term
        EISIF (K,K)=EISIF (K,K)+PSCV   ! add 1/ omega term
      ENDDO

      
      

      DEALLOCATE(CWORK1,CWORK2)
      RETURN
      END SUBROUTINE

!************************ SUBROUTINE STREHAR ***************************
!
! this calculates the correction to the stress
! +   if the harris-functional is used                   LPAR=.FALSE.
!        (CHGGA = gradient of Harris functional)
! +   stress due to partial core corrections
!        (CHGGA = CVXC)                                  LPAR=.TRUE.
!
! LUNIF determines if the 1/volume  dependence of the charge is taken
!       into account
! for LPAR=.TRUE. LUNIF=.F. because this part is calculated in POTEX
!
!***********************************************************************

      SUBROUTINE STREHAR(GRIDC,P,T_INFO,LATT_CUR,LPAR,CHGGA,CSTRF, HARSIF)
      USE prec


      USE charge
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      COMPLEX(q) CHGGA(GRIDC%MPLWV)
      REAL(q)    HARSIF(3,3)
      LOGICAL LUNIF,LPAR
! work arrays
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK2(:)

      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
!=======================================================================
! first calculate the chargedensity and its derivative
!=======================================================================
      LUNIF=.NOT. LPAR
      CALL RHOATO(.FALSE.,LPAR,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CWORK1,CWORK2)
!=======================================================================
! now sum up all components
! as in STRELO there is (1._q,0._q) 1/ volume term plus other terms
!=======================================================================

      HARSIF=0
      PSCV  =0       ! required to calc 1/ volume term

      NG=1
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        FACTM=1
        

        GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
        GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
        GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

        GSQU=GX**2+GY**2+GZ**2
        GMOD=SQRT(GSQU)
        PSCV=PSCV+ REAL( CONJG(CHGGA(NG))*CWORK1(NG) ,KIND=q)
!=======================================================================
! avoid G=0
!=======================================================================
        IF(GMOD>1E-7_q) THEN
         PSCD=  REAL( CONJG(CHGGA(NG))*CWORK2(NG) ,KIND=q)
         HARSIF(1,1)=HARSIF(1,1)+PSCD*GX*GX/GMOD
         HARSIF(2,2)=HARSIF(2,2)+PSCD*GY*GY/GMOD
         HARSIF(3,3)=HARSIF(3,3)+PSCD*GZ*GZ/GMOD
         HARSIF(1,2)=HARSIF(1,2)+PSCD*GX*GY/GMOD
         HARSIF(2,3)=HARSIF(2,3)+PSCD*GY*GZ/GMOD
         HARSIF(3,1)=HARSIF(3,1)+PSCD*GZ*GX/GMOD
        ENDIF
        NG=NG+1
      ENDDO row
      ENDDO col

!=======================================================================
! scale the forces on the unit cell by 2 Pi
! and add the  contribution to the total force from the
! 1/(unit cell volume)  dependence of the energy if LUNIF is TRUE
!=======================================================================
      HARSIF(2,1)=HARSIF(1,2)
      HARSIF(3,2)=HARSIF(2,3)
      HARSIF(1,3)=HARSIF(3,1)

      DO K=1,3
        DO M=1,3
          HARSIF(M,K)=HARSIF(M,K)*TPI
        ENDDO
      IF (LUNIF) HARSIF (K,K)=HARSIF (K,K)+PSCV
      ENDDO

      

      DEALLOCATE(CWORK1,CWORK2)
      RETURN
      END SUBROUTINE

!************************ SUBROUTINE STRKIN ****************************
!
! this subroutine calculates the stress resulting
! from the change in the kinetic
! energy of the plane wave basis states as the size of the cell changes
! formula (10.50) in thesis gK
!***********************************************************************

      SUBROUTINE STRKIN(W,WDES, B,SIKEF)
      USE prec


      USE constant
      USE wave
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W

      DIMENSION SIKEF(3,3)
      DIMENSION B(3,3)
      DIMENSION QSP(3)
! work arrays

      SIKEF=0
      spin: DO ISP=1,WDES%ISPIN
!=======================================================================
! loop over all k-points and bands
!=======================================================================
      band_k: DO NCOUNT=0,WDES%NBANDS*WDES%NKPTS-1

      NK=NCOUNT/WDES%NBANDS+1
      N=MOD(NCOUNT,WDES%NBANDS)+1
      NPL=WDES%NGVECTOR(NK)

      SIKE1=0
      SIKE2=0
      SIKE3=0
      SIKE12=0
      SIKE23=0
      SIKE31=0
!-MM- changes to accommodate spin spirals      
      QSP=WDES%QSPIRAL/2
!-MM- end of addition
      spinor: DO ISPINOR=0,WDES%NRSPINORS-1
      DO M=1,NPL
        MM=M+NPL*ISPINOR
        G1=WDES%IGX(M,NK)+WDES%VKPT(1,NK)-QSP(1)
        G2=WDES%IGY(M,NK)+WDES%VKPT(2,NK)-QSP(2)
        G3=WDES%IGZ(M,NK)+WDES%VKPT(3,NK)-QSP(3)

        GX= (G1*B(1,1)+G2*B(1,2)+G3*B(1,3)) *TPI
        GY= (G1*B(2,1)+G2*B(2,2)+G3*B(2,3)) *TPI
        GZ= (G1*B(3,1)+G2*B(3,2)+G3*B(3,3)) *TPI

        CPT   =W%CPTWFP(MM,N,NK,ISP)
        WFMAG =CPT*CONJG(CPT)

        SIKE1 =SIKE1 + GX*GX*HSQDTM *WFMAG
        SIKE2 =SIKE2 + GY*GY*HSQDTM *WFMAG
        SIKE3 =SIKE3 + GZ*GZ*HSQDTM *WFMAG
        SIKE12=SIKE12+ GX*GY*HSQDTM *WFMAG
        SIKE23=SIKE23+ GY*GZ*HSQDTM *WFMAG
        SIKE31=SIKE31+ GZ*GX*HSQDTM *WFMAG

      ENDDO
!-MM- changes to accommodate spin spirals      
      QSP=-QSP
!-MM- end of addition
      ENDDO spinor
!=======================================================================
! sum the contributions to the force on the unit cell from each band
! multiplying by the
! k point weight and fermi-weight, by WDES%RSPIN for electron spins and
! and additional factor 2
!=======================================================================
      SCALE=  WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*2._q*WDES%RSPIN

      SIKEF(1,1)=SIKEF(1,1)+SIKE1 *SCALE
      SIKEF(2,2)=SIKEF(2,2)+SIKE2 *SCALE
      SIKEF(3,3)=SIKEF(3,3)+SIKE3 *SCALE
      SIKEF(1,2)=SIKEF(1,2)+SIKE12*SCALE
      SIKEF(2,3)=SIKEF(2,3)+SIKE23*SCALE
      SIKEF(3,1)=SIKEF(3,1)+SIKE31*SCALE
      SIKEF(2,1)=SIKEF(2,1)+SIKE12*SCALE
      SIKEF(3,2)=SIKEF(3,2)+SIKE23*SCALE
      SIKEF(1,3)=SIKEF(1,3)+SIKE31*SCALE

      ENDDO band_k
      ENDDO spin
! reduction of SIKEF
      

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE CHGGRA ****************************
!
! this subroutine calculates the gradient vector associated with
! Harris-functional
! i.e.
!  CHTOTL contains the input-chargedenisty
!  CHTOT  the output-chargedenisty
!  CHGGA  is the final result i.e
!***********************************************************************

      SUBROUTINE CHGGRA(GRIDC,LATT_CUR,EXCTAB, &
     &     CHGGA,CHTOT,CHTOTL,DENCOR)
      USE prec

      USE mpimy
      USE mgrid
      USE lattice
      USE setexm
      USE xcgrad
      USE charge
      USE pot
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (exctable)    EXCTAB
      TYPE (latt)        LATT_CUR

      COMPLEX(q)  CHTOT(GRIDC%RC%NP),CHTOTL(GRIDC%RC%NP)
      COMPLEX(q)  CHGGA(GRIDC%MPLWV)
      COMPLEX(q)       DENCOR(GRIDC%RL%NP)
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWORK1(:),CWORK2(:)
      REAL(q), ALLOCATABLE  :: DWORK2(:)
      REAL(q) TMPSIF(3,3)

      ALLOCATE (DWORK2(GRIDC%RL%NP),CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
!-----------------------------------------------------------------------
!  we need the  input-charge in real space so transform CHTOTL in
!  real space -> use here CHGGA as a work array ...
!-----------------------------------------------------------------------
      DO N=1,GRIDC%RC%NP
        CHGGA(N)=CHTOTL(N)
      ENDDO
      CALL FFT3RC(CHGGA,GRIDC,1)
!-----------------------------------------------------------------------
!  calculate the gradient of the exchange correlation potential
!-----------------------------------------------------------------------
      CALL FEXCP(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
             CHGGA,DENCOR,CWORK1,DWORK2,CVZERO,EXC,XCENC,TMPSIF,.FALSE.)
!-----------------------------------------------------------------------
!  calculate difference between output and input charge-density in
!  real space -> store result in CHGGA
!-----------------------------------------------------------------------
      DO N=1,GRIDC%RC%NP
        CHGGA(N)=CHTOT(N)-CHTOTL(N)
      ENDDO
      CALL FFT3RC(CHGGA,GRIDC,1)
!-----------------------------------------------------------------------
! multiply by gradient of xc-potential
! in real space and FFT this to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      CALL RLR_MUL(CHGGA,DWORK2,RINPL/LATT_CUR%OMEGA,CHGGA,GRIDC)
      CALL FFT3RC(CHGGA,GRIDC,-1)
      CALL SETUNB(CHGGA,GRIDC)
!-----------------------------------------------------------------------
! add hartree-term
!-----------------------------------------------------------------------
      DO N=1,GRIDC%RC%NP
        CWORK1(N)=CHTOT(N)-CHTOTL(N)
      ENDDO

      CALL POTHAR(GRIDC,LATT_CUR, CWORK1,CWORK2,DENC)

      DO N=1,GRIDC%RC%NP
       CHGGA(N)=CHGGA(N)+CWORK2(N)
      ENDDO
      DEALLOCATE (DWORK2,CWORK1,CWORK2)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE FORLOC ****************************
!
! this subroutine calculates the hellmann-feynman forces exerted on the
! ions by the electrons which equals the sum over reciprocal lattice
! vectors of  c.c. of the charge density at wavevector G
!  *IG*EXP(+IG.R)* pseudopotential at wavevector G
!***********************************************************************

      SUBROUTINE FORLOC(GRIDC,P,T_INFO,LATT_CUR, &
              CHTOT,EIFOR)
      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%RC%NP)
      REAL(q)    EIFOR(3,T_INFO%NIONS)
! work arrays
      REAL(q), ALLOCATABLE :: WORK(:)

      ALLOCATE(WORK(GRIDC%RC%NP))

      NIS=1
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
!=======================================================================
       NIADD=T_INFO%NITYP(NT)
!=======================================================================
! interpolate the pseudopotential on the grid of reciprocal
! lattice-vectors
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS

      ZZ=  -4*PI*P(NT)%ZVALF*FELECT

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)

        FACTM=1
        
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI
        IF (G/=0 .AND. G <PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal lattice vector to a position
! in the pseudopotential arrays and interpolate the pseudopotential and
! its derivative
!=======================================================================
        I  =INT(G*ARGSC)+1
        REM=G-P(NT)%PSP(I,1)
        VPST =(P(NT)%PSP(I,2)+REM*(P(NT)%PSP(I,3)+ &
     &                     REM*(P(NT)%PSP(I,4)  +REM*P(NT)%PSP(I,5))))
        WORK (N)=( VPST+ ZZ / G**2) /LATT_CUR%OMEGA

        ELSE
          WORK(N)=0
        ENDIF

      ENDDO

      ion: DO NI=NIS,NIADD+NIS-1
!=======================================================================
! initialise the force on the ion to (0._q,0._q)
!=======================================================================
         FOR1=0
         FOR2=0
         FOR3=0
!=======================================================================
! CGXDX,Y,Z = I* the changes in the phase factor g.r on moving (1._q,0._q)
! reciprocal lattice vector in the x,y,z directions, respectively
!=======================================================================
!=======================================================================
! calculate the total force on the ions by summing over reciprocal
! lattice vectors
! first calculate phase factor:
! there are two version for calculating the phase factor
! on vector machines you might try the first version
! (see stufak.F)
!=======================================================================
         CX =EXP(-CITPI*T_INFO%POSION(1,NI))
         G1 =T_INFO%POSION(1,NI)*(-(GRIDC%NGX/2-1))

         DO NC=1,GRIDC%RC%NCOL
           NGP=(NC-1)*GRIDC%RC%NROW+1

           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)
           G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
           G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
           CE=EXP(-CITPI*(G3+G2+G1))

           DO N1P=0,GRIDC%RC%NROW-1
           N1=MOD(N1P+(-(GRIDC%NGX/2-1))+GRIDC%NGX,GRIDC%NGX)
           NG=NGP+N1
           N1=N1+1

           FACTM=1
           
           CEXPF=CE
           CE=CE*CX
!=======================================================================
! add the contribution to the force from the present reciprocal lattice
! vector  and multiply by i (ie take imaginary part)
!=======================================================================
           FOR=WORK(NG)* AIMAG(CONJG(CHTOT(NG))*CEXPF)
           FOR1=FOR1-GRIDC%LPCTX_(N1)*FOR
           FOR2=FOR2-GRIDC%LPCTY_(N2)*FOR
           FOR3=FOR3-GRIDC%LPCTZ_(N3)*FOR
         ENDDO
         ENDDO
!=======================================================================
! multiply forces by 2*Pi
!=======================================================================
         EIFOR(1,NI)=FOR1*TPI
         EIFOR(2,NI)=FOR2*TPI
         EIFOR(3,NI)=FOR3*TPI

      ENDDO ion
      NIS=NIS+NIADD
      ENDDO typ

!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
      

      CALL  DIRKAR(T_INFO%NIONS,EIFOR,LATT_CUR%B)
      DEALLOCATE(WORK)


      RETURN
      END SUBROUTINE


!************************ SUBROUTINE FORHAR ****************************
!
! this subroutine calculates the correction to the
! hellman feynman forces
! +   if the harris-functional is used
!        CHGGA = gradient of Harris functional,
!           PSPRHO = atom. charge                LPAR=.FALSE.
! +   forces due to partial core corrections
!        CHGGA = CVXC(xc-potential),
!           PSPRHO=partial core                  LPAR=.TRUE.
!***********************************************************************


      SUBROUTINE FORHAR(GRIDC,P,T_INFO,LATT_CUR, &
              CHGGA,HARFOR,LPAR)
      USE prec

      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR
      COMPLEX(q)   CHGGA(GRIDC%MPLWV)
      REAL(q)      HARFOR(3,T_INFO%NIONS)
      LOGICAL   LPAR
! work arrays
      REAL(q), POINTER :: PRHO(:)
      REAL(q), ALLOCATABLE :: WORK(:)

      ALLOCATE(WORK(GRIDC%RC%NP))

      HARFOR(1:3,1:T_INFO%NIONS)=0
!=======================================================================
! loop over all types of atoms
!=======================================================================
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
      NIADD=T_INFO%NITYP(NT)
       IF (LPAR) THEN
         PRHO=>P(NT)%PSPCOR
       ELSE
         PRHO=>P(NT)%PSPRHO
       ENDIF
       IF (.NOT.ASSOCIATED(PRHO)) GOTO 200
!=======================================================================
! iterpolate the pseudopotential on the grid of reciprocal
! lattice-vectors
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)

        FACTM=1
        
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/.6_q
          WORK(N)=T0+REM*(T1+REM*(T2+REM*T3))
        ELSE
          WORK(N)=0
        ENDIF

      ENDDO

      ion: DO NI=NIS,NIADD+NIS-1
!=======================================================================
! initialise the force on the ion to (0._q,0._q)
!=======================================================================
         FOR1=0._q
         FOR2=0._q
         FOR3=0._q
!=======================================================================
! calculate the total force on the ions by summing over reciprocal
! lattice vectors
! first calculate phase factor:
! there are two version for calculating the phase factor
! on vector machines you might try the first version
! but usually the second version is much faster (less calls to CEXP)
!=======================================================================

         CX =EXP(-CITPI*T_INFO%POSION(1,NI))
         G1 =T_INFO%POSION(1,NI)*(-(GRIDC%NGX/2-1))

         DO NC=1,GRIDC%RC%NCOL
           NGP=(NC-1)*GRIDC%RC%NROW+1

           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)
           G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
           G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
           CE=EXP(-CITPI*(G3+G2+G1))

           DO N1P=0,GRIDC%RC%NROW-1
           N1=MOD(N1P+(-(GRIDC%NGX/2-1))+GRIDC%NGX,GRIDC%NGX)
           NG=NGP+N1
           N1=N1+1
           FACTM=1
           
           CEXPF=CE
           CE=CE*CX

!=======================================================================
! add the contribution to the force from the present reciprocal lattice
! vector  and multiply by i (ie take imaginary part)
!=======================================================================
           FOR=WORK(NG)* AIMAG(CONJG(CHGGA(NG))*CEXPF)
           FOR1=FOR1-GRIDC%LPCTX_(N1)*FOR
           FOR2=FOR2-GRIDC%LPCTY_(N2)*FOR
           FOR3=FOR3-GRIDC%LPCTZ_(N3)*FOR
         ENDDO

         ENDDO

!=======================================================================
! multiply forces by 2*Pi
!=======================================================================
         HARFOR(1,NI)=FOR1*TPI
         HARFOR(2,NI)=FOR2*TPI
         HARFOR(3,NI)=FOR3*TPI

      ENDDO ion
  200 NIS=NIS+NIADD
      ENDDO typ
!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
      
      CALL  DIRKAR(T_INFO%NIONS,HARFOR,LATT_CUR%B)
      DEALLOCATE(WORK)

      RETURN
      END SUBROUTINE
