#include "symbol.inc"

#define pro_details

      MODULE msphpro
      USE prec
      CHARACTER (LEN=3), ALLOCATABLE, SAVE :: LMTABLE(:, :)

      CONTAINS
!************************ SUBROUTINE SPHPRO ****************************
! RCS:  $Id: sphpro.F,v 1.7 2003/06/27 13:22:23 kresse Exp kresse $
!
! SPHPRO calculates the projection of the wavefunctions onto spherical
! waves and from that the local charge on each ion and the 
! partial density of states
!
!***********************************************************************

      SUBROUTINE SPHPRO( &
          GRID,LATT_CUR,LATT_INI, P,T_INFO,W,WDES, IUP,IU6, &
          LOVERL,LMDIM,CQIJ, LPAR, LDIMP,LMDIMP,LTRUNC, LORBIT,PAR)
      USE prec
      USE main_mpi
      USE constant
      USE wave
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE nonl
      USE relativistic

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NIONS)
      TYPE (potcar)      BET(1)
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES,WDES_1K
      TYPE (nonl_struct) NONL_S

      LOGICAL LOVERL
      INTEGER LORBIT
! OLebacq:begin      
      COMPLEX(q), PARAMETER :: ci=(0._q,1._q)
      INTEGER LMINDp,LMSp
      REAL(q), ALLOCATABLE ::  PAR_lmom(:,:,:,:,:)
      REAL(q) SUMAUG_lmom(LMDIMP,T_INFO%NIONP,WDES%NCDIJ),PARAUG_lmom(LDIMP,T_INFO%NIONP,WDES%NCDIJ)  
      REAL(q) ION_SUM_lmom(LDIMP,T_INFO%NIONP,WDES%NCDIJ) 
      REAL(q) SUMION_lmom(LDIMP),SUMTOT_lmom(LDIMP)
      COMPLEX(q) CSUM_ABS_lmom(LMDIMP,WDES%NCDIJ)
      COMPLEX(q), ALLOCATABLE :: L_OP(:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: DUMMY(:,:,:,:)
   
! OLebacq:end      

      OVERLAP   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q) PAR(WDES%NBANDS,WDES%NKDIM,LPAR ,T_INFO%NIONP,WDES%NCDIJ)
!  allocate dynamically
      REAL(q) WORKR(WDES%NRPLWV), WORKI(WDES%NRPLWV)
      COMPLEX(q) CSUM_ABS(LMDIMP,WDES%NCDIJ),CSUM_PHASE(LMDIMP,WDES%NCDIJ)
      COMPLEX(q) CSUM(LMDIMP*LTRUNC,WDES%NCDIJ)
      REAL(q)  SUMAUG(LMDIMP,T_INFO%NIONP,WDES%NCDIJ), &
               PARAUG(LPAR,T_INFO%NIONP,WDES%NCDIJ)
      REAL(q) SUMION(LPAR),SUMTOT(LPAR)

      REAL(q) ION_SUM(LPAR,T_INFO%NIONP,WDES%NCDIJ)
      REAL(q) ION_SUM2(LDIMP,T_INFO%NIONP,WDES%NCDIJ)

      COMPLEX(q),ALLOCATABLE :: PHAS(:,:,:,:,:)
      REAL(q) QZERO(LTRUNC)
      CHARACTER*(8) :: STR
      CHARACTER(LEN=3) :: LCHAR(4)=(/'  s','  p','  d','  f'/)
      CHARACTER(LEN=3) :: LMCHAR(16)=(/'  s',' py',' pz',' px','dxy', &
         'dyz','dz2','dxz','dx2','f-3','f-2','f-1',' f0',' f1',' f2',' f3'/)

      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME= WDES%COMM%NODE_ME
      IONODE = WDES%COMM%IONODE
#endif
!=======================================================================
! some initialization
!=======================================================================

! Get me the L operator
      ALLOCATE(L_OP(2*LDIMP+1,2*LDIMP+1,3,0:LDIMP),DUMMY(2*LDIMP+1,2*LDIMP+1,4,0:LDIMP))

      IF (WDES%LNONCOLLINEAR) THEN
         ALLOCATE(PAR_lmom(WDES%NBANDS,WDES%NKDIM,LDIMP,T_INFO%NIONP,WDES%NCDIJ))
      ENDIF
           
      L_OP=(0._q,0._q)      
      
      DO L=1,LDIMP
         CALL SETUP_LS(L,0._q,0._q,L_OP(1:2*L+1,1:2*L+1,1:3,L),DUMMY(1:2*L+1,1:2*L+1,1:4,L))
      ENDDO

      IF (LORBIT==2) THEN
         ALLOCATE(PHAS(LMDIMP,T_INFO%NIONP,WDES%NKDIM,WDES%NBANDS,WDES%ISPIN))
         PHAS=0
      ENDIF

      RSPIN=WDES%RSPIN

!     compute volume per type
      io_begin
      DO NT=1,T_INFO%NTYPP
         IF ( NT<=T_INFO%NTYP ) THEN
            WRITE(IU6,2220) NT,100*T_INFO%NITYP(NT)*2*TPI*T_INFO%RWIGS(NT)**3/3/LATT_CUR%OMEGA
         ELSE
            WRITE(IU6,2221) NT,100*T_INFO%NITYP(NT)*2*TPI*T_INFO%RWIGS(NT)**3/3/LATT_CUR%OMEGA
         ENDIF

 2220  FORMAT(/'volume of typ          ',I3,':  ',F6.1,' %')
 2221  FORMAT(/'volume of empty sphere ',I3,':  ',F6.1,' %')
      ENDDO
!   write header fo file PROOUT

      IF (LORBIT==5) THEN
         DO ISP=1,WDES%NCDIJ
           WRITE(STR,'(A,I1)') "PROOUT.",ISP
           OPEN(UNIT=IUP+ISP-1,FILE=DIR_APP(1:DIR_LEN)//STR,STATUS='UNKNOWN')
         ENDDO
         DO ISP=1,WDES%NCDIJ
           WRITE(IUP+ISP-1,*) 'PROOUT'
           WRITE(IUP+ISP-1,3200)    WDES%NKPTS,WDES%NB_TOT,T_INFO%NIONP
           WRITE(IUP+ISP-1,'(9I4)') T_INFO%NTYPP,T_INFO%NTYP,(T_INFO%NITYP(I),I=1,T_INFO%NTYPP)
           WRITE(IUP+ISP-1,'(9F7.3)') ((W%FERTOT(NB,NK,ISP),NK=1,WDES%NKPTS),NB=1,WDES%NB_TOT)
         ENDDO
      ELSE
         OPEN(UNIT=IUP,FILE=DIR_APP(1:DIR_LEN)//'PROCAR',STATUS='UNKNOWN')
         IF (LORBIT==1) THEN
            WRITE(IUP,'(A)')'PROCAR lm decomposed'
         ELSE  IF (LORBIT==2) THEN
            WRITE(IUP,'(A)')'PROCAR lm decomposed + phase factor'
         ELSE
            WRITE(IUP,'(A)')'PROCAR new format'
         ENDIF
      ENDIF
      io_end


!   allocate descriptor for 1 kpoints
      CALL CREATE_SINGLE_KPOINT_WDES(WDES,WDES_1K,1)
!   allocate BET
      ALLOCATE(BET(1)%PSPNL(0:NPSNL,LDIMP*LTRUNC), &
               BET(1)%LPS(LDIMP*LTRUNC))
      BET(1)%LMDIM   =LMDIMP*LTRUNC
      BET(1)%LMAX    =LDIMP*LTRUNC
      BET(1)%LMMAX   =LMDIMP*LTRUNC

      CALL NONL_ALLOC_SPHPRO(NONL_S,BET(1),WDES_1K)
!=======================================================================
      NIS=1

      typ: DO NT=1,T_INFO%NTYPP
!-----------------------------------------------------------------------
!    set up table with BETAs:
!
!    BET(1)%PSPNL=  4 Pi \int_0^R_{cut}
!                j_l(qr) j_l(q_n r) r^2 dr / Sqrt(A(q,qp))
!           q_n is chosen so that j_l(q_n R_{cut}) = 0
!     we use the relation
!     A(q,qp) = \int_0^R j_l(qr) j_l(qp r) r^2 dr
!             = R^2/(q^2-qp^2) [ j(qR) j(qpR)' - j(qp R) j(qR)')
!-----------------------------------------------------------------------
      GMAX =SQRT(WDES%ENMAX/HSQDTM)*1.2_q

      BET(1)%PSMAXN = GMAX

      PSTEP=GMAX/(NPSNL-1)

      INDEX=0
    setbet: DO LL= 0,LDIMP-1
      CALL BEZERO(QZERO,LL,LTRUNC)
      DO I=1,LTRUNC
        INDEX=INDEX+1

        Q1=QZERO(I)/T_INFO%RWIGS(NT)
        QR=QZERO(I)
        CALL SBESSE2(QR, BQ, BQP, LL)
        A= 1/Q1* T_INFO%RWIGS(NT)**2/2*(BQ*BQP+QR*BQP**2+ &
                            QR*BQ**2- (LL+1)*LL/QR*BQ**2)
        SQRTIA=1/SQRT(A)
          DO N = 0,NPSNL-1
            QQ=PSTEP*N
            IF (QQ==0) QQ=1E-5_q
            CALL SBESSE2(T_INFO%RWIGS(NT)*QQ, BQQ , BQQP, LL)
            IF (ABS(QQ-Q1)<1E-5_q) THEN
            A= 1/Q1* T_INFO%RWIGS(NT)**2/2*(BQ*BQP+QR*BQP**2+ &
                            QR*BQ**2- (LL+1)*LL/QR*BQ**2)
            ELSE
            A= T_INFO%RWIGS(NT)**2/(Q1**2-QQ**2)*(BQ*BQQP*QQ-BQQ*BQP*Q1)
            ENDIF
            BET(1)%PSPNL(N+1,INDEX)=TPI*2*A*SQRTIA
        ENDDO
      IF (MOD(LL,2)==0) THEN
        BET(1)%PSPNL(0,INDEX)=BET(1)%PSPNL(2,INDEX)
      ELSE
        BET(1)%PSPNL(0,INDEX)=BET(1)%PSPNL(2,INDEX)
      ENDIF
      BET(1)%LPS(INDEX)=LL
      ENDDO
      ENDDO setbet
!=======================================================================
      kpoint: DO NK=1,WDES%NKPTS
!=======================================================================
      CALL CREATE_SINGLE_KPOINT_WDES(WDES,WDES_1K,NK)

      IZERO =1
      CALL SPHER(GRID,NONL_S,BET,WDES_1K,LATT_CUR,  IZERO,LATT_INI%B)
!=======================================================================
      ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
!=======================================================================
      NONL_S%POSION=>T_INFO%POSION(:,NI:NI)
      CALL PHASE(WDES_1K,NONL_S,0)  ! reset phase factor
      CALL PHASE(WDES_1K,NONL_S,1)  ! and force calculation
! phase factor e(i k R)
      GXDX=T_INFO%POSION(1,NI)
      GYDY=T_INFO%POSION(2,NI)
      GZDZ=T_INFO%POSION(3,NI)
      CGDR=EXP(CITPI*(WDES%VKPT(1,NK)*GXDX+WDES%VKPT(2,NK)*GYDY+WDES%VKPT(3,NK)*GZDZ))

      band: DO NB=1,WDES%NBANDS
!=======================================================================
! multiply with phasefactor and divide into real and imaginary part
!=======================================================================
      NPL = WDES%NGVECTOR(NK)

      CSUM_PHASE=0
      CSUM_ABS  =0
      CSUM_ABS_lmom =0

!-MM- spin spiral stuff
      ISPIRAL = 1
!-MM- end of addition
      spin: DO ISP=1,WDES%ISPIN
      DO ISPINOR=0,WDES%NRSPINORS-1

      DO K=1,NPL
         KK=K+ISPINOR*NPL
         CTMP=    NONL_S%CREXP(K,1) * W%CPTWFP(KK,NB,NK,ISP)
         WORKR(K) = REAL( CTMP ,KIND=q)
         WORKI(K) = AIMAG(CTMP)
      ENDDO
!=======================================================================
! loop over indices L,M,N and calculate 
! CSUM(lmn,alp) = < phi(alpha) | beta lmn >
!=======================================================================

      LMS=0
      DO LL=0,LDIMP-1
         DO I=1,LTRUNC
            DO M=0,2*LL
               LMS=LMS+1
               LM=LL*LL+M+1
               SUMR=0
               SUMI=0
!DIR$ IVDEP
!OCL NOVREC
              ! loop over G-vectors
               DO K=1,NPL
!-MM- changes to accomodate spin spirals
! original statements
!                 SUMR = SUMR + WORKR(K) * NONL_S%QPROJ(K,LMS,1,1)
!                 SUMI = SUMI + WORKI(K) * NONL_S%QPROJ(K,LMS,1,1)
                  SUMR = SUMR + WORKR(K) * NONL_S%QPROJ(K,LMS,1,1,ISPIRAL)
                  SUMI = SUMI + WORKI(K) * NONL_S%QPROJ(K,LMS,1,1,ISPIRAL)
!-MM- end of alterations
               ENDDO
               CTMP=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *NONL_S%CQFAK(LMS,1)*CGDR)
               CSUM(LMS,ISP+ISPINOR)     =CTMP
               CSUM_PHASE(LM,ISP+ISPINOR)=CSUM_PHASE(LM,ISP+ISPINOR)+CTMP*ABS(CTMP)
            ENDDO
         ENDDO
      ENDDO
!-MM- spin spiral stuff
      IF (NONL_S%LSPIRAL) ISPIRAL=2
!-MM- end of addition
      ENDDO
      ENDDO spin
      CALLMPI( M_sum_z( WDES%COMM_INB, CSUM, LMDIMP*LTRUNC*(WDES%ISPIN+WDES%NRSPINORS-1)))
      CALLMPI( M_sum_z( WDES%COMM_INB, CSUM_PHASE, LMDIMP*(WDES%ISPIN+WDES%NRSPINORS-1)))

      IF (LORBIT==2) THEN
         IF(WDES%LNONCOLLINEAR)  THEN
            ! the phase factor is only qualitative, we just sum over up and down
            ! for the non collinear case
            PHAS(:,NI,NK,NB,1)= & 
             CSUM_PHASE(:,1)*SQRT(ABS(CSUM_ABS(:,1)))/ABS(CSUM_PHASE(:,1))+ &
             CSUM_PHASE(:,2)*SQRT(ABS(CSUM_ABS(:,2)))/ABS(CSUM_PHASE(:,2))
         ELSE
            PHAS(:,NI,NK,NB,:)=CSUM_PHASE(:,1:WDES%ISPIN)* &
               SQRT(ABS(CSUM_ABS(:,1:WDES%ISPIN)))/ABS(CSUM_PHASE(:,1:WDES%ISPIN))
         ENDIF
      ENDIF

      IF (LMS /= LMDIMP*LTRUNC) THEN
         WRITE(*,*) 'internal error 1:',LMS,LMDIMP*LTRUNC
         STOP
      ENDIF

      !
      ! now calculate rho(lm,alp,bet) = 
      ! sum_n < phi(alp) | beta lmn > <beta lmn | phi(bet) >
      !

      DO ISP=1,WDES%ISPIN
      DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1
      
      II=ISP+ISPINOR_+2*ISPINOR

      LMS=0
      DO LL=0,LDIMP-1
         DO I=1,LTRUNC
            DO M =1,2*LL+1
               LMS=LMS+1
               LM=LL*LL+M
               CSUM_ABS(LM,II)  = CSUM_ABS(LM,II)+ &
                    CSUM(LMS,ISP+ISPINOR)*CONJG(CSUM(LMS,ISP+ISPINOR_))               
            ENDDO
         ENDDO
      ENDDO

! OLebacq : calculate the orbital moment outside the augmentation region
      IF (WDES%LNONCOLLINEAR) THEN
      LMS=0
      DO LL=0,LDIMP-1
         DO I=1,LTRUNC
            LMS_BASE=LMS
            DO M =1,2*LL+1
               LMS=LMS+1
               LM=LL*LL+M
               LMS_=LMS_BASE
               DO M_=1,2*LL+1
                  LMS_=LMS_+1
                  CPRODCSUM_x=0._q;CPRODCSUM_y=0._q;CPRODCSUM_z=0._q
                  IF (ISPINOR==ISPINOR_) THEN
                     CPRODCSUM_x=CSUM(LMS,ISP+ISPINOR)*CONJG(CSUM(LMS_,ISP+ISPINOR_))*L_OP(M_,M,1,LL)
                     CPRODCSUM_y=CSUM(LMS,ISP+ISPINOR)*CONJG(CSUM(LMS_,ISP+ISPINOR_))*L_OP(M_,M,2,LL)
                     CPRODCSUM_z=CSUM(LMS,ISP+ISPINOR)*CONJG(CSUM(LMS_,ISP+ISPINOR_))*L_OP(M_,M,3,LL) 
                  ENDIF
                  CSUM_ABS_lmom(LM,1)  = CSUM_ABS_lmom(LM,1)+ CPRODCSUM_x
                  CSUM_ABS_lmom(LM,2)  = CSUM_ABS_lmom(LM,2)+ CPRODCSUM_y
                  CSUM_ABS_lmom(LM,3)  = CSUM_ABS_lmom(LM,3)+ CPRODCSUM_z
               ENDDO         
            ENDDO
         ENDDO
      ENDDO   
      ENDIF
! OLebacq : end orbital moment part outside the spheres         

      ENDDO
      ENDDO
      ENDDO
      

      IF (WDES%LNONCOLLINEAR) THEN
      CALL C_FLIP(CSUM_ABS,LMDIMP,LMDIMP,WDES%NCDIJ,.FALSE.)
      ENDIF

      DO ISP=1,WDES%NCDIJ
         DO LL=0,LDIMP-1
            IF (LORBIT==1.OR.LORBIT==2) THEN
               DO M=0,2*LL
                  LM=LL*LL+M+1
                  PAR(NB,NK,LM,NI,ISP)=CSUM_ABS(LM,ISP)
               ENDDO
            ELSE
               SUML=0
               SUML2=0
               DO M=0,2*LL
                  LM=LL*LL+M+1
                  SUML=SUML  +CSUM_ABS(LM,ISP)
                  SUML2=SUML2+CSUM_ABS_lmom(LM,ISP)
               ENDDO
               PAR(NB,NK,LL+1,NI,ISP)=SUML
               IF (WDES%LNONCOLLINEAR) PAR_lmom(NB,NK,LL+1,NI,ISP)=SUML2
            ENDIF
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
      io_begin
      IF (LORBIT==5) THEN
         CSUM_PHASE=CSUM_PHASE* &
         SQRT(ABS(CSUM_ABS))/ABS(CSUM_PHASE)
         DO ISP=1,WDES%NCDIJ
           WRITE(IUP+ISP-1,'(9F12.6)')  CSUM_PHASE(:,ISP)
         ENDDO
      ENDIF
      io_end

      ENDDO band
      ENDDO ion
      ENDDO kpoint

      NIS=NIS+T_INFO%NITYP(NT)
      ENDDO typ

      CALL NONL_DEALLOC_SPHPRO(NONL_S)
      DEALLOCATE(BET(1)%PSPNL,BET(1)%LPS)
!-----------------------------------------------------------------------
      io_begin
      IF (LORBIT==5) THEN
        DO ISP=1,WDES%NCDIJ
         WRITE(IUP+ISP-1,*) 'augmentation part'
        ENDDO
      ENDIF
      io_end
!=======================================================================
! calculate contribution from augmentation-part
!=======================================================================
      overl: IF (LOVERL) THEN
      SUMALL_AUG=0

    kpoint_aug: DO NK=1,WDES%NKPTS
    band_aug:   DO N=1 ,WDES%NBANDS

      NIS=1
      NPRO=0
      SUMAUG=0
      PARAUG=0
      SUMAUG_lmom=0
      PARAUG_lmom=0

    typ_aug:    DO NT= 1,T_INFO%NTYP
    ion_aug:    DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
      NIP=NI_LOCAL(NI, WDES%COMM_INB)     !  local storage index
      IF (NIP==0) CYCLE ion_aug
!-----------------------------------------------------------------------
! find blocks with same quantum number L
! assuming that the block is continously arranged in the arrays
!-----------------------------------------------------------------------
      LOW=1
      LMBASE=1

    block: DO

      LL=P(NT)%LPS(LOW)
      DO LHI=LOW,P(NT)%LMAX
        IF (LL/=P(NT)%LPS(LHI)) EXIT
      ENDDO

      LHI=LHI-1
!-----------------------------------------------------------------------
! only terms with equal L L' and M M' contribute
!-----------------------------------------------------------------------
      MMAX=2*LL+1
      CSUM_ABS=0
      CSUM_ABS_lmom=0

      DO ISP=1,WDES%ISPIN
      DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1
      II=ISP+ISPINOR_+2*ISPINOR

      DO L =LOW,LHI
      DO LP=LOW,LHI
      DO M =0,MMAX-1
         LMIND = NPRO+LMBASE+(L -LOW)*MMAX+M+ISPINOR *WDES%NPRO/2
         LMIND_= NPRO+LMBASE+(LP-LOW)*MMAX+M+ISPINOR_*WDES%NPRO/2
         CTMP=W%CPROJ(LMIND ,N,NK,ISP)*GCONJG(W%CPROJ(LMIND_,N,NK,ISP))*P(NT)%QION(LP,L)

         CSUM_ABS(LL*LL+M+1,II)=CSUM_ABS(LL*LL+M+1,II)+CTMP

!  OLebacq : orbital moments in the augmentation spheres
         IF(WDES%LNONCOLLINEAR)  THEN
         DO M_=0,MMAX-1
            LMINDp = NPRO+LMBASE+(LP-LOW)*MMAX+M_+ISPINOR_*WDES%NPRO/2
            CTMP_x=0._q;CTMP_y=0._q;CTMP_z=0._q                    
            IF(ISPINOR==ISPINOR_.AND.L==LP) THEN            
               CTMP_x=W%CPROJ(LMIND,N,NK,ISP)*GCONJG(W%CPROJ(LMINDp,N,NK,ISP)) * &
              &                    P(NT)%QION(LP,L) * L_OP(M_+1,M+1,1,LL)
               CTMP_y=W%CPROJ(LMIND,N,NK,ISP)*GCONJG(W%CPROJ(LMINDp,N,NK,ISP)) * &
              &                    P(NT)%QION(LP,L) * L_OP(M_+1,M+1,2,LL)
               CTMP_z=W%CPROJ(LMIND,N,NK,ISP)*GCONJG(W%CPROJ(LMINDp,N,NK,ISP)) * &
              &                    P(NT)%QION(LP,L) * L_OP(M_+1,M+1,3,LL)
            ENDIF
            CSUM_ABS_lmom(LL*LL+M+1,1)=CSUM_ABS_lmom(LL*LL+M+1,1)+CTMP_x 
            CSUM_ABS_lmom(LL*LL+M+1,2)=CSUM_ABS_lmom(LL*LL+M+1,2)+CTMP_y
            CSUM_ABS_lmom(LL*LL+M+1,3)=CSUM_ABS_lmom(LL*LL+M+1,3)+CTMP_z
         ENDDO
         ENDIF
!  OLebacq : End Orbital moment part        
      ENDDO
      ENDDO
      ENDDO

      ENDDO
      ENDDO
      ENDDO

      IF (WDES%LNONCOLLINEAR) THEN
         CALL C_FLIP(CSUM_ABS,LMDIMP,LMDIMP,WDES%NCDIJ,.FALSE.)
      ENDIF

      DO M =0,MMAX-1
         SUMAUG(LL*LL+M+1,NI,:)= CSUM_ABS(LL*LL+M+1,:)
      ENDDO
      SUMAUG_lmom(:,NI,:)= CSUM_ABS_lmom

      IF (LORBIT==1.OR.LORBIT==2) THEN
         DO ISP=1,WDES%NCDIJ
            DO M =0,MMAX-1
               PARAUG(LL*LL+M+1,NI,ISP)=PARAUG(LL*LL+M+1,NI,ISP)+REAL(CSUM_ABS(LL*LL+M+1,ISP),q)
            ENDDO
         ENDDO         
      ELSE
         DO ISP=1,WDES%NCDIJ
            PARAUG(LL+1,NI,ISP)=PARAUG(LL+1,NI,ISP)+SUM(REAL(CSUM_ABS(:,ISP),q))
!  OL : Orb mom part      
            PARAUG_lmom(LL+1,NI,ISP)=PARAUG_lmom(LL+1,NI,ISP)+SUM(REAL(CSUM_ABS_lmom(:,ISP),q))
         ENDDO
      ENDIF

      LMBASE=LMBASE+(LHI-LOW+1)*MMAX
      LOW=LHI+1
      IF (LOW >  P(NT)%LMAX) EXIT block
    ENDDO block

      NPRO=NPRO+LMBASE-1
      ENDDO ion_aug

      NIS=NIS+T_INFO%NITYP(NT)
      ENDDO typ_aug

      CALLMPI (M_sum_d( WDES%COMM, PARAUG, LPAR*T_INFO%NIONP*WDES%NCDIJ))
      PAR(N,NK,:,:,:)=PAR(N,NK,:,:,:)+PARAUG(:,:,:)

      CALLMPI (M_sum_d( WDES%COMM_INB, SUMAUG, LMDIMP*T_INFO%NIONP*WDES%NCDIJ))

!  OL : Orb mom part      
      IF(WDES%LNONCOLLINEAR)  THEN
         CALLMPI (M_sum_d( WDES%COMM,PARAUG_lmom, LDIMP*T_INFO%NIONP*WDES%NCDIJ))
         PAR_lmom(N,NK,:,:,:)=PAR_lmom(N,NK,:,:,:)+PARAUG_lmom(:,:,:)
         CALLMPI (M_sum_d( WDES%COMM_INB, SUMAUG_lmom,LMDIMP*T_INFO%NIONP*WDES%NCDIJ))
      ENDIF

      IF (LORBIT==5) THEN
         io_begin
         DO ISP=1,WDES%NCDIJ
           WRITE(IUP+ISP-1,'(9F12.6)') SUMAUG(:,:,ISP)
         ENDDO
         io_end
      ENDIF

      ENDDO band_aug
      ENDDO kpoint_aug

      ENDIF overl

!
! calculate the ionic occupancies
!
      ION_SUM=0
      ION_SUM_lmom=0._q
      DO ISP=1,WDES%NCDIJ
      ISP_=MIN(ISP,WDES%ISPIN)
      DO NK=1 ,WDES%NKPTS
      DO NB=1 ,WDES%NBANDS
         ION_SUM(:,:,ISP)=ION_SUM(:,:,ISP)+PAR(NB,NK,:,:,ISP)*RSPIN* &
              WDES%WTKPT(NK)*W%FERWE(NB,NK,ISP_)
!  OL : Orb mom part
         IF (WDES%LNONCOLLINEAR) &
           ION_SUM_lmom(:,:,ISP)=ION_SUM_lmom(:,:,ISP)+PAR_lmom(NB,NK,:,:,ISP)*RSPIN* &
                WDES%WTKPT(NK)*W%FERWE(NB,NK,ISP_)
      ENDDO
      ENDDO
      ENDDO

      ND=LPAR*T_INFO%NIONP
      IF (.NOT.WDES%LNONCOLLINEAR) THEN
           CALL R_FLIP(ION_SUM,ND,ND,WDES%NCDIJ,.FALSE.)
!  OL : Orb mom part           
!          CALL R_FLIP(ION_SUM_lmom,ND,ND,WDES%NCDIJ,.FALSE.)      
      ENDIF
      io_begin

      IF (LORBIT==1.OR.LORBIT==2) THEN
         ION_SUM2=0
         DO LL=0,LDIMP-1
            DO M=1,2*LL+1
               LM=LL*LL+M
               ION_SUM2(LL+1,:,:)=ION_SUM2(LL+1,:,:)+ION_SUM(LM,:,:)
            ENDDO
         ENDDO
      ELSE
         ION_SUM2=ION_SUM
      ENDIF
!=======================================================================
!   write PAR on file PROCAR
!=======================================================================
      IF (LORBIT /=5) THEN

      DO ISP=1,WDES%ISPIN

      WRITE(IUP,3200) WDES%NKPTS,WDES%NBANDS,T_INFO%NIONP
      DO NK=1,WDES%NKPTS
      WRITE(IUP,3201) NK,WDES%VKPT(1,NK),WDES%VKPT(2,NK),WDES%VKPT(3,NK),WDES%WTKPT(NK)
      DO NB=1,WDES%NB_TOT
      NI=1

      WRITE(IUP,3203) NB,REAL( W%CELTOT(NB,NK,ISP) ,KIND=q),2*W%FERTOT(NB,NK,ISP)
      WRITE(IUP,*)
      
      WRITE(IUP,'(A3)',ADVANCE='No') "ion"
      IF (LORBIT==1.OR.LORBIT==2) THEN
         DO NL=1,LPAR
            WRITE(IUP,'(A7)',ADVANCE='No') LMCHAR(NL)
         ENDDO
      ELSE
         DO NL=1,LPAR
            WRITE(IUP,'(A7)',ADVANCE='No') LCHAR(NL)
         ENDDO
      ENDIF
      WRITE(IUP,'(A7)',ADVANCE='yes') "tot"

      DO II=0,WDES%NRSPINORS*WDES%NRSPINORS-1
      PARSUM=0
      SUMION=0
      DO NI=1,T_INFO%NIONP
         S=0
         PARSUM=0
            DO NL=1,LPAR
               PARSUM=PARSUM+PAR(NB,NK,NL,NI,ISP+II)
               SUMION(NL)=SUMION(NL)+PAR(NB,NK,NL,NI,ISP+II)
               S=S+SUMION(NL)
            ENDDO
         WRITE(IUP,3204) NI,(PAR(NB,NK,NL,NI,ISP+II),NL=1,LPAR),PARSUM
      ENDDO
      IF (T_INFO%NIONP>1) THEN
         WRITE(IUP,3205) (SUMION(NL),NL=1,LPAR),S
      ENDIF
      ENDDO

      IF (LORBIT==2) THEN
         WRITE(IUP,'(A3)',ADVANCE='No') "ion"
         DO NL=1,LMDIMP
            WRITE(IUP,'(A7)',ADVANCE='No') LMCHAR(NL)
         ENDDO
         WRITE(IUP,*)

         DO NI=1,T_INFO%NIONP
            WRITE(IUP,3204) NI, (REAL (PHAS(M,NI,NK,NB,ISP)),M=1,LMDIMP)
            WRITE(IUP,3204) NI, (AIMAG(PHAS(M,NI,NK,NB,ISP)),M=1,LMDIMP)
         ENDDO
      ENDIF
      WRITE(IUP,*)

      ENDDO  ! bands
      ENDDO  ! kpoints
      ENDDO  ! spin

      ENDIF
!=======================================================================
!   set the LMTABLE
!=======================================================================
! allocate the LMTABLE
      IF (ALLOCATED(LMTABLE)) THEN
         DEALLOCATE(LMTABLE)
      ENDIF
      ALLOCATE(LMTABLE(LPAR, T_INFO%NIONP))
      LMTABLE="   "

      DO NI=1,T_INFO%NIONP
         IF (LORBIT==1.OR.LORBIT==2) THEN
            DO NL=1,LPAR
               LMTABLE(NL,NI)=LMCHAR(NL)
            ENDDO
         ELSE
            DO NL=1,LPAR
               LMTABLE(NL,NI)=LCHAR(NL)
            ENDDO
         ENDIF
      ENDDO
!=======================================================================
!   write condensed form of PAR on file OUTCAR
!=======================================================================
      DO ISP=1,WDES%NCDIJ
      WRITE(IU6,*)
      IF (IU6>=0) THEN
         SELECT CASE (ISP)
         CASE (1)
            WRITE(IU6,2011) 'total charge     '
         CASE (2)
            WRITE(IU6,2011) 'magnetization (x)'
         CASE (3)
            WRITE(IU6,2011) 'magnetization (y)'
         CASE (4)
            WRITE(IU6,2011) 'magnetization (z)'
         END SELECT
      ENDIF
 2011 FORMAT(//A18)
      WRITE(IU6,*)

      IF (LDIMP==4) THEN
         WRITE(IU6,12212)
      ELSE
         WRITE(IU6,2212)
      ENDIF

      NI=1
      PARSUM=0
      DO NL=1,LDIMP
         SUMION(NL)=ION_SUM2(NL,1,ISP)
         PARSUM    =PARSUM+SUMION(NL)
         SUMTOT(NL)=SUMION(NL)
      ENDDO
      WRITE(IU6,2213) NI,(SUMION(NL),NL=1,LDIMP),PARSUM
      IF (T_INFO%NIONS>1) THEN
        DO NI=2,T_INFO%NIONS
        S=0
        PARSUM=0
        DO NL=1,LDIMP
          SUMION(NL)=ION_SUM2(NL,NI,ISP)
          PARSUM    =PARSUM+SUMION(NL)
          SUMTOT(NL)=SUMTOT(NL)+SUMION(NL)
          S         =S  +SUMTOT(NL)
        ENDDO
        WRITE(IU6,2213) NI,(SUMION(NL),NL=1,LDIMP),PARSUM
        ENDDO
        WRITE(IU6,2215) (SUMTOT(NL),NL=1,LDIMP),S
      ENDIF

      ENDDO
!=======================================================================
!  OLebacq :   Print Orbital moments in OUTCAR
!=======================================================================      
      IF(WDES%LNONCOLLINEAR)  THEN
      
      DO ISP=1,WDES%NCDIJ-1
      WRITE(IU6,*)
      IF (IU6>=0) THEN
         SELECT CASE (ISP)
         CASE (1)
            WRITE(IU6,2011) 'orbital moment (x)'
         CASE (2)
            WRITE(IU6,2011) 'orbital moment (y)'
         CASE (3)
            WRITE(IU6,2011) 'orbital moment (z)'
         END SELECT
      ENDIF
      WRITE(IU6,*)

      IF (LDIMP==4) THEN
          WRITE(IU6,12212)
      ELSE
          WRITE(IU6,2212)
      ENDIF
         NI=1 
         PARSUM2=0 
         DO NL=1,LDIMP
            SUMION_lmom(NL)=ION_SUM_lmom(NL,1,ISP) 
            PARSUM2    =PARSUM2+SUMION_lmom(NL) 
            SUMTOT_lmom(NL)=SUMION_lmom(NL) 
         ENDDO

         WRITE(IU6,2216) &
             NI,(SUMION_lmom(NL),NL=1,LDIMP),PARSUM2 

         IF (T_INFO%NIONS>1) THEN   
             DO NI=2,T_INFO%NIONS 
                S2=0  
                PARSUM2=0
                DO NL=1,LDIMP   
                   SUMION_lmom(NL)=ION_SUM_lmom(NL,NI,ISP)  
                   PARSUM2    =PARSUM2+SUMION_lmom(NL)  
                   SUMTOT_lmom(NL)=SUMTOT_lmom(NL)+SUMION_lmom(NL) 
                   S2         =S2  +SUMTOT_lmom(NL) 
                ENDDO   
                WRITE(IU6,2216) &
                NI,(SUMION_lmom(NL),NL=1,LDIMP),PARSUM2
             ENDDO
             WRITE(IU6,2215) (SUMTOT_lmom(NL),NL=1,LDIMP),S2
          ENDIF
      ENDDO   
 2216 FORMAT(I3,5X,5(F8.3))
! OLebacq : end orbital moments calculation      

      ENDIF

      io_end

      IF (WDES%LNONCOLLINEAR) THEN
         DEALLOCATE(PAR_lmom)
      ENDIF


 2201 FORMAT(/' k-point ',I4,' :',3X,3F10.4,'     weight = ',F8.6/)
 2203 FORMAT(I3,5X,9(1X,F6.3),1X,F8.3)
 2204 FORMAT(I3,5X,9(1X,F6.3),1X,F8.3)
 2205 FORMAT(' ----------------------------------------------', &
     &          '-------------------------------------',/ &
     &          'tot',5X,9(1X,F6.3),1X,F8.3/)
 2212 FORMAT('# of ion     s       p ', &
     &       '      d       tot'/ &
     &          '----------------------------------------')
12212 FORMAT('# of ion     s       p ', &
     &       '      d       f       tot'/ &
     &          '------------------------------------------------')

 2213 FORMAT(I3,5X,5(F8.3))
 2215 FORMAT('------------------------------------------------',/ &
     &          'tot',5X,5(F8.2)/)

 3200 FORMAT('# of k-points:',I5,9X,'# of bands:',I4,9X,'# of ions:',I4)
 3201 FORMAT(/' k-point ',I4,' :',3X,3F11.8,'     weight = ',F10.8/)

 3203 FORMAT('band ',I3,' # energy',F14.8,' # occ.',F12.8)

 3204 FORMAT(I3,17(1X,F6.3))
 3205 FORMAT('tot',17(1X,F6.3)/)

      do_io CLOSE (IUP)
      IF (IU6>=0) WRITE(IU6,*)
      IF (LORBIT==2) DEALLOCATE(PHAS)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE SPHPRO_FAST************************
!
! fast partial density of states using the build in projector functions
! of the pseudopotentials
!
!***********************************************************************

      SUBROUTINE SPHPRO_FAST( &
          GRID,LATT_CUR,LATT_INI, P,T_INFO,W,WDES, IUP,IU6, &
          LOVERL,LMDIM,CQIJ,LPAR,LDIMP,LMDIMP,LFINAL, LORBIT,PAR)
      USE prec
      USE main_mpi
      USE constant
      USE wave
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE nonl

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR,LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (potcar)      BET(1)
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES_1K
      TYPE (nonl_struct) NONL_S

      LOGICAL LOVERL,LFINAL
      INTEGER LORBIT

      OVERLAP   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NRSPINORS*WDES%NRSPINORS)
      REAL(q) PAR(WDES%NB_TOT,WDES%NKDIM,LPAR, T_INFO%NIONP,WDES%NCDIJ)
!  allocate dynamically
      REAL(q) SUMION(LPAR),SUMTOT(LPAR)
      REAL(q) ION_SUM(LDIMP,T_INFO%NIONS,WDES%NCDIJ)
      REAL(q) ION_SUM_DETAIL(LMDIMP,T_INFO%NIONS,WDES%NCDIJ)
      COMPLEX(q) :: CSUM(WDES%NBANDS,WDES%NKPTS,WDES%NCDIJ)
      COMPLEX(q),ALLOCATABLE :: PHAS(:,:,:,:,:)
      CHARACTER(LEN=3) :: LCHAR(4)=(/'  s','  p','  d','  f'/)
      CHARACTER(LEN=3) :: LMCHAR(16)=(/'  s',' py',' pz',' px', &
      'dxy','dyz','dz2','dxz','dx2','f-3','f-2','f-1',' f0',' f1',' f2',' f3'/)

      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME= WDES%COMM%NODE_ME
      IONODE = WDES%COMM%IONODE
#endif
      RSPIN=WDES%RSPIN
!=======================================================================
! some initialization
!=======================================================================
      io_begin
!   write header fo file PROCAR and 

      IF (LFINAL) THEN
         OPEN(UNIT=IUP,FILE=DIR_APP(1:DIR_LEN)//'PROCAR',STATUS='UNKNOWN')
         IF (LORBIT==11) THEN
            WRITE(IUP,'(A)')'PROCAR lm decomposed'
         ELSE  IF (LORBIT==12) THEN
            WRITE(IUP,'(A)')'PROCAR lm decomposed + phase'
         ELSE
            WRITE(IUP,'(A)')'PROCAR new format'
         ENDIF
      ENDIF
      io_end

      IF (LFINAL) THEN
         PAR=0

         IF (LORBIT==12) THEN
            ALLOCATE(PHAS(LMDIMP,T_INFO%NIONS,WDES%NKDIM,WDES%NB_TOT,WDES%ISPIN))
            PHAS=0
         ENDIF
      ENDIF

      ION_SUM=0
      ION_SUM_DETAIL=0

      LMBASE =0
      NIS=1

      typ: DO NT=1,T_INFO%NTYP
      ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
!=======================================================================
      NIP=NI_LOCAL(NI, WDES%COMM_INB)
      IF (NIP==0) CYCLE ion ! not on local node

      LOW=1
      block: DO
      LL=P(NT)%LPS(LOW)
      DO LHI=LOW,P(NT)%LMAX
         IF (LL/=P(NT)%LPS(LHI)) EXIT
      ENDDO
      LHI=LHI-1

      MMAX=2*LL+1

      DO L =LOW,LHI
      DO LP=LOW,LHI

      DO M=1,MMAX
      LM=LL*LL+M

      IF (LFINAL .AND. (LORBIT==11.OR.LORBIT==12) .AND. LM > LPAR) THEN
         WRITE(0,*) 'internal ERROR: LPAR is too small in SPHPRO_FAST (LM)',LM,LPAR
         STOP
      ELSE IF (LFINAL .AND. LL+1 > LPAR) THEN
         WRITE(0,*) 'internal ERROR: LPAR is too small in SPHPRO_FAST (LL)',LL+1,LPAR
         STOP
      ENDIF

      IF (ASSOCIATED (P(NT)%QTOT) ) THEN

      DO ISP=1,WDES%ISPIN
      DO ISPINOR=0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMIND  =LMBASE +(L -LOW) *MMAX+M + ISPINOR *WDES%NPRO/2
      LMIND_ =LMBASE +(LP-LOW) *MMAX+M + ISPINOR_*WDES%NPRO/2
      II=ISP+ISPINOR_+2*ISPINOR

         DO NK=1,WDES%NKPTS
         DO NB=1,WDES%NBANDS
            CSUM(NB,NK,II)=  &
                 W%CPROJ(LMIND,NB,NK,ISP)*P(NT)%QTOT(LP,L)*GCONJG(W%CPROJ(LMIND_,NB,NK,ISP))
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ND=WDES%NBANDS*WDES%NKPTS

      IF (WDES%LNONCOLLINEAR) &
      CALL C_FLIP(CSUM,ND,ND,WDES%NCDIJ,.FALSE.)

      DO ISP=1,WDES%NCDIJ
      ISP_=MIN(ISP,WDES%ISPIN)
      DO NK=1 ,WDES%NKPTS

      CALL SETWDES(WDES,WDES_1K,NK)
      DO NB=1,WDES%NB_TOT
         NB_=NB_LOCAL(NB,WDES_1K)
         IF(NB_==0) CYCLE
         ION_SUM(LL+1,NI,ISP)=ION_SUM(LL+1,NI,ISP)+CSUM(NB_,NK,ISP)*RSPIN* &
              WDES%WTKPT(NK)*W%FERWE(NB_,NK,ISP_)
         ION_SUM_DETAIL(LM,NI,ISP)=ION_SUM_DETAIL(LM,NI,ISP)+CSUM(NB_,NK,ISP)*RSPIN* &
              WDES%WTKPT(NK)*W%FERWE(NB_,NK,ISP_)

         IF (LFINAL) THEN
          IF (LORBIT==11.OR.LORBIT==12) THEN 
            PAR(NB,NK,LM,NI,ISP)=PAR(NB,NK,LM,NI,ISP)+CSUM(NB_,NK,ISP)
          ELSE
            PAR(NB,NK,LL+1,NI,ISP)=PAR(NB,NK,LL+1,NI,ISP)+CSUM(NB_,NK,ISP)
          ENDIF
         ENDIF

      ENDDO
      ENDDO
      ENDDO
         
      ENDIF

      ENDDO
      ENDDO
      ENDDO

      IF (LORBIT==12 .AND. LFINAL) THEN
         DO ISP=1,WDES%ISPIN
         DO II=0,WDES%NRSPINORS-1
         DO L =LOW,LHI
         DO M=1,MMAX
         LM=LL*LL+M
         DO NK=1,WDES%NKPTS
         DO NB=1,WDES%NB_TOT
            NB_=NB_LOCAL(NB,WDES_1K)
            IF(NB_==0) CYCLE
            CTMP= W%CPROJ(LMBASE+M,NB_,NK,ISP+II)
            PHAS(LM,NI,NK,NB,ISP)=PHAS(LM,NI,NK,NB,ISP)+CTMP
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

!-----------------------------------------------------------------------
      LMBASE =LMBASE +(LHI-LOW+1)*MMAX
      LOW=LHI+1
      IF (LOW > P(NT)%LMAX) EXIT block
      ENDDO block

      ENDDO ion
      NIS = NIS+T_INFO%NITYP(NT)
      ENDDO typ

      CALLMPI( M_sum_d( WDES%COMM, ION_SUM, LDIMP*T_INFO%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_d( WDES%COMM, ION_SUM_DETAIL, LMDIMP*T_INFO%NIONS*WDES%NCDIJ))
      IF (LFINAL) THEN
        CALLMPI( M_sum_d( WDES%COMM, PAR(1,1,1,1,1), WDES%NB_TOT*WDES%NKDIM*LPAR*T_INFO%NIONP*WDES%NCDIJ))
        IF (LORBIT==12) THEN
           CALLMPI( M_sum_d( WDES%COMM, PHAS(1,1,1,1,1), LMDIMP*T_INFO%NIONS*WDES%NKDIM*WDES%NB_TOT*WDES%ISPIN*2))
        ENDIF
      ENDIF

      ND=LDIMP*T_INFO%NIONS
      IF (.NOT.WDES%LNONCOLLINEAR) &
           CALL R_FLIP(ION_SUM,ND,ND,WDES%NCDIJ,.FALSE.)
      io_begin
!=======================================================================
!   write PAR on file PROCAR
!=======================================================================
      IF (LFINAL) THEN

      DO ISP=1,WDES%ISPIN

      WRITE(IUP,3200) WDES%NKPTS,WDES%NB_TOT,T_INFO%NIONP
      DO NK=1,WDES%NKPTS
      WRITE(IUP,3201) NK,WDES%VKPT(1,NK),WDES%VKPT(2,NK),WDES%VKPT(3,NK),WDES%WTKPT(NK)
      DO NB=1,WDES%NB_TOT
      NI=1

      WRITE(IUP,3203) NB,REAL( W%CELTOT(NB,NK,ISP) ,KIND=q),2*W%FERTOT(NB,NK,ISP)
      WRITE(IUP,*)
      
      WRITE(IUP,'(A3)',ADVANCE='No') "ion"
      IF (LORBIT==11.OR.LORBIT==12) THEN
         DO NL=1,LPAR
            WRITE(IUP,'(A7)',ADVANCE='No') LMCHAR(NL)
         ENDDO
      ELSE
         DO NL=1,LPAR
            WRITE(IUP,'(A7)',ADVANCE='No') LCHAR(NL)
         ENDDO
      ENDIF
      WRITE(IUP,'(A7)',ADVANCE='yes') "tot"

      DO II=0,WDES%NRSPINORS*WDES%NRSPINORS-1
      PARSUM=0
      SUMION=0
      DO NI=1,T_INFO%NIONP
         S=0
         PARSUM=0
            DO NL=1,LPAR
               PARSUM=PARSUM+PAR(NB,NK,NL,NI,ISP+II)
               SUMION(NL)=SUMION(NL)+PAR(NB,NK,NL,NI,ISP+II)
               S=S+SUMION(NL)
            ENDDO
         WRITE(IUP,3204) NI,(PAR(NB,NK,NL,NI,ISP+II),NL=1,LPAR),PARSUM
      ENDDO
      IF (T_INFO%NIONP>1) THEN
         WRITE(IUP,3205) (SUMION(NL),NL=1,LPAR),S
      ENDIF
      ENDDO

      IF (LORBIT==12) THEN
         WRITE(IUP,'(A3)',ADVANCE='No') "ion"
         DO NL=1,LMDIMP
            WRITE(IUP,'(A7)',ADVANCE='No') LMCHAR(NL)
         ENDDO
         WRITE(IUP,*)

         DO NI=1,T_INFO%NIONP
            WRITE(IUP,3204) NI, (REAL (PHAS(M,NI,NK,NB,ISP)),M=1,LMDIMP)
            WRITE(IUP,3204) NI, (AIMAG(PHAS(M,NI,NK,NB,ISP)),M=1,LMDIMP)
         ENDDO
      ENDIF
      WRITE(IUP,*)

      ENDDO
      ENDDO
      ENDDO

      ENDIF
!=======================================================================
!   set the LMTABLE
!=======================================================================
      IF (LFINAL) THEN
         IF (ALLOCATED(LMTABLE)) THEN
            DEALLOCATE(LMTABLE)
         ENDIF
         ALLOCATE(LMTABLE(LPAR, T_INFO%NIONP))
         LMTABLE="   "

         DO NI=1,T_INFO%NIONP
            IF (LORBIT==11.OR.LORBIT==12) THEN
               DO NL=1,LPAR
                  LMTABLE(NL,NI)=LMCHAR(NL)
               ENDDO
            ELSE
               DO NL=1,LPAR
                  LMTABLE(NL,NI)=LCHAR(NL)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!=======================================================================
!   write condensed form of PAR on file OUTCAR
!=======================================================================
      DO ISP=1,WDES%NCDIJ
      WRITE(IU6,*)
      IF (IU6>=0) THEN
         SELECT CASE (ISP)
         CASE (1)
            WRITE(IU6,2011) 'total charge     '
         CASE (2)
            WRITE(IU6,2011) 'magnetization (x)'
         CASE (3)
            WRITE(IU6,2011) 'magnetization (y)'
         CASE (4)
            WRITE(IU6,2011) 'magnetization (z)'
         END SELECT
      ENDIF
 2011 FORMAT(//A18)
      WRITE(IU6,*)

      IF (LDIMP==4) THEN
         WRITE(IU6,12212)
      ELSE
         WRITE(IU6,2212)
      ENDIF

      NI=1
      PARSUM=0
      DO NL=1,LDIMP
        SUMION(NL)=ION_SUM(NL,1,ISP)
        PARSUM    =PARSUM+SUMION(NL)
        SUMTOT(NL)=SUMION(NL)
      ENDDO
      WRITE(IU6,2213) NI,(SUMION(NL),NL=1,LDIMP),PARSUM
      IF (T_INFO%NIONS>1) THEN
        DO NI=2,T_INFO%NIONS
        S=0
        PARSUM=0
        DO NL=1,LDIMP
          SUMION(NL)=ION_SUM(NL,NI,ISP)
          PARSUM    =PARSUM+SUMION(NL)
          SUMTOT(NL)=SUMTOT(NL)+SUMION(NL)
          S         =S  +SUMTOT(NL)
        ENDDO
        WRITE(IU6,2213) NI,(SUMION(NL),NL=1,LDIMP),PARSUM
        ENDDO
        WRITE(IU6,2215) (SUMTOT(NL),NL=1,LDIMP),S
      ENDIF

      ENDDO
      io_end


 2201 FORMAT(/' k-point ',I4,' :',3X,3F10.4,'     weight = ',F8.6/)
 2203 FORMAT(I3,5X,9(1X,F6.3),1X,F8.3)
 2204 FORMAT(I3,5X,9(1X,F6.3),1X,F8.3)
 2205 FORMAT(' ----------------------------------------------', &
     &          '-------------------------------------',/ &
     &          'tot',5X,9(1X,F6.3),1X,F8.3/)
 2212 FORMAT('# of ion     s       p ', &
     &       '      d       tot'/ &
     &          '----------------------------------------')
12212 FORMAT('# of ion     s       p ', &
     &       '      d       f       tot'/ &
     &          '------------------------------------------------')

 2213 FORMAT(I3,5X,5(F8.3))
 2215 FORMAT('------------------------------------------------',/ &
     &          'tot',5X,5(F8.3)/)


 3200 FORMAT('# of k-points:',I5,9X,'# of bands:',I4,9X,'# of ions:',I4)
 3201 FORMAT(/' k-point ',I4,' :',3X,3F11.8,'     weight = ',F10.8/)

 3203 FORMAT('band ',I3,' # energy',F14.8,' # occ.',F12.8)

 3204 FORMAT(I3,17(1X,F6.3))
 3205 FORMAT('tot',17(1X,F6.3)/)

 3207 FORMAT(' ion   s      py     pz    px     dxy    dyz    dz2    dxz  dx2-y2')

      do_io CLOSE (IUP)
      IF (IU6>=0) WRITE(IU6,*)
      IF (LORBIT==12 .AND. LFINAL) DEALLOCATE(PHAS)

      RETURN
      END SUBROUTINE

!*******************************************************************
!  SUBROUTINE BEZERO
!  searches for NMAX zeros in the sperical Bessel functions
!  i/o:
!         XNULL(NMAX) result
!         L           quantum number l
!  great full spaghetti code (writen by gK)
!********************************************************************

      SUBROUTINE BEZERO(XNULL,L,NMAX)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (STEP=.1_q, BREAK= 1E-10_q )
      DIMENSION XNULL(NMAX)
! initialization
      LBES = L
      X=STEP
      N=0
! entry point for next q_n
  30  CALL SBESSEL(X, BJ1, L)
! coarse search
  10  X=X+STEP
      CALL SBESSEL(X, BJ2, L)
! found one point
      IF(BJ1*BJ2 <0) THEN
        ETA=0.0_q
! intervall bisectioning
        SSTEP=STEP
        XX   =X
  20    SSTEP=SSTEP/2
        IF (BJ1*BJ2<0) THEN
          XX=XX-SSTEP
        ELSE
          XX=XX+SSTEP
        ENDIF
        CALL SBESSEL(XX, BJ2, L)
        IF (SSTEP>BREAK) GOTO 20

        N=N+1
        XNULL(N)=XX
        IF (N==NMAX) RETURN
        GOTO 30
      ENDIF
      GOTO 10

      END SUBROUTINE
      END MODULE msphpro


!************************ SUBROUTINE C_FLIP ***************************
!
! rearranges the storage mode of an array from 
! (up, down) (i.e. spinor representation) to (total,magnetization)
! also the reverse operation is possible if setting LBACK=.TRUE.
!
!***********************************************************************

      SUBROUTINE C_FLIP(C,NDIM,NELM,NCDIJ,LBACK)
      USE prec
      IMPLICIT NONE

      LOGICAL LBACK
      INTEGER NCDIJ,NDIM,NELM,N
      COMPLEX(q) :: C(NDIM,NCDIJ)
      REAL(q) FAC
      COMPLEX(q) :: CQU,CQD,C01,C10
      COMPLEX(q) :: C11,C00,CX,CY,CZ

      IF (NCDIJ==2 ) THEN
!=======================================================================
         FAC=1._q
         IF (LBACK) FAC=0.5_q
      
         DO N=1,NELM
            CQU=C(N,1)
            CQD=C(N,2)
            C(N,1)=FAC*(CQU+CQD)
            C(N,2)=FAC*(CQU-CQD)
         ENDDO

      ELSE IF ( NCDIJ==4 .AND. .NOT. LBACK) THEN
!=======================================================================
         DO N=1,NELM
            C00=C(N,1)
            C01=C(N,2)
            C10=C(N,3)
            C11=C(N,4)

            C(N,1)= C00+C11
            C(N,2)= C01+C10
            C(N,3)=(C01-C10)*(0._q,1._q)
            C(N,4)= C00-C11             
         ENDDO
      ELSE IF ( NCDIJ==4 .AND. LBACK) THEN
!=======================================================================
         FAC=0.5_q
         DO N=1,NELM
            C00=C(N,1)
            CX =C(N,2)
            CY =C(N,3)
            CZ =C(N,4)
            
            C(N,1)= (C00+CZ)*FAC
            C(N,2)= (CX-CY*(0._q,1._q))*FAC
            C(N,3)= (CX+CY*(0._q,1._q))*FAC
            C(N,4)= (C00-CZ)*FAC
         ENDDO
      ENDIF

      END SUBROUTINE
 

!************************ SUBROUTINE R_FLIP ***************************
!
! rearranges the storage mode of an array from 
! (up, down) (i.e. spinor representation) to (total,magnetization)
! also the reverse operation is possible if setting LBACK=.TRUE.
! (collinear version only)
!***********************************************************************

      SUBROUTINE R_FLIP(C,NDIM,NELM,NCDIJ,LBACK)
      USE prec
      IMPLICIT NONE

      LOGICAL LBACK
      INTEGER NCDIJ,NDIM,NELM,N
      REAL(q) :: C(NDIM,NCDIJ)
      REAL(q) FAC
      REAL(q) :: CQU,CQD

      IF (NCDIJ==2 ) THEN
         FAC=1._q
         IF (LBACK) FAC=0.5_q
      
         DO N=1,NELM
            CQU=C(N,1)
            CQD=C(N,2)
            C(N,1)=FAC*(CQU+CQD)
            C(N,2)=FAC*(CQU-CQD)
         ENDDO
      ENDIF

      END SUBROUTINE
      


!************************ SUBROUTINE SPHPRO_DESCRIPTION ****************
!
! this subroutine can be used to get a content description
! of the partial DOSCAR array
! it returns an empty string if the table is not yet set up
!
!***********************************************************************
    
      SUBROUTINE SPHPRO_DESCRIPTION(ni, l, lmtype)
      USE prec
      USE msphpro
      INTEGER ni, nis
      INTEGER l
      CHARACTER (LEN=3) :: lmtype

      IF (.NOT.ALLOCATED(LMTABLE)) THEN
        lmtype="err"
        RETURN
      ENDIF

      IF (ni > SIZE( LMTABLE, 2) .OR.  l > SIZE( LMTABLE, 1)) THEN
         lmtype="err"
         RETURN
      ENDIF


      IF (ni==0) THEN
         lmtype="   "
         DO nis=1,SIZE( LMTABLE, 2)
            IF (LMTABLE(l,nis) /= "   ") THEN
               lmtype=LMTABLE(l,nis)
               EXIT
            ENDIF
         ENDDO
      ELSE
         lmtype=LMTABLE(l, ni)
      ENDIF

      END SUBROUTINE
