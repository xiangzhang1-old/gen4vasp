!=======================================================================
!
!  this module writes the file that is required to evaluate
!  STM pictures using Bardeens approach
!  five parameters are supplied via the STM line in the INCAR file
!  STM(1)   minimum z position
!  STM(2)   maximum z position
!  STM(3)   delta z (recommended)
!           or refine (how many points are to be interpolated)
!  STM(4)   energy cutoff for G vectors  in ev
!           if hbar^2 G^2 / m_e < ABS (STM(4)) the corresponding G
!           vector is written to the file
!           if this is negative a window function is applied
!             psi(G) = psi(G) * 0.5_q*(Cos(PI*(G-GMIN)/(GMAX-GMIN))+1._q)
!             GMIN = 0 at present
!  STM(5)   energy window around fermi energy
!
!  a typical line in the INCAR file might be:
!  STM = 7.276 9.900 0.052918 -80 0.30
!  at present the STM program requires a spacing of 0.1 a.u. and 50 data
!  points
!  set STM(2) to STM(1)+49.5*STM(3)
!
!  we recommend to increase the cutoff ENMAX by 30 % (with respect to the
!  default one). 
!  PREC=Med is sufficient (although PREC= High yields slightly
!  better results).
!  We recommend a negative setting for STM(4) to remove any wiggles
!  in the corrugation.
!  
!
!=======================================================================

  SUBROUTINE WRT_STM_FILE(GRID, WDES, WUP, EFERMI, LATT_CUR, STM, T_INFO)
    USE prec
    USE lattice
    USE poscar
    USE mpimy
    USE mgrid
    USE wave
    USE constant

    IMPLICIT NONE

    TYPE (type_info)   T_INFO
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (wavedes)     WDES       ! description of wavefunctions
    TYPE (wavefun)     WUP        ! wavefunction (up)
    REAL(q)            EFERMI     ! fermi energy
    TYPE (latt)        LATT_CUR

    REAL(q) :: STM(5)             ! all stm paramters
    REAL(q) :: A(3,3)             ! lattice

    ! local

    COMPLEX(q),ALLOCATABLE :: CWORK(:,:,:)
    COMPLEX(q),ALLOCATABLE :: C1(:),C2(:),CWRK(:)
    REAL(q), ALLOCATABLE :: PREAL(:,:),PCMPLX(:,:)
    REAL(q) :: ENMAX_STM          ! energy cutoff for STM
    LOGICAL :: LWINDOW            ! use a window function before FFT
    INTEGER :: REFINE             ! how many points must be interpolated
    REAL(q) :: DELTAZ             ! distance between the points z direction
    REAL(q) :: DELTAE             ! energy interval around Fermi energy
    INTEGER,ALLOCATABLE :: COARSE_TO_FINE(:)
    REAL(q),ALLOCATABLE :: GZ(:)              ! 2 Pi G/ || a(3) ||
    LOGICAL,ALLOCATABLE :: DO_IT(:,:)
    REAL(q) :: EIGENVAL
    INTEGER :: NBANDS_WRITE       ! how many bands must be written for one k-point
    INTEGER :: MAX_BANDS
    INTEGER :: MAX_NG_WRITE
    INTEGER :: NG_WRITE           ! how many G vectors do we have
    INTEGER :: N1,N2,N3,N,I,NZ_START,NZ_END,NK,II
    REAL(q) :: ZSTART,ZEND,Z,D,WSCALE,WR,WI
    REAL(q) :: GIX,GIY,GIZ,G1,G2,G3,ENERGI,SCALE, GMAX, GMIN
    REAL(q) :: DISTANCE
    ! fft 
    INTEGER :: NFFT,NI
    REAL(q),ALLOCATABLE:: TRIG(:)
    INTEGER IFAC(19)
    LOGICAL, EXTERNAL :: FFTCHK_FURTH
    ! integer
    INTEGER, PARAMETER :: IU=77
    INTEGER ::  NWRT = 0

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------
 
    ! no STM output required, do quick return
    IF (STM(1) >= STM(2) .OR. STM(3) == 0 ) THEN
    RETURN
    ENDIF
    ! start counting the number of loops
    NWRT = NWRT + 1

    ! set the energy cutoff for STM calculation
    DELTAZ   =0
    REFINE   =STM(3)
    ENMAX_STM=ABS(STM(4))
    IF (STM(4)<0) THEN
       LWINDOW=.TRUE.
    ELSE
       LWINDOW=.FALSE.
    ENDIF
    DELTAE   =STM(5)

    ! set the scale for the amplitudes
    WSCALE = LATT_CUR%A(3,3) / LATT_CUR%OMEGA
    WSCALE = WSCALE * AUTOA**1.5

    IF (REFINE==0) THEN
       DELTAZ=STM(3)
       ZSTART=MAX(STM(1),0._q)
       ZEND  =MIN(STM(2),LATT_CUR%ANORM(3))
       IF (DELTAZ==0) DELTAZ=0.2
       NFFT=LATT_CUR%ANORM(3)/DELTAZ
    ELSE
       DELTAZ=0
       NFFT=REFINE*GRID%NGZ
    ENDIF


    DO 
       IF (FFTCHK_FURTH(NFFT)) THEN
	EXIT
	ENDIF
       NFFT=NFFT+1
    ENDDO

    ! allocate required work arrays
    ALLOCATE( DO_IT(GRID%NGX, GRID%NGY), &
              CWORK(GRID%NGZ, GRID%NGX, GRID%NGY), &
              C1(NFFT), C2(NFFT), CWRK(20*NFFT), TRIG(2*NFFT), &
              COARSE_TO_FINE(GRID%NGZ),GZ(GRID%NGZ), PREAL(NFFT,5), PCMPLX(NFFT,5) )

    
    CALL CFTTAB(NFFT,IFAC,TRIG)

    ! setup gradient
    DO I=1,GRID%NGZ
       COARSE_TO_FINE(I)=MOD(GRID%LPCTZ(I)+NFFT,NFFT)+1
       GZ(I)=GRID%LPCTZ(I)*TPI/LATT_CUR%ANORM(3)
    ENDDO

    ! figure out which z-values have to be written to the file
    IF ( DELTAZ==0 ) THEN
       NZ_START=    MIN(MAX(INT(STM(1)/LATT_CUR%ANORM(3)*NFFT+1.5),1),NFFT)
       NZ_END  =    MIN(MAX(INT(STM(2)/LATT_CUR%ANORM(3)*NFFT+1.5),1),NFFT)
    ELSE
       NZ_START=0
       NZ_END=(ZEND-ZSTART)/DELTAZ
    ENDIF
!-----------------------------------------------------------------------
! open the file and write the most important information
! lattice, maximal number of bands, number of kpoints
! only write for magnetic systems if we are in first round
!-----------------------------------------------------------------------
    IF (NWRT.EQ.1) &
     OPEN(UNIT=IU,FILE="STM",STATUS="UNKNOWN",POSITION="APPEND")

    A= LATT_CUR%A
    ! A(3,3) is defined to be the distance from between the
    ! topmost atom and the z-starting value
    DISTANCE=1

    DO NI=1,T_INFO%NIONS
       DISTANCE=MIN(ABS(T_INFO%POSION(3,NI)- ZSTART/LATT_CUR%ANORM(3)),DISTANCE)
    ENDDO
    
    A(3,3) = DISTANCE*LATT_CUR%ANORM(3)
    A=A / AUTOA

    IF (NWRT.EQ.1) THEN
     WRITE(IU,10) WSCALE
     WRITE(IU,11) A
    END IF

    ! count the number of bands
    MAX_NG_WRITE=0
    MAX_BANDS=0
    DO NK=1,WDES%NKPTS
       NBANDS_WRITE=0
       DO N=1,WDES%NBANDS
          EIGENVAL=WUP%CELEN(N,NK)
          IF ( EIGENVAL> EFERMI-DELTAE .AND. EIGENVAL< EFERMI+DELTAE) THEN
             NBANDS_WRITE=NBANDS_WRITE+1
          ENDIF
       ENDDO
       MAX_BANDS=MAX( MAX_BANDS, NBANDS_WRITE)

       NG_WRITE=1
       DO N1=1,GRID%NGX
       DO N2=1,GRID%NGY
          N3=1
          
          G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
          G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))

          GIX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)) *TPI
          GIY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)) *TPI
          GIZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)) *TPI

          ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
            
          IF(ENERGI< ENMAX_STM) THEN
             NG_WRITE=NG_WRITE+1
          ENDIF
       ENDDO
       ENDDO
       MAX_NG_WRITE=MAX(MAX_NG_WRITE,NG_WRITE)

    ENDDO

    IF (NWRT.EQ.1) &
     WRITE(IU,12) EFERMI/RYTOEV/2,WDES%ISPIN, &
                  WDES%NKPTS, NZ_END+1-NZ_START, MAX_NG_WRITE, MAX_BANDS 

    kpoints: DO NK=1,WDES%NKPTS
!-----------------------------------------------------------------------
! first find out the grid points G(x,y) for which the
! wavefunction has to be written to the file
!-----------------------------------------------------------------------
       NG_WRITE=0

       DO N1=1,GRID%NGX
       DO N2=1,GRID%NGY
          DO_IT(N1,N2)=.FALSE.
          N3=1
          
          G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
          G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
          G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))

          GIX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)) *TPI
          GIY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)) *TPI
          GIZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)) *TPI

          ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))
            
          IF(ENERGI< ENMAX_STM) THEN
             DO_IT(N1,N2)=.TRUE.
             NG_WRITE=NG_WRITE+1
          ENDIF
       ENDDO
       ENDDO
!-----------------------------------------------------------------------
! now for all bands that are within a certain interval around the Fermi
! level:
! extract the wavefunction and perform an FFT in z direction
!-----------------------------------------------------------------------
     ! first count the number of bands that are inside a certain intervall
     ! of the Fermi-level
     NBANDS_WRITE=0
     DO N=1,WDES%NBANDS
        EIGENVAL=WUP%CELEN(N,NK)
        IF ( EIGENVAL> EFERMI-DELTAE .AND. EIGENVAL< EFERMI+DELTAE) THEN
           NBANDS_WRITE=NBANDS_WRITE+1
        ENDIF
     ENDDO

     WRITE(IU,1) NK, NBANDS_WRITE, NG_WRITE, WDES%VKPT(1:2,NK),WDES%WTKPT(NK)
     
     bands:   DO N=1,WDES%NBANDS
       EIGENVAL=WUP%CELEN(N,NK)
       IF ( EIGENVAL> EFERMI-DELTAE .AND. EIGENVAL< EFERMI+DELTAE) THEN

       WRITE(IU,2) EIGENVAL/RYTOEV/2,WUP%FERWE(N,NK)*WDES%WTKPT(NK)
      
       ! transfer the points from the array to a 3d grid
       CWORK=0
       DO I=1,WDES%NPLWKP(NK)
          N1=MOD( WDES%IGX(I,NK)+ GRID%NGX, GRID%NGX)+1
          N2=MOD( WDES%IGY(I,NK)+ GRID%NGY, GRID%NGY)+1
          N3=MOD( WDES%IGZ(I,NK)+ GRID%NGZ, GRID%NGZ)+1

          CWORK(N3, N1, N2)=WUP%CPTWFP(I, N, NK)
       ENDDO

       ! now perform an fft of the relevant data
       ! and writen them to the file
       DO N1=1,GRID%NGX
       DO N2=1,GRID%NGY
          IF (DO_IT(N1,N2)) THEN
             C1=0
             C2=0

             WRITE(IU,3) GRID%LPCTX(N1),GRID%LPCTY(N2)


             DO N3=1,GRID%NGZ
                C1(COARSE_TO_FINE(N3))=CWORK(N3, N1, N2)
                C2(COARSE_TO_FINE(N3))=CWORK(N3, N1, N2)*GZ(N3)*(0._q,1._q)
             ENDDO

             IF (LWINDOW) THEN
                GMAX=SQRT(WDES%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
                GMIN=0
                DO N3=1,GRID%NGZ
                   I=GRID%LPCTZ(N3)
                   IF ( ABS(I) < GMIN) THEN
                      SCALE=1
                   ELSE IF ( ABS(I) > GMAX) THEN
                      SCALE=0
                   ELSE
                      SCALE=0.5_q*(COS(PI*(ABS(I)-GMIN)/(GMAX-GMIN))+1._q)
                   ENDIF
                   C1(COARSE_TO_FINE(N3))=CWORK(N3, N1, N2)*SCALE
                   C2(COARSE_TO_FINE(N3))=CWORK(N3, N1, N2)*GZ(N3)*(0._q,1._q)*SCALE
                ENDDO
             ENDIF

             CALL FFT_ONE(C1, CWRK, NFFT, TRIG, IFAC)
             CALL FFT_ONE(C2, CWRK, NFFT, TRIG, IFAC)

             IF (DELTAZ==0) THEN
                WRITE(IU,4) C1(NZ_START:NZ_END)
!               WRITE(IU,*)
!               WRITE(IU,4) C2(NZ_START:NZ_END)
!                DO I=NZ_START,NZ_END
!                   WRITE(IU,'(3F14.7)') (I-1)*LATT_CUR%ANORM(3)/NFFT,C1(I)
!                ENDDO
!                WRITE(IU,*)
!                DO I=NZ_START,NZ_END
!                   WRITE(IU,'(5F14.7)') (I-1)*LATT_CUR%ANORM(3)/NFFT,C2(I), &
!                        (C1(MIN(I+1,NFFT))-C1(MAX(1,I-1)))/(2*LATT_CUR%ANORM(3)/NFFT)
!                ENDDO
             ELSE
                DO I=1,NFFT
                   PREAL (I,1)=(I-1)*LATT_CUR%ANORM(3)/NFFT
                   PCMPLX(I,1)=(I-1)*LATT_CUR%ANORM(3)/NFFT
                   PREAL (I,2)=REAL(C1(I))
                   PCMPLX(I,2)=AIMAG(C1(I))
                ENDDO
                CALL SPLCOF(PCMPLX,NFFT,NFFT,1E30_q)
                CALL SPLCOF(PREAL ,NFFT,NFFT,1E30_q)

                DO II=0,NZ_END-NZ_START
                   Z=ZSTART+II*DELTAZ
                   I=MIN(INT(Z/LATT_CUR%ANORM(3)*NFFT+1),NFFT-1)
                   D=Z-PREAL(I,1)
                   WR = ((PREAL(I,5)*D+PREAL(I,4))*D+PREAL(I,3))*D+PREAL(I,2)
                   WI = ((PCMPLX(I,5)*D+PCMPLX(I,4))*D+PCMPLX(I,3))*D+PCMPLX(I,2)
                   WRITE (IU,4) WR * WSCALE, WI * WSCALE
!                   WRITE(IU,4) ((PREAL(I,5)*D+PREAL(I,4))*D+PREAL(I,3))*D+PREAL(I,2), &
!                                ((PCMPLX(I,5)*D+PCMPLX(I,4))*D+PCMPLX(I,3))*D+PCMPLX(I,2)
!                   WRITE(IU,'(5F14.7)') Z,((PREAL(I,5)*D+PREAL(I,4))*D+PREAL(I,3))*D+PREAL(I,2), &
!                        ((PCMPLX(I,5)*D+PCMPLX(I,4))*D+PCMPLX(I,3))*D+PCMPLX(I,2)
                ENDDO
                      
             ENDIF
          ENDIF
       ENDDO
       ENDDO
 
       
       ENDIF
    ENDDO bands
    ENDDO kpoints


!-----------------------------------------------------------------------
! clean up the mess
!-----------------------------------------------------------------------
    DEALLOCATE( DO_IT, CWORK, C1, C2, CWRK, TRIG, COARSE_TO_FINE, GZ )

1   FORMAT('k-point ',I6,' bands ',I6,'  G-vectors ',I6/, &
           'k-point ',3F16.10)
2   FORMAT('eigenenergy ',F16.10,'  occupancy ',F16.10)
3   FORMAT('G-vector: ',2I6)
4   FORMAT('(',E13.5,',',E13.5,')')

10  FORMAT('Scale for VASP output:',F14.7,/)
11  FORMAT((3F14.7))
12  FORMAT('fermi-energy:',F14.7,/ &
           'ispin:',I4,'  k-points: ',I4,'  z-values: ',I4,'  G-vectors: ',I4,'  max-eigenval',I4/)
  END SUBROUTINE WRT_STM_FILE


!=======================================================================
! one dimensional FFT
!=======================================================================

  SUBROUTINE FFT_ONE(C1, CWRK, NFFT, TRIG, IFAC)
    USE prec
    IMPLICIT NONE
    INTEGER NFFT
    REAL(q) C1(2*NFFT)
    REAL(q) CWRK(4*NFFT)
    REAL(q) TRIG(2*NFFT)
    INTEGER IFAC(19)
    INTEGER ISIGN

    ISIGN=+1
    CALL  CFFTML(C1(1),C1(2),CWRK(1),TRIG(1),IFAC(1),2,2*NFFT,NFFT,1,ISIGN,4*NFFT)

  END SUBROUTINE FFT_ONE
