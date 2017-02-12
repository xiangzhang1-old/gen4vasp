!*******************************************************************
!
! due to a design fault presently the PAW potentials to not
! contain information about the local pseudopotential in
! real space 
! this subroutine extracts the required information
! by means of a sine transform of the local pseudopotential
! from reciprocal to real space
!
! the subroutine calculates either the charge density rho(r) 
! or the potential V(r) 
! that corresponds to the local pseudopotential V(q)
! i.e.
!
! rho(r) = 1/r/(2 pi^2) \int_0^Infty sin(qr) q^2/ ( 4 Pi e^2) V(q) q dq
! r V(r) = 1/  (2 pi^2) \int_0^Infty sin(qr) V(q) q dq
! 
! on entry: 
! Z valence
! V stores the potential
!  V(q) = VQ -  8 \pi Z/q^2
!
!*******************************************************************

      SUBROUTINE POTTORHO( Z, NQL, VQ, DELQL, LPOT , NMAX, R, POTPSC ) 
      USE prec
      USE constant
      IMPLICIT NONE
      
      REAL(q) Z          ! valence
      INTEGER NQL        ! number of grid points in local potential          
      REAL(q) VQ(NQL)    ! local potential
      REAL(q) DELQL      ! step width for the array VQ
      LOGICAL LPOT       ! transform to get potential
                         ! or charge
      INTEGER NMAX       ! number of grid points
      REAL(q) R(NMAX)    ! r of each grid point
      REAL(q) POTPSC(NMAX)! potential
! local variable
      INTEGER, PARAMETER :: NFFT=32768
      REAL(q) RMAXF      ! maximum r values after FFT
      REAL(q) RDEL       ! spacing in real space
      REAL(q) FFT(NFFT+2)! array for FFT
      REAL(q) SPL(NFFT/2,5)
                         ! spline array for potential/charge in real space
      REAL(q) QQ

      INTEGER K
      REAL(q) :: FAKT,DUMMY,RHO,RR,ALPC,ALP,ALPSQRT
      REAL(q) :: ERRF
      REAL(q) :: FE=2
      INTEGER I
      
      RMAXF= 2*PI/DELQL

      ALPC=-(NQL*DELQL)**2/4/LOG(1.D-1)
      ALP=1/AUTOA**2

! first set up the values in the reciprocal space array

      DO I=1,NFFT
         QQ=DELQL*(I-1)
! fft of charge
         IF (I<=NQL .AND. LPOT) THEN
!            FFT(I)=VQ(I)*QQ
! fft of potential
! 
! rho(q) = q^2/ ( 4 Pi e^2) V(q) 
            IF (I/=1) THEN
              FFT(I)=(VQ(I)-Z*(1-EXP(-QQ*QQ/4/ALP))*4*PI*FELECT/ &
     &              (QQ*QQ))*QQ !*EXP(-QQ*QQ/4/ALPC)
            ELSE
              FFT(I)=0
            ENDIF
         ELSE IF (I<=NQL) THEN
!            FFT(I)=VQ(I)*QQ
! fft of potential
! the charge is given by
! rho(q) = q^2/ ( 4 Pi e^2) V(q) 
            FFT(I)=(VQ(I)*(QQ*QQ/(4*PI*FELECT))-Z)*QQ*EXP(-QQ*QQ/4/ALPC)
         ELSE
            FFT(I)=0
         ENDIF
      ENDDO
      
      CALL REALFT(FFT,NFFT/2,1)

      RDEL = RMAXF/NFFT
! f(r) = int sin_0^Infty  f(q) sin qr dq = del q sum_0^NFFT  f(q) sin qr

      DO K=1,NFFT/2
         SPL(K,2)=FFT(K*2)*DELQL
      ENDDO
! setup linear r mesh after FFT
      DO K=1,NFFT/2
         SPL(K,1)=(K-1)*RDEL
      ENDDO
!
! inverse transform of f(q) = 4 pi / q       \int sin (qr) f(r) r dr
! is                   f(r) = 1 /(2 pi^2) /r \int sin (qr) f(r) q dr 
! since 2 / pi \int_0^Infty sin(qr) dr =  delta(q)
      DO K=2,NFFT/2
         SPL(K,2)=SPL(K,2)/SPL(K,1)/2/PI/PI
      ENDDO

      IF (LPOT) THEN
         ALPSQRT = SQRT(ALP)
         DO K=2,NFFT/2
            SPL(K,2)=SPL(K,2)-FELECT*Z*ERRF(ALPSQRT*SPL(K,1))/SPL(K,1)
         ENDDO
         SPL(1,1)=SPL(1,1)-FELECT*Z*2*ALPSQRT/SQRT(PI)
         
      ENDIF
!
! now we need to go from the linear grid to a logarithmic one
!
      CALL SPLCOF(SPL ,NFFT/2, NFFT/2, 1.E30_q)

      POTPSC=0
      DO K=1,NMAX
         IF (R(K)/RDEL>NFFT) EXIT
         CALL SPLVAL(R(K), RHO, DUMMY, SPL, NFFT/2, NFFT/2)
         IF (LPOT) THEN
            POTPSC(K)=RHO
         ELSE
         ! set radially integrated charge
            POTPSC(K)=RHO*4*PI*R(K)**2
         ENDIF
      ENDDO
!      WRITE(95,"(2F14.7)") (R(K), POTPSC(K)*R(K)/RYTOEV/AUTOA,K=1,NMAX)

      END


!*******************************************************************
!
! the following subroutines are from Numerical Recepies
! they perform real to complex transforms
! or sinus transforms
!
!*******************************************************************

      SUBROUTINE SINFT(Y,N)
      USE prec
      IMPLICIT NONE
      REAL(q) WR,WI,WPR,WPI,WTEMP,THETA,Y1,Y2,SUM
      INTEGER N,M,J
      REAL(q) Y(N)
      THETA=3.14159265358979D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.0
      M=N/2
      DO 11 J=1,M
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
        Y1=WI*(Y(J+1)+Y(N-J+1))
        Y2=0.5*(Y(J+1)-Y(N-J+1))
        Y(J+1)=Y1+Y2
        Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL REALFT(Y,M,+1)
      SUM=0.0
      Y(1)=0.5*Y(1)
      Y(2)=0.0
      DO 12 J=1,N-1,2
        SUM=SUM+Y(J)
        Y(J)=Y(J+1)
        Y(J+1)=SUM
12    CONTINUE
      RETURN
      END



      SUBROUTINE REALFT(DATA,N,ISIGN)
      USE prec
      IMPLICIT NONE
      REAL(q) WR,WI,WPR,WPI,WTEMP,THETA,WRS,WIS
      REAL(q) H2R,H2I,H1R,H1I,C1,C2
      REAL(q) DATA(*)
      INTEGER N
      INTEGER I3,I4,I1,I2,N2P3,ISIGN,I
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
        DATA(2*N+1)=DATA(1)
        DATA(2*N+2)=DATA(2)
      ELSE
        C2=0.5
        THETA=-THETA
        DATA(2*N+1)=DATA(2)
        DATA(2*N+2)=0.0
        DATA(2)=0.0
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      N2P3=2*N+3
      DO 11 I=1,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        DATA(2)=DATA(2*N+1)
      ELSE
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      USE prec
      IMPLICIT NONE
      REAL(q) WR,WI,WPR,WPI,WTEMP,THETA,TEMPI,TEMPR
      REAL(q) DATA(*)
      INTEGER M,ISTEP,MMAX,ISIGN,NN,N,I,J
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END


