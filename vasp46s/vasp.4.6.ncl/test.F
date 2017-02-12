#include "symbol.inc"

#define WRITE1 IF(IONODE == NODE_ME) WRITE(0,*)

PROGRAM TEST
  USE prec
  USE mpimy
  USE jacobi
  USE scala
  IMPLICIT NONE
  TYPE (communic)    :: COMM
  INTEGER, PARAMETER :: ND=20 ! ND=320
  REAL(q), PARAMETER :: PERT=10.0,SPC=4.0
  GDEF :: CH(ND,ND), CU(ND,ND), CA(ND,ND)
  INTEGER, PARAMETER :: LWORK=20
  COMPLEX(q) :: CTMP(LWORK*ND)
  REAL(q)    :: R(ND),RWORK(3*ND)

  INTEGER NODE_ME,IONODE,N,N1,N2,INFO,IFAIL
  REAL(q) :: TV,TV0,TC,TC0
  REAL(q), EXTERNAL  :: RANE
  INTEGER, EXTERNAL  :: MYID,MROOTID
! init mpi
   CALL INITPM(COMM)
   CALL INITCOMM(COMM)
   NODE_ME= MYID()
   IONODE = MROOTID()
   WRITE1 "hello"
!-----------------------------------------------------------------
! setup some diagonal dominant matrix
!-----------------------------------------------------------------

  DO N1=1,ND
  DO N2=1,N1-1
!   CH(N1,N2)=CMPLX(RANE(),RANE())*PERT
    CH(N1,N2)=RANE()*PERT
    CH(N2,N1)=GCONJG(CH(N1,N2))
  ENDDO
  ENDDO
  WRITE1 "random done"
  DO N=1,ND
    CH(N,N)=N+SPC*RANE()
  ENDDO
! some positive devinite matrix

  CU=MATMUL(TRANSPOSE(GCONJG(CH)),CH)
  WRITE1 "matmul done"
  IF (IONODE == NODE_ME) CALL MATDUMP(CH,ND)

!-----------------------------------------------------------------
! time gsum 4 times
!-----------------------------------------------------------------

  CA=CH
  CALL VTIME(TV0,TC0)
  CALLMPI( M_sum_g(COMM,CA, ND*ND))
  CALLMPI( M_sum_g(COMM,CA, ND*ND))
  CALLMPI( M_sum_g(COMM,CA, ND*ND))
  CALLMPI( M_sum_g(COMM,CA, ND*ND))

  CALL VTIME(TV,TC);  WRITE1 "gsum",  TV-TV0
!-----------------------------------------------------------------
! time M_alltoall_d 4 times
!-----------------------------------------------------------------
  CALL VTIME(TV0,TC0)
  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
  CALL VTIME(TV,TC);  WRITE1 "gsum_inter",  TV-TV0
!-----------------------------------------------------------------
! time  DPOTRF scaLAPACK
!-----------------------------------------------------------------

  CA=CU
  CALL VTIME(TV0,TC0)

  #ifdef gammareal
     CALL DPOTRF &
  #else
     CALL ZPOTRF &
  #endif
      & ('U',ND,CA(1,1),ND,INFO)
     IF (INFO/=0) THEN
        WRITE(*,*) 'LAPACK: Routine ZPOTRF failed!',INFO
        STOP
     ENDIF
  #ifdef gammareal
     CALL DTRTRI &
  #else
     CALL ZTRTRI &
  #endif
      & ('U','N',ND,CA(1,1),ND,INFO)
     IF (INFO/=0) THEN
        WRITE(*,*) 'LAPACK: Routine ZTRTRI failed!',INFO
        STOP
     ENDIF
  CALL VTIME(TV,TC);  WRITE1 "POTRF_TRTRI",  TV-TV0

! time POTRF_TRITRI scalapack

  CA=CU
  CALL VTIME(TV0,TC0)
  CALL pPOTRF_TRTRI (COMM, CA,ND)
  CALL VTIME(TV,TC);  WRITE1 "scaPOTRF_TRTRI",  TV-TV0
  IF (IONODE == NODE_ME) CALL MATDUMP(CA,ND)

!-----------------------------------------------------------------
! time ZHEEV
!-----------------------------------------------------------------

  CA=CH
  CALL VTIME(TV0,TC0)
  #ifdef  gammareal
      CALL DSYEV &
     &  ('V','U',ND,CA(1,1),ND,R, &
     &           CTMP,LWORK*ND, IFAIL)
  #else
      CALL ZHEEV &
     &  ('V','U',ND,CA(1,1),ND,R, &
     &           CTMP,LWORK*ND, RWORK,  IFAIL)
  #endif
      IF (IFAIL/=0) THEN
         WRITE(*,*) 'ERROR EDDIAG: Call to routine ZHEEV failed! '// &
     &              'Error code was ',IFAIL
         STOP
      ENDIF

  CALL VTIME(TV,TC);  WRITE1 "ZHEEV",  TV-TV0


  IF (NODE_ME==IONODE) WRITE(0,'(8F10.4)') R
  IF (IONODE == NODE_ME) CALL MATDUMP(CA,ND)

!-----------------------------------------------------------------
! time ZHEEV scalapack
!-----------------------------------------------------------------

  CA=CH
  CALL VTIME(TV0,TC0)
  CALL pDSSYEX_ZHEEVX(COMM, CA, R, ND)
  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
  CALL VTIME(TV,TC);  WRITE1 "scaZHEEV",  TV-TV0

  IF (NODE_ME==IONODE) WRITE(0,'(8F10.4)') R
  IF (IONODE == NODE_ME) CALL MATDUMP(CA,ND)
  STOP
!-----------------------------------------------------------------
! time ZHEEV jacobi
!-----------------------------------------------------------------

  CA=CH
  CALL VTIME(TV0,TC0)
!  CALL jacDSSYEV(COMM, CA, R, ND)
!  CALLMPI( M_sumf_g(COMM,CA, ND*ND))
!  CALLMPI( M_sumf_g(COMM,R , ND))

!  IF (NODE_ME==IONODE) WRITE(0,'(8F10.4)') R
!  IF (IONODE == NODE_ME) CALL MATDUMP(CA,ND)
  CALL VTIME(TV,TC);  WRITE1 "jacZHEEV",  TV-TV0

  STOP


CONTAINS
   SUBROUTINE MATDUMP(CH,ND)
   USE prec
   USE mpimy
   IMPLICIT NONE
   INTEGER :: ND
   GDEF :: CH(ND,ND)
   INTEGER :: NODE_ME,IONODE,N1,N2,N
  INTEGER, EXTERNAL  :: MYID,MROOTID

   NODE_ME= MYID()
   IONODE = MROOTID()

      io_begin
        N=MIN(10,ND)
        DO N1=1,N
          WRITE(6,1)N1,(REAL( CH(N1,N2) ,KIND=q) ,N2=1,N)
        ENDDO
         WRITE(6,*)
     #ifndef gammareal
        DO N1=1,N
          WRITE(6,2)N1,(AIMAG(CH(N1,N2)),N2=1,N)
        ENDDO
        WRITE(6,*)
     #endif
      1 FORMAT(1I2,3X,20F9.5)
      2 FORMAT(1I2,3X,20E9.1)
      io_end

   END SUBROUTINE
END PROGRAM
