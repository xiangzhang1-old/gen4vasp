#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: chain.F,v 1.2 2002/08/14 13:59:37 kresse Exp $
!
! Module which implements the elastic band and the nudged
! elastic band method (for references see below)
! module becomes active if IMAGES tag is read from the INCAR files
!
!**********************************************************************

  MODULE chain
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    IMPLICIT NONE

    REAL(q),ALLOCATABLE,SAVE :: posion_all(:,:,:)
    REAL(q) :: spring=10

!**********************************************************************
!
!  routine for forces between the images on the elastic band
!
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions,posion,toten,force,a,b,iu6)
      INTEGER :: nions
      INTEGER :: iu6
      REAL(q) :: posion(3,nions)
      REAL(q) :: force(3,nions),toten
      REAL(q) :: tangent(3,nions)
      REAL(q) :: a(3,3),b(3,3)
! local variables
      REAL(q) :: x(3),x1(3),x2(3)
      REAL(q) :: e_chain,e_image,norm,norm1,norm2,proj,proj2,d1,d2
      REAL(q) :: force_chain(3,nions),d,const,ftangent
      INTEGER node,i,ni

      IF (images==0) RETURN
      IF (spring==-1000) RETURN
#if defined(MPI) || defined (MPI_CHAIN)
!
! communicate all positions to all nodes (brute force, but simple)
!
      node=comm_chain%node_me

      posion_all(:,:,1:images)=0
      posion_all(:,:,node)    =posion
      CALLMPI_C( M_sum_d( comm_chain, posion_all(1,1,1), nions*3*images))

!======================================================================
!
!  calculate the tangent in each point
!  and the distance to left and right image
!======================================================================
      i=node
      norm=0
      d1=0
      d2=0

      DO ni=1,nions

     ! distance to left image
         x1(:)=posion_all(:,ni,i-1)-posion_all(:,ni,i)
         x1(:)=MOD(x1(:)+100.5_q,1._q)-0.5_q
         CALL  dirkar(1,x1,a)
         norm1=sum(x1**2)
         d1=d1+norm1

     ! distance to right image
         x2(:)=posion_all(:,ni,i+1)-posion_all(:,ni,i)
         x2(:)=MOD(x2(:)+100.5_q,1._q)-0.5_q
         CALL  dirkar(1,x2,a)
         norm2=sum(x2**2)
         d2=d2+norm2
         
      ENDDO    

      d1=sqrt(d1)
      d2=sqrt(d2)

      ! local tangent for each ion

      DO ni=1,nions
         x(:)=(MOD(posion_all(:,ni,i+1)-posion_all(:,ni,i)+100.5_q,1._q)-0.5_q)*(d1**2) &
             -(MOD(posion_all(:,ni,i-1)-posion_all(:,ni,i)+100.5_q,1._q)-0.5_q)*(d2**2)
         x=x/((d1**2)*d2+d1*(d2**2))
         CALL  dirkar(1,x,a)
         norm=norm+sum(x**2)
         tangent(:,ni)=x
      ENDDO
         
      IF (norm/=0) tangent=tangent/sqrt(norm)

      e_chain=0
      e_image=0
      force_chain=0
      ftangent=0
!======================================================================
!
!  plain elastic band method
!
!======================================================================
      IF (spring/= 0) THEN
      const=ABS(spring)

      DO i=1,images+1
         DO ni=1,nions
            x(1)=posion_all(1,ni,i)-posion_all(1,ni,i-1)
            x(2)=posion_all(2,ni,i)-posion_all(2,ni,i-1)
            x(3)=posion_all(3,ni,i)-posion_all(3,ni,i-1)
            ! minimum image convention
            x(1)=MOD(x(1)+100.5_q,1._q)-0.5_q
            x(2)=MOD(x(2)+100.5_q,1._q)-0.5_q
            x(3)=MOD(x(3)+100.5_q,1._q)-0.5_q

            CALL  dirkar(1,x,a)
            d = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)

            e_chain=e_chain+ d*const/2

            ! add force to force_chain
            IF (i-1==node) THEN
               force_chain(1,ni)=force_chain(1,ni)+x(1)*const
               force_chain(2,ni)=force_chain(2,ni)+x(2)*const
               force_chain(3,ni)=force_chain(3,ni)+x(3)*const
               e_image=e_image+d*const/4
            ENDIF
            IF (i==node) THEN
               force_chain(1,ni)=force_chain(1,ni)-x(1)*const
               force_chain(2,ni)=force_chain(2,ni)-x(2)*const
               force_chain(3,ni)=force_chain(3,ni)-x(3)*const
               e_image=e_image+d*const/4
            ENDIF
         ENDDO
      ENDDO
      ENDIF
!======================================================================
!
!  nudged elastic band method
!   Hannes Jonsson and Greg Mills, preprint
!    (sub to Journal of Chem Phys.)
!  this is a little bit tricky, and the method might not work
!  with CG or molecular dynamics because no correction to energy
!  is calculated
!======================================================================
      IF ( spring <= 0 ) THEN
         ! force_chain contains already the spring force
         ! e_image must be reset to 0
         e_image=0
         e_chain=0
         ! project force and chain force onto local tangent
         proj =0
         proj2=0
         DO ni=1,nions
            proj =proj +tangent(1,ni)*force(1,ni) &
                       +tangent(2,ni)*force(2,ni) &
                       +tangent(3,ni)*force(3,ni)
            proj2=proj2+tangent(1,ni)*force_chain(1,ni) &
                       +tangent(2,ni)*force_chain(2,ni) &
                       +tangent(3,ni)*force_chain(3,ni)
         ENDDO
         ftangent   =(proj-proj2)
         ! force into direction of local tangent
         force_chain=tangent* (-proj+proj2)
      ENDIF
!======================================================================
!
! report important results to unit 6
!
!======================================================================
90    FORMAT( '  energy of chain is (eV) ',F16.6,' for this image ',F16.6,/ &
           '  tangential force (eV/A) ',F16.6, / &
           '  left and right image ',2F10.6,' A')

100   FORMAT( ' TANGENT     ',35X,'CHAIN-FORCE (eV/Angst)'/ &
           ' ----------------------------------------------', &
           '-------------------------------------')
110   FORMAT((3F13.5,3X,3F14.6))
120   FORMAT( ' ----------------------------------------------', &
           '-------------------------------------')

101   FORMAT( ' CHAIN + TOTAL  (eV/Angst)'/ &
           ' ----------------------------------------------')
111   FORMAT((3F13.5))
121   FORMAT( ' ----------------------------------------------')

      toten=toten+e_chain/images

      IF (iu6 >=0) THEN
         WRITE(iu6,90 ) e_chain,e_image,ftangent,d1,d2
         WRITE(iu6,100)
         WRITE(iu6,110) (tangent(:,ni),force_chain(:,ni),ni=1,nions)
         WRITE(iu6,120)
      ENDIF

      force=force+force_chain
      IF (iu6 >=0) THEN
         WRITE(iu6,101)
         WRITE(iu6,111) (force(:,ni),ni=1,nions)
         WRITE(iu6,121)
      ENDIF



#endif
    END SUBROUTINE chain_force

!**********************************************************************
!
! there are several points where a global sum between chains
! is required in order to get correct results
! terms like the total energy / or the kinetic energy
! make only sense in a global and not per image sense
!
!**********************************************************************

    SUBROUTINE sum_chain( value )
      REAL(q) :: value
      IF (images==0) RETURN
      IF (spring==-1000) RETURN

      CALLMPI_C( M_sum_d( comm_chain, value, 1 ))
    END SUBROUTINE sum_chain

!**********************************************************************
!
! also the logical break conditionions must be
!
!**********************************************************************

    SUBROUTINE and_chain( value )
      LOGICAL :: value
      REAL(q) :: sum

      IF (images==0) RETURN
      IF (spring==-1000) RETURN

! if one node is .FALSE., .FALSE. is returned on all nodes
      IF (value) THEN
         sum=0
      ELSE
         sum=1
      ENDIF
      CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))

      IF (sum>=1) THEN
         value=.FALSE.
      ELSE
         value=.TRUE.
      ENDIF
    END SUBROUTINE and_chain

!**********************************************************************
!
! initialize the chain (repeated image mode)
! read the spring constant
! and  the two outer images, these images are kept fixed
! during the entire simulation
!
!**********************************************************************

    SUBROUTINE chain_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

! needed only temporarily
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics)  :: DYN
      INTEGER     IDUM,IERR,N,idir,node
      CHARACTER*1   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!R.S
      integer tiu6, tiu0
        tiu6 = IO%IU6
        tiu0 = IO%IU0

! quick return, if we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
! read the spring constant
      spring=-5.

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPRING','=','#',';','F', &
     &            IDUM,spring,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''BMIX'' from file INCAR.'
         STOP
      ENDIF
      CALL XML_INCAR('SPRING','F',IDUM,spring,CDUM,LDUM,CHARAC,N)

! allocate the positions
      ALLOCATE(posion_all(3,t_info%nions,0:images+1))

! read 00/POSCAR file, a little bit of fiddling is required
      idir=0
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, &
           IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,0)= DYN%POSION
! read images+1/POSCAR file
      idir=images+1
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, &
           IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,images+1)= DYN%POSION

      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)
#endif
    END SUBROUTINE chain_init

!**********************************************************************
!
! returns true if hyper nudged elastic band method is used
!
!**********************************************************************

    FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
      IF (images==0 .OR. spring /= 0 ) THEN
         LHYPER_NUDGE=.FALSE.
      ELSE
         LHYPER_NUDGE=.TRUE.
      ENDIF

    END FUNCTION LHYPER_NUDGE
END MODULE chain
