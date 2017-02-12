#include "symbol.inc"
!********** TUTORIAL PACKAGE -- GIVE WARNINGS AND ADVICES **************
! RCS:  $Id: tutor.F,v 1.6 2003/06/27 13:22:23 kresse Exp kresse $
!                                                                      *
      SUBROUTINE VTUTOR(TYPE,WTOPIC,RDAT,NR,IDAT,NI, &
     &                  CDAT,NC,LDAT,NL,IU,IDIOT)
      USE prec
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      IMPLICIT COMPLEX(q) (C)
!                                                                      *
!***********************************************************************
!                                                                      *
!  This routine gives warnings, advices, error messages on very very   *
!  important things which are often done wrong by bloody newbies ...   *
!                                                                      *
!  Variables:                                                          *
!                                                                      *
!    TYPE   tells us what it is (fatal E rror, W arning, A dvice       *
!           or fatal error which makes it necessary to S top ...       *
!    TOPIC  contains some identifier string about what to talk ...     *
!    RDAT,IDAT,CDAT and LDAT are arrays containing possible data       *
!    being necessary for the one or the other message (type Real,      *
!    Integer, Complex, Logical ...) with NR,NI,NC and NL being the     *
!    dimensions of these arrays ... .                                  *
!    IDIOT  flags the 'expert level' of the user:                      *
!           0: 'complete expert' (no messages at all)                  *
!           1: 'almost complete expert' (only hard errors)             *
!           2: 'somehow experienced user' (only warnings and errors)   *
!           3: 'complete idiot' (all kind of messages ...)             *
!                                                                      *
!  On output you should receive messages on I/O-unit IU telling you    *
!  what to do or better not to do (i.e. what you have done wrong ...)  *
!                                                                      *
!***********************************************************************

      CHARACTER*1   TYPE
      CHARACTER*255 TOPIC
      CHARACTER*(*) WTOPIC
      LOGICAL       LDAT, LIO
      DIMENSION RDAT(NR),IDAT(NI),CDAT(NC),LDAT(NL)

      LIO=.TRUE.
      IF (IU<0) LIO=.FALSE.

      TOPIC=WTOPIC
      CALL STRIP(TOPIC,LTOPIC,'B')

! Header of the message (if there is a message at all ...)
      IF ((TYPE=='E').OR.(TYPE=='S')) THEN
! the 'complete expert' needs no error messages, warnings or advices ...
         IF (IDIOT<=0) RETURN
         IF (IU>=0) WRITE(IU,1)
      ELSE IF ((TYPE=='W').OR.(TYPE=='U')) THEN
! the 'complete expert' needs no error messages, warnings or advices ...
! the 'almost complete expert' needs no warnings or advices ...
         IF ((TYPE=='W').AND.(IDIOT<=1)) RETURN
! ... but 'U'rgent warning (almost like error!) shall be given to all
! except for 'the complete expert' needing no help at all  :-)
         IF ((TYPE=='U').AND.(IDIOT<=0)) RETURN
         IF (LIO) WRITE(IU,2)
      ELSE
! the 'complete expert' needs no error messages, warnings or advices ...
! the 'almost complete expert' needs no warnings or advices ...
! the 'somehow experienced user' needs no advices ...
         IF (IDIOT<=2) RETURN
         IF (LIO) WRITE(IU,3)
      ENDIF
! ... but the 'complete idiot' needs all ...                   :-)

    1 FORMAT(/' ------------------------------------------------', &
     & '----------------------------- '/, &
     & '|                                                ', &
     & '                             |'/, &
     & '|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR ', &
     & '     ###     ###     ###     |'/, &
     & '|     E        R     R  R     R  O     O  R     R', &
     & '     ###     ###     ###     |'/, &
     & '|     E        R     R  R     R  O     O  R     R', &
     & '     ###     ###     ###     |'/, &
     & '|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR ', &
     & '      #       #       #      |'/, &
     & '|     E        R   R    R   R    O     O  R   R  ', &
     & '                             |'/, &
     & '|     E        R    R   R    R   O     O  R    R ', &
     & '     ###     ###     ###     |'/, &
     & '|     EEEEEEE  R     R  R     R  OOOOOOO  R     R', &
     & '     ###     ###     ###     |'/, &
     & '|                                                ', &
     & '                             |')


    2 FORMAT(/ &
     & ' -----------------------------------------------', &
     & '------------------------------ '/, &
     & '|                                               ', &
     & '                              |'/, &
     & '|           W    W    AA    RRRRR   N    N  II  ', &
     & 'N    N   GGGG   !!!           |'/, &
     & '|           W    W   A  A   R    R  NN   N  II  ', &
     & 'NN   N  G    G  !!!           |'/, &
     & '|           W    W  A    A  R    R  N N  N  II  ', &
     & 'N N  N  G       !!!           |'/, &
     & '|           W WW W  AAAAAA  RRRRR   N  N N  II  ', &
     & 'N  N N  G  GGG   !            |'/, &
     & '|           WW  WW  A    A  R   R   N   NN  II  ', &
     & 'N   NN  G    G                |'/, &
     & '|           W    W  A    A  R    R  N    N  II  ', &
     & 'N    N   GGGG   !!!           |'/, &
     & '|                                               ', &
     & '                              |')

    3 FORMAT(/ &
     & ' -----------------------------------------------', &
     & '------------------------------ '/, &
     & '|                                               ', &
     & '                              |'/, &
     & '|  ADVICE TO THIS USER RUNNING ''VASP/VAMP''  ', &
     & ' (HEAR YOUR MASTER''S VOICE ...):  |'/, &
     & '|                                               ', &
     & '                              |')

    4 FORMAT( &
     & '|                                               ', &
     & '                              |'/, &
     & ' -----------------------------------------------', &
     & '------------------------------ '/)

    5 FORMAT( &
     & '|                                               ', &
     & '                              |'/, &
     & '|      ---->  I REFUSE TO CONTINUE WITH THIS ', &
     & 'SICK JOB ..., BYE!!! <----       |'/, &
     & '|                                               ', &
     & '                              |'/, &
     & ' -----------------------------------------------', &
     & '------------------------------ '/)

! Now following the long long long long code printing all messages ...:
! =====================================================================

      IF (LIO.AND. TOPIC(1:LTOPIC)=='REAL-SPACE WITHOUT OPTIMIZATION') THEN
         WRITE(IU,'(A)')'|     The real-space-projection scheme for '// &
     &                  'the treatment of the nonlocal      |'
         WRITE(IU,'(A)')'|     pseudopotentials has been switched on'// &
     &                  ' -- but on file POTCAR I have      |'
         WRITE(IU,'(A)')'|     not found any entries signalling me '// &
     &                  'that you have ever performed a      |'
         WRITE(IU,'(A)')'|     real-space optimization of all '// &
     &                  'nonlocal projectors! BE WARNED that      |'
         WRITE(IU,'(A)')'|     a calculation using the real-space-'// &
     &                  'projection scheme together with      |'
         WRITE(IU,'(A)')'|     nonlocal projectors which were not '// &
     &                  'real-space-optimized (according      |'
         WRITE(IU,'(A)')'|     to the plane-wave-cutoff used in this '// &
     &                  'calculation) might give more      |'
         WRITE(IU,'(A)')'|     or less inaccurate results         '// &
     &                  '                                     |'
         WRITE(IU,'(A)')'|     I hope you know what you are doing   '// &
     &                  '                                   |'
         GOTO 99999
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='VASP.4.4') THEN
         WRITE(IU,'(A)')'|     You are running vasp.4.5 in the vasp.4.4 compatibility mode             |'
         WRITE(IU,'(A)')'|     vasp.4.5 has some numerical improvements, which are not applied in this |'
         WRITE(IU,'(A)')'|     mode (charge at unbalanced lattice vectors are no longer zeroed,        |'
         WRITE(IU,'(A)')'|           PAW augmentation charges are integrated more accurately)          |'
        GOTO 99999
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='NO REAL-SPACE AND YOU COULD') THEN
         WRITE(IU,'(A)')'|      You have a (more or less) ''large '// &
     &              'supercell'' and for larger cells       |'
         WRITE(IU,'(A)')'|      it might be more efficient to '// &
     &           'use real space projection operators      |'
         WRITE(IU,'(A)')'|      So try LREAL= Auto  in the INCAR   '// &
     &                 'file.                               |'
         WRITE(IU,'(A)')'|      Mind: If you want to do an extremely'// &
     &                 ' accurate calculations keep the    |'
         WRITE(IU,'(A)')'|      reciprocal projection scheme         ' &
     &                 //' (i.e. LREAL=.FALSE.)             |'
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='NO REAL-SPACE AND YOU SHOULD') THEN
         WRITE(IU,'(A)')'|      You have a (more or less) ''large '// &
     &              'supercell'' and for larger cells       |'
         WRITE(IU,'(A)')'|      it might be more efficient to '// &
     &           'use real space projection opertators     |'
         WRITE(IU,'(A)')'|      So try LREAL= Auto  in the INCAR   '// &
     &                 'file.                               |'
         WRITE(IU,'(A)')'|      Mind: At the moment your POTCAR file'// &
     &                 ' does not contain real space       |'
         WRITE(IU,'(A)')'|       projectors, and has to be modified,'// &
     &                 '  BUT if you                       |'
         WRITE(IU,'(A)')'|      want to do an extremely '// &
     &                ' accurate calculation you might also keep the  |'
         WRITE(IU,'(A)')'|      reciprocal projection scheme         ' &
     &                 //' (i.e. LREAL=.FALSE.)             |'

      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='IALGO8') THEN
         WRITE(IU,'(A)')'|      Recently Corning got a patent for the'// &
     &                  ' Teter Allan Payne algorithm      |'
         WRITE(IU,'(A)')'|      therefore VASP.4.5 does not support IALG'// &
     &                  'O=8 any longer                 |'
         WRITE(IU,'(A)')'|      a much faster algorithm, IALGO=38, is '// &
     &                  ' now implemented in VASP         |'
         WRITE(IU,'(A)')'|      this algorithm is a blocked Davidson'// &
     &                  ' like method and as reliable as    |'
         WRITE(IU,'(A)')'|      IALGO=8 used to be                  '// &
     &                  '                                   |'
         WRITE(IU,'(A)')'|      for ultimate performance IALGO=48 is'// &
     &                  ' still the method of choice        |' 
         WRITE(IU,'(A)')'|      -- SO MUCH ABOUT PATENTS :)         '// &
     &                  '                                   |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='REAL-SPACE NOMORE RECOMMENDED') THEN
         WRITE(IU,'(A)')'|      You have a (more or less) ''small '// &
     &              'supercell'' and for smaller cells      |'
         WRITE(IU,'(A)')'|      it is recommended  to use the '// &
     &            'reciprocal-space projection scheme!      |'
         WRITE(IU,'(A)')'|      The real space optimization is not  '// &
     &                  'efficient for small cells and it   |'
         WRITE(IU,'(A)')'|      is also less accurate ...            ' &
     &                 //'                                  |'
         WRITE(IU,'(A)')'|      Therefore set LREAL=.FALSE. in the  '// &
     &                  'INCAR file                         |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='WRONG OPTIMZATION REAL-SPACE') THEN
         WRITE(IU,'(A)') '|      One real space projector is optimized' &
     &                 //' for                              |'
         WRITE(IU,'(A,F10.2,A)') &
     &              '|      E    =',RDAT(1),', eV  but you are using a ' &
     &               //' cutoff of                   |'
         WRITE(IU,'(A,F10.2,A,F10.2,A)')  &
     &              '|      ENMAX=',RDAT(2),' eV  ( QCUT=',RDAT(3), &
     &              ' a.u.)                           |'
         WRITE(IU,'(A)') '|      This makes no sense reoptimize the  ' &
     &                 //' projector                         |'
         WRITE(IU,'(A)') '|      with the a.u. value given above      ' &
     &                 //'                                  |'
      ENDIF


      IF (LIO.AND. TOPIC(1:LTOPIC)=='DIFFERENT XCGRAD TYPES') THEN
         WRITE(IU,'(A)') '|      You have build up your multi-ion-type' &
     &                 //' POTCAR file out of POTCAR        |'
         WRITE(IU,'(A)') '|      files with incompatible specifications' &
     &                  //' for the XC-types used to        |'
         WRITE(IU,'(A)') '|      generate the pseudopotential. This '// &
     &                 'makes no sense at all!! What        |'
         WRITE(IU,'(A,I2,A,I3,A,I2,A)') '|      I found is LEXCHG = ', &
     &             IDAT(1),'  for atom types <= ',IDAT(3)-1, &
     &                    ' but LEXCHG = ',IDAT(2),'        |'
         WRITE(IU,'(A,I3,A)') '|      was found for atom type = ', &
     &            IDAT(3),'. Use identical XC-functionals for        |'
         WRITE(IU,'(A)') '|      the pseudopotential generation for '// &
     &                 'all atom types, please ... !        |'
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='DIFFERENT LDA-XC TYPES') THEN
         WRITE(IU,'(A)') '|      You have build up your multi-ion-type' &
     &                 //' POTCAR file out of POTCAR        |'
         WRITE(IU,'(A)') '|      files with incompatible specifications' &
     &                  //' for the XC-types used to        |'
         WRITE(IU,'(A)') '|      generate the pseudopotential. This '// &
     &                 'makes no sense at all!! What        |'
         WRITE(IU,'(A,I2,A,I3,A,I2,A)') '|      I found is LEXCH  = ', &
     &             IDAT(1),'  for atom types <= ',IDAT(3)-1, &
     &                    ' but LEXCH  = ',IDAT(2),'        |'
         WRITE(IU,'(A,I3,A)') '|      was found for atom type = ', &
     &            IDAT(3),'. Use identical XC-functionals for        |'
         WRITE(IU,'(A)') '|      the pseudopotential generation for '// &
     &                 'all atom types, please ... !        |'
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='DIFFERENT REL-XC TYPES') THEN
         WRITE(IU,'(A)') '|      You have build up your multi-ion-type' &
     &                 //' POTCAR file out of POTCAR        |'
         WRITE(IU,'(A)') '|      files with incompatible specifications' &
     &                  //' for the XC-types used to        |'
         WRITE(IU,'(A)') '|      generate the pseudopotential. This '// &
     &                 'makes no sense at all!! What        |'
         WRITE(IU,'(A)') '|      I found is that the flag which '// &
     &             'switches on/off the relativistic        |'
         WRITE(IU,'(A,L1,A,I3,A)') '|      corrections has been set .', &
     &                        LDAT(1),'. for atom types <= ', &
     &                        IDAT(1)-1,' but it was        |'
         WRITE(IU,'(A,L1,A,I3,A)') '|      set .',LDAT(2), &
     &           '. for atom type no. ',IDAT(1), &
     &           '. Use identical XC-functionals for        |'
         WRITE(IU,'(A)') '|      the pseudopotential generation for '// &
     &                 'all atom types, please ... !        |'
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='DIFFERENT XC FILES') THEN
         WRITE(IU,'(A)') '|      The XC type on the POTCAR file is not' &
     &                 //' compatibel with the XC type on   |'
         WRITE(IU,'(A)') '|      the EXHCAR file                       ' &
     &                  //'                                 |'
         GOTO 99999
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='NOLDAU') THEN
         WRITE(IU,'(A)') '|      VASP supports LDA+U only for PAW pote' &
     &                 //'ntials and not for US or NC       |'
         WRITE(IU,'(A)') '|      pseudopotentials                      ' &
     &                  //'                                 |'
         WRITE(IU,'(A)') '|      please restart with with the approriat' &
     &                  //'e potentials                     |'
         GOTO 99999
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='NOWANNIER') THEN
         WRITE(IU,'(A)') '|      VASP supports the calculation of maxi' &
     &                 //'mally localized wannier functions |'
         WRITE(IU,'(A)') '|      only in the Gamma-only version of the' &
     &                 //' code.                            |'
         GOTO 99999
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='LREALA') THEN
         WRITE(IU,'(A)') '|      LREAL=A is not well tested yet, pleas' &
     &                 //'e use LREAL=O                     |'
         WRITE(IU,'(A)') '|      or test your results very carefully   ' &
     &                  //'                                 |'
         GOTO 99999
      ENDIF

      IF (LIO.AND. TOPIC(1:LTOPIC)=='DIFFERENT SLATER-XC TYPES') THEN
         WRITE(IU,'(A)') '|      You have build up your multi-ion-type' &
     &                 //' POTCAR file out of POTCAR        |'
         WRITE(IU,'(A)') '|      files with incompatible specifications' &
     &                  //' for the XC-types used to        |'
         WRITE(IU,'(A)') '|      generate the pseudopotential. This '// &
     &                 'makes no sense at all!! What        |'
         WRITE(IU,'(A,F9.6,A,I3,A)') '|      I found is slater '// &
     &                       'parameter = ',RDAT(1),' for atom '// &
     &                       'types <= ',IDAT(1)-1,'        |'
         WRITE(IU,'(A,F9.6,A,I3,A)') '|      but slater parameter = ', &
     &               RDAT(2),' was found for atom type = ', &
     &                                    IDAT(1),'.        |'
         WRITE(IU,'(A)') '|      Use identical slater-parameters '// &
     &              'in the exchange functionals for        |'
         WRITE(IU,'(A)') '|      the pseudopotential generation for '// &
     &                 'all atom types, please ... !        |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='partial DOS') THEN
         WRITE(IU,'(A)') '|      The partial DOS and the PROCAR file are' &
     &                 //' not evaluated for NPAR/=1      |'
         WRITE(IU,'(A)') '|      please rerun with NPAR=1             '&
     &                  //'                                  |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='nooptics') THEN
         WRITE(IU,'(A)') '|      The optical properties can be evaluated' &
     &                 //' only for NPAR=1                |'
         WRITE(IU,'(A)') '|      please rerun with NPAR=1             '&
     &                  //'                                  |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='POSITION') THEN
         WRITE(IU,'(A)') '|      The distance between some ions is very'&
     &                 //' small                           |'
         WRITE(IU,'(A)') '|      please check the nearest neigbor list'&
     &                  //' in the OUTCAR file               |'
         WRITE(IU,'(A)') '|          I HOPE YOU KNOW, WHAT YOU ARE ' &
     &                  //' DOING                               |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='DEGREES OF FREEDOM') THEN
         WRITE(IU,'(A,I6,A)') '|      VASP found ',IDAT(1),' degrees of freedom  '&
     &                 //'                                 |'
         WRITE(IU,'(A)') '|      the temperature will equal 2*E(kin)/ '&
     &                  //'(degrees of freedom)              |'
         WRITE(IU,'(A)') '|      this differs from previous rel' &
     &                  //'eases, where T was 2*E(kin)/(3 NIONS).   |'
         WRITE(IU,'(A)') '|      The new definition is more con' &
     &                  //'sistent                                  |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='CENTER OF MASS DRIFT') THEN
         WRITE(IU,'(A)') '|      The initial velocities result in a cen'&
     &                 //'ter of mass drift but            |'
         WRITE(IU,'(A)') '|      there must be no drift.              '&
     &                  //'                                  |'
         WRITE(IU,'(A)') '|      The drift will be removed !          '&
     &                  //'                                  |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='ENFORCED LDA') THEN
         WRITE(IU,'(A)') '|      You enforced a specific xc-type in the' &
     &                 //' INCAR file,                     |'
         WRITE(IU,'(A)') '|      a different type was found on the ' &
     &                  //'POTCAR file                          |'
         WRITE(IU,'(A)') '|          I HOPE YOU KNOW, WHAT YOU ARE ' &
     &                  //' DOING                               |'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='LINTET') THEN
         WRITE(IU,'(A)') '|      The linear tetrahedron method can not ' &
     &                 //' be used with the KPOINTS file   |'
         WRITE(IU,'(A)') '|      (generation of strings of k-points' &
     &                  //')                                    |'
 
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='GAMMAK') THEN
         WRITE(IU,'(A)') '|      You are using the Gamma-point only ver' &
     &                 //'sion with more than one k-point  |'
         WRITE(IU,'(A)') '|      or some other non Gamma k-point   ' &
     &                  //')                                    |'
 
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='KPOINTS') THEN
         WRITE(IU,'(A)') '|      Error reading KPOINTS file            ' &
     &                 //'                                 |'
         WRITE(IU,'(A,I5,A)') '|      the error occured at line:       ',IDAT(1), &
     &                 '                                 |'
 
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='KPOINTS') THEN
         WRITE(IU,'(A)') '|      Error reading KPOINTS file            ' &
     &                 //'                                 |'
         WRITE(IU,'(A,I5,A)') '|      the error occured at line:       ',IDAT(1), &
     &                 '                                 |'
 
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='LNOINVERSION') THEN
         WRITE(IU,'(A)') '|      Full k-point grid generated           ' &
     &                 //'                                 |'
         WRITE(IU,'(A)') '|      Inversion symmetry is not applied     ' &
     &                 //'                                 |'
 
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='FFT-GRID IS NOT SUFFICIENT') THEN
         WRITE(IU,'(A)') '|      Your FFT grids (NGX,NGY,NGZ) are not ' &
     &                 //'sufficient for an accurate        |'
         WRITE(IU,'(A)') '|      calculation.                         ' &
     &                 //'                                  |'
         WRITE(IU,'(A)') '|      The results might be wrong            ' &
     &                  //'                                 |'
         WRITE(IU,'(A)') '|      good settings for NGX NGY and ' &
     &               //' NGZ are                                 |'
         WRITE(IU,'(A,2I4,A,I4,A)') '|                      ', &
     &             IDAT(1),IDAT(2),'  and',IDAT(3), &
     &             '                                      |'
         WRITE(IU,'(A)') '|     Mind: This setting results in a small' &
     &                 //' but reasonable wrap around error  |'
         WRITE(IU,'(A)') '|     It is also necessary to adjust these ' &
     &                 //' values to the FFT routines you use|'
      ENDIF
      IF (LIO.AND. TOPIC(1:LTOPIC)=='METAGGA and forces') THEN
         WRITE(IU,'(A)') '|      You have switched METAGGA and/or ASPH' &
     &                 //'ERICAL radial calculations on     |'
         WRITE(IU,'(A)') '|      but ions are moved based on forces.  ' &
     &                 //'                                  |'
         WRITE(IU,'(A)') '|      METAGGA and ASPH change only the exch' &
     &                 //'ange and correlation energy.      |'
         WRITE(IU,'(A)') '|      You might reconsider your choices.   ' &
     &                 //'                                  |'
         WRITE(IU,'(A,L4,A,L4,A)') &
                         '|      Current values: LASPH=',LDAT(1), &
     &                   ', and LMETAGGA=',LDAT(2), &
     &             '                          |'
      ENDIF
      
      IF (LIO.AND. TOPIC(1:LTOPIC)=='HIGHEST BANDS OCCUPIED') THEN
         WRITE(IU,'(A)') '|      Your highest band is occupied at some k-points! Unless you are         |'
         WRITE(IU,'(A)') '|      performing a calculation for an insulator or semiconductor, without    |'
         WRITE(IU,'(A)') '|      unoccupied bands, you have included TOO FEW BANDS!! Please increase    |'
         WRITE(IU,'(A)') "|      the parameter NBANDS in file 'INCAR' to ensure that the highest band   |"
         WRITE(IU,'(A)') '|      is unoccupied at all k-points. It is always recommended that one       |'
         WRITE(IU,'(A)') '|      include a few unoccupied bands to accellerate the convergence of       |'
         WRITE(IU,'(A)') '|      molecular dynamics runs (even for insulators or semiconductors)        |'
         WRITE(IU,'(A)') '|      Because the presence of unoccupied bands improves wavefunction         |'
         WRITE(IU,'(A)') "|      prediction, and helps to suppress 'band-crossings.'                    |"

         IF (NI==NR) THEN
            WRITE(IU,'(A)') '|      Following all k-points will be '// &
     &                'listed (with the Fermi weights of       |'
            WRITE(IU,'(A)') '|      the highest band given in '// &
     &           'paranthesis) ... :                           |'
            WRITE(IU,'(A)') '|                                     '// &
     &                '                                        |'
            DO 100 I=1,NI-1
               WRITE(IU,'(A,I5,A,F8.5,A)') &
     &              '|                      ',IDAT(I),'       (', &
     &              RDAT(I),')                                 |'
  100       CONTINUE
            WRITE(IU,'(A)') '|                                     '// &
     &                '                                        |'
            WRITE(IU,'(A,I5,A,F8.5,A)') &
     &              '|      The total occupancy of band no. ',IDAT(NI), &
     &                  ' is  ',RDAT(NI),' electrons ...       |'
         ENDIF
      ENDIF

! Here we have told the user all we can tell, bring it to an end ...:
99999 CONTINUE
! Hmmmm ..., very very fatal error --->  S T O P  !!
      IF (TYPE=='S') THEN
! Well announce the brute end and exit ...
         IF (LIO) WRITE(IU,5)
         STOP
      ENDIF
! Normal end of message ...
         IF (LIO) WRITE(IU,4)
! Bye ...
      RETURN
      END
