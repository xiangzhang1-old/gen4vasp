!  Some important Parameters, to convert to a.u.
!  - AUTOA  = 1. a.u. in Angstroem
!  - RYTOEV = 1 Ry in Ev
!  - EVTOJ  = 1 eV in Joule
!  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
!  - BOLKEV = Boltzmanns constant in eV/K
!  - BOLK   = Boltzmanns constant in Joule/K

      REAL(q), PARAMETER :: AUTOA=0.529177249_q,RYTOEV=13.605826_q
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q, &
     &           BOLKEV=8.6173857E-5_q,BOLK=BOLKEV*EVTOJ

      REAL(q), PARAMETER :: EVTOKCAL=23.06
! FELECT = (the electronic charge)/(4*pi*the permittivity of free space)
!         in atomic units that is just e^2
! EDEPS = electron charge divided by the permittivity of free space
!         in atomic units that is just 4 pi e^2
! HSQDTM = (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
!
      REAL(q),PARAMETER  :: PI =3.141592653589793238_q,TPI=2*PI
      COMPLEX(q),PARAMETER  :: CITPI = (0._q,1._q)*TPI
      REAL(q),PARAMETER  :: FELECT = 2*AUTOA*RYTOEV,EDEPS=4*PI*2*RYTOEV*AUTOA,&
     &                   HSQDTM = RYTOEV*AUTOA*AUTOA

      REAL(q),PARAMETER :: AUTOA2=AUTOA *AUTOA
      REAL(q),PARAMETER :: AUTOA3=AUTOA2*AUTOA
      REAL(q),PARAMETER :: AUTOA4=AUTOA2*AUTOA2
      REAL(q),PARAMETER :: AUTOA5=AUTOA3*AUTOA2
