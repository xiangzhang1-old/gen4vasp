!
! poscar description input file
! only included if MODULES are not supported
!
      TYPE type_info
!only T_INFO
        CHARACTER*40 SZNAM2           ! name of poscar file
        INTEGER NTYPD                 ! dimension for types
        INTEGER NTYP                  ! number of types
        INTEGER NTYPPD                ! dimension for types inc. empty spheres
        INTEGER NTYPP                 ! number of types empty spheres
        INTEGER NIOND                 ! dimension for ions
        INTEGER NIONPD                ! dimension for ions inc. empty spheres
        INTEGER NIONS                 ! actual number of ions
        INTEGER NIONP                 ! actual number of ions inc. empty spheres
        LOGICAL LSDYN                 ! selective dynamics (yes/ no)
        LOGICAL LDIRCO                ! positions in direct/recproc. lattice
        REAL(q), POINTER :: POSION(:,:)  ! positions usually same as DYN%POSION
        LOGICAL,POINTER ::  LSFOR(:,:) ! selective dynamics
        INTEGER, POINTER :: ITYP(:)   ! type for each ion
        INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
        REAL(q), POINTER :: POMASS(:) ! mass for each ion type
        REAL(q), POINTER :: RWIGS(:)  ! wigner seitz radius for each ion type
        REAL(q), POINTER :: ROPT(:)   ! optimization radius for each type
        REAL(q), POINTER :: ATOMOM(:) ! initial local spin density for each ion
      END TYPE


      TYPE dynamics
!only DYN
        REAL(q), POINTER :: POSION(:,:) ! positions
        REAL(q), POINTER :: POSIOC(:,:) ! old positions
        REAL(q), POINTER :: VEL(:,:)  ! velocities
        REAL(q), POINTER :: D2(:,:)   ! predictor corrector/coordinates
        REAL(q), POINTER :: D2C(:,:)  ! predictor corrector/coordinates
        REAL(q), POINTER :: D3(:,:)   ! predictor corrector/coordinates
        REAL(q) A(3,3)                ! current lattice (presently unused)
        REAL(q) AC(3,3)               ! old lattice (presently unused)
        REAL(q) SNOSE(4)              ! nose thermostat
        INTEGER IBRION                ! mode for relaxation
        INTEGER ISIF                  ! mode for stress/ ionic relaxation
        REAL(q) POTIM                 ! time step
        REAL(q) EDIFFG                ! accuracy for ionic relaxation
        REAL(q), POINTER :: POMASS(:) ! mass of each ion for dynamics
        REAL(q) SMASS                 ! mass of nose thermostat
        REAL(q) PSTRESS               ! external pressure
        REAL(q) TEBEG, TEEND          ! temperature during run
        REAL(q) TEMP                  ! current temperature
        INTEGER NSW                   ! number of ionic steps
        INTEGER NBLOCK,KBLOCK         ! blocks
        INTEGER INIT                  ! predictore corrector initialized
	INTEGER NFREE                 ! estimated ionic degrees of freedom
      END TYPE





