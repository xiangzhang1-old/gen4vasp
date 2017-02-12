!
!  structure required to support storage of wavefunction for
!  several kpoints and spins
!
      TYPE wavedes
!only WDES
         REAL(q)    RSPIN              ! spin multiplicity
         REAL(q)    ENMAX              ! energy cutoff
         INTEGER NRSPINORS             ! number of spinors (1 for collinear, 2 for non collinear)
         INTEGER NGDIM                 ! first dimension of any array related to the basis
         INTEGER NRPLWV                ! first dimension of wavefunction array
                                       ! collinear:  NRPLWV=NGDIM, noncollinear:  NRPLWV=2*NGDIM
         INTEGER NPROD                 ! first dimension of projected wave array
         INTEGER NPRO                  ! local number of elements in projected wave array
         INTEGER NPRO_TOT              ! total number of elements
                                       ! NPRO, NPROD, and NPRO_TOT are all doubled in the non collinear version
         INTEGER NBANDS                ! local number of bands
         INTEGER NB_TOT                ! total number bands
         INTEGER NB_PAR                ! number of bands treated simultaneosly
         INTEGER NSIM                  ! band blocking (mainly for seriel version)
         INTEGER NB_LOW                ! lowest band index in global
         INTEGER NKDIM                 ! dim for number of k-points
         INTEGER NKPTS                 ! number of k-points
         INTEGER ISPIN                 ! number of spins
	 INTEGER NCDIJ                 ! dimension of arrays like CDIJ, CQIJ
         INTEGER NIONS                 ! number of ions required for the layout of CPROJ
         INTEGER NTYP                  ! number of types

         INTEGER,POINTER :: NPLWKP(:)  ! number of coefficients for each k-point and band per node
	 INTEGER,POINTER :: NGVECTOR(:)! number of G-vectors in the basis for each k-point per node
                                       ! collinear: NPLWKP= NGVECTOR, noncollinear NPLWKP = 2*NGVECTOR 
                                       ! NGVECTOR is the same for collinear and non collinear calculations
         INTEGER,POINTER :: NPLWKP_TOT(:) ! total number of coefficients in plane wave array 
                                       ! (summed over nodes, doubled in the non collinear case)
      ! few things which are only required in parallel version
         INTEGER         :: NCOL       ! number of columns
         INTEGER,POINTER :: PL_INDEX(:,:) ! index a column would have in serial version
         INTEGER,POINTER :: PL_COL(:,:)! number of plane wave in this column
         INTEGER,POINTER ::NPRO_POS(:) ! index CPROJ would have in serial version

         INTEGER,POINTER :: LMMAX(:)   ! required for the layout of CPROJ
         INTEGER,POINTER :: NITYP(:)   ! required for the layout of CPROJ
         REAL(q),POINTER :: VKPT(:,:)  ! coordinate of k-point
         REAL(q),POINTER :: WTKPT(:)   ! symmetry weight-factor for each k-point
         INTEGER,POINTER :: NINDPW(:,:)! index to the FFT box for each pw comp and k-point
         INTEGER,POINTER :: IGX(:,:)   ! x index of each pw comp and k-point
         INTEGER,POINTER :: IGY(:,:)   ! y index of each pw comp and k-point
         INTEGER,POINTER :: IGZ(:,:)   ! z index of each pw comp and k-point
!-MM- changes to accommodate spin spirals
! original statement
!        REAL(q),POINTER :: DATAKE(:,:) ! kinetic energy for each plane wave
         REAL(q),POINTER :: DATAKE(:,:,:) ! kinetic energy for each plane wave
                                          ! last index labels up and down components
                                          ! of the spinor in case of spin spirals
         LOGICAL LSPIRAL               ! calculate spin spirals?
         LOGICAL LZEROZ                ! set m_z to zero in SET_CHARGE?
         REAL(q) QSPIRAL(3)            ! propagation vector of spin spiral
!-MM- end of alteration+addition
         TYPE(communic),POINTER  :: COMM,COMM_INTER,COMM_INB
         LOGICAL LCRITICAL_MEM         ! save some memory by accepting slow down
         LOGICAL LNONCOLLINEAR         ! allows one to turn on noncollinear calculations
	 LOGICAL LSORBIT               ! spin orbit coupling
	 REAL(q) SAXIS(3)              ! quantisation axis of the spin operator
      END TYPE

!
! description for one k point
! contains also all information required for simple calculations in real space
!
      TYPE wavedes1
!only WDES1
         REAL(q)    RSPIN              ! spin multiplicity
         INTEGER NRSPINORS             ! number of spinors (1 for collinear, 2 for non collinear)
         INTEGER NGDIM                 ! first dimension of any array related to the basis
         INTEGER NRPLWV                ! first dimension of wavefunction array
                                       ! collinear:  NRPLWV=NGDIM, noncollinear:  NRPLWV=2*NGDIM
         INTEGER NPROD                 ! first dimension of projected wave array
         INTEGER NPRO                  ! local number of elements in projected wave array
         INTEGER NPRO_TOT              ! total number of elements
                                       ! NPRO, NPROD, and NPRO_TOT are all doubled in the non collinear version
         INTEGER NBANDS                ! bands
         INTEGER NB_TOT                ! total number bands
         INTEGER NB_PAR                ! number of bands treated simultaneosly
         INTEGER NSIM                  ! band blocking (mainly for serial version)
         INTEGER NB_LOW                ! lowest band index in global 
         INTEGER NPL                   ! number of plane waves coefficients (local)
         INTEGER NGVECTOR              ! number of G-vectors in the basis (local)
                                       ! collinear: NGVECTOR == NPL, noncollinear 2*NGVECTOR == NPL
         INTEGER NPL_TOT               ! total number of plane waves (global)
!        INTEGER ISPIN                 ! number of spins	      
! 	 INTEGER NCDIJ                 ! dimension of arrays like CDIJ, CQIJ
         INTEGER NIONS                 ! number of ions required for the layout of CPROJ
         INTEGER NTYP                  ! number of types
         INTEGER NPLWV                 ! total number of plane waves in real space
         INTEGER NPLWVL                ! local number of plane waves in real space
         INTEGER MPLWV                 ! dimension of arrays in real space
         REAL(q) RINPL                 ! inverse of total number of plane waves
         INTEGER NK                    ! k-point number (required for non-local proj.)
      ! few things which are only required in parallel version
         INTEGER         :: NCOL       ! number of columns
         INTEGER,POINTER :: PL_INDEX(:)! index a column would have in serial version
         INTEGER,POINTER :: PL_COL(:)  ! number of plane waves in this column
         INTEGER,POINTER :: NPRO_POS(:)! index CPROJ would have in serial version
         INTEGER,POINTER :: LMMAX(:)   ! required for the layout of CPROJ
         INTEGER,POINTER :: NITYP(:)   ! required for the layout of CPROJ
         INTEGER,POINTER :: NINDPW(:)  ! index to the FFT box for each pw comp
         INTEGER,POINTER :: IGX(:)     ! x index of each pw comp and k-point
         INTEGER,POINTER :: IGY(:)     ! y index of each pw comp and k-point
         INTEGER,POINTER :: IGZ(:)     ! z index of each pw comp and k-point
!-MM- changes to accommodate spin spirals
! original statement
!        REAL(q),POINTER :: DATAKE(:)  ! kinetic energy for each plane wave
         REAL(q),POINTER :: DATAKE(:,:) ! kinetic energy for each plane wave
                                        ! last index labels up and down components
                                        ! of the spinor in case of spin spirals
         LOGICAL LSPIRAL               ! calculate spin spirals?
         LOGICAL LZEROZ                ! set m_z to zero in SET_CHARGE?
         REAL(q) QSPIRAL(3)            ! propagation vector of spin spiral
!-MM- end of alteration+addition
         TYPE(communic),POINTER  :: COMM,COMM_INTER,COMM_INB
         LOGICAL LCRITICAL_MEM         ! save some memory by accepting slow down
         LOGICAL LNONCOLLINEAR         ! allows one to turn on noncollinear calculations
	 LOGICAL LSORBIT               ! spin orbit coupling
	 REAL(q) SAXIS(3)              ! quantisation axis of the spin operator
         END TYPE

!
! structure required to store a set of wavefunctions (bands & kpoints)
! we use two pointers for the projected wavefunctions
! one for complex and one for the real version
! symbol.inc (CPP) switches for real version from CPROJ -> GPROJ
! defined: WUP WDW
      TYPE wavefun
!only WUP
         TYPE(wavedes),POINTER:: WDES  ! descriptor for wavefunction
         COMPLEX(q),POINTER:: CPTWFP(:,:,:) ! wavefunction
         COMPLEX(q),POINTER:: CPROJ(:,:,:)  ! projector (complex)
         REAL(q),   POINTER:: GPROJ(:,:,:)  ! projector (real)
         REAL(q),   POINTER:: FERWE(:,:)    ! fermi-weight for each band
         COMPLEX(q),POINTER:: CELEN(:,:)    ! eigenvalues
         REAL(q),   POINTER:: FERTOT(:,:)   ! global array for fermi-weights
         COMPLEX(q),POINTER:: CELTOT(:,:)   ! global array for eigenvalues
	 LOGICAL           :: OVER_BAND     ! distribution over bands or not
      END TYPE

!
! structure required to store one wavefunctions (no spin no bands)
!
! defined: WUP WDW
      TYPE wavefun1
!only W1
         TYPE(wavedes),POINTER:: WDES1 ! descriptor for wavefunction
         COMPLEX(q),POINTER:: CPTWFP(:)   ! wavefunction
         COMPLEX(q),POINTER:: CPROJ(:)    ! projector (complex)
         REAL(q),   POINTER:: GPROJ(:)    ! projector (real)
         COMPLEX(q),POINTER:: CR(:)       ! wavefunction in real space
         REAL(q)           :: FERWE       ! fermi-weight for each band
         COMPLEX(q)        :: CELEN       ! eigenvalues
      END TYPE


!
! structure required to store a set of wavefunctions including spin
! define: W
!
      TYPE wavespin
!only W
         TYPE(wavedes),POINTER:: WDES   ! descriptor for wavefunction
         COMPLEX(q),POINTER:: CPTWFP(:,:,:,:) ! wavefunction
         COMPLEX(q),POINTER:: CPROJ(:,:,:,:)  ! projector (complex)
         REAL(q),   POINTER:: GPROJ(:,:,:,:)  ! projector (real)
         REAL(q),   POINTER:: FERWE(:,:,:)    ! fermi-weight for each band
         COMPLEX(q),POINTER:: CELEN(:,:,:)    ! eigenvalues
         REAL(q),   POINTER:: FERTOT(:,:,:)   ! global array for fermi-weights
         COMPLEX(q),POINTER:: CELTOT(:,:,:)   ! global array for eigenvalues
	 LOGICAL           :: OVER_BAND       ! distribution over bands or not
      END TYPE























