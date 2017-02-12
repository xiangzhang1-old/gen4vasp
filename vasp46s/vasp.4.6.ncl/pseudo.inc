!
! pseudo potential discription include file
! only included if MODULES are not supported
!only P
      INTEGER, PARAMETER:: NEKERR=100
      INTEGER, PARAMETER:: NPSPTS=1000
      INTEGER, PARAMETER:: NPSNL=100,NPSRNL=100
!
! structure to support radial grid
!
      TYPE potcar
        REAL(q) ZVALF,POMASS,RWIGS  ! valence, mass, wigner seitz radius
        REAL(q) ZVALF_ORIG          ! original valence
                                    ! ZVAL might be reset when core level shifts are calculated
        REAL(q) EATOM               ! atomic energy
        REAL(q) EGGA(4)             ! atomic energy for GGA
        REAL(q) ENMAXA,ENMINA,EAUG  ! energy cutoffs
        REAL(q) QOPT1,QOPT2         ! real space optimized for this cutoffs
        REAL(q) EKECUT(NEKERR),EKEERR(NEKERR) ! description of cutofferrors
        INTEGER LEXCH   	    ! exchange type
        REAL(q) PSGMAX              ! maximal G for local potential
        REAL(q) PSMAXN              ! maximal G for non local potential
        REAL(q) PSRMAX              ! maximal r for non local contrib.
        REAL(q) PSDMAX              ! maximal r for augmentation charge
        REAL(q),POINTER :: PSP(:,:) ! local pseudopotential in rec. space
        REAL(q),POINTER :: PSPCOR(:)! partial core information in rec. space
        REAL(q),POINTER :: PSPRHO(:)! atomic pseudo charge density in rec. space
        REAL(q),POINTER :: PSPNL(:,:)  ! non local proj. rec. space
        REAL(q),POINTER :: PSPRNL(:,:,:) ! non local proj. real space
        REAL(q),POINTER :: DION(:,:)! non local strength
        REAL(q),POINTER :: QION(:,:)! spherical augmentation charge
        REAL(q),POINTER :: QTOT(:,:)! total charge in each channel
        REAL(q),POINTER :: QPAW(:,:,:)! aspherical full augmentation charge
        REAL(q),POINTER :: QATO(:,:)! initial occupancies (in atom)
        REAL(q),POINTER :: QDEP(:,:,:) ! augmentation charge
        REAL(q),POINTER :: NABLA(:,:,:)! atomic augmentation matrix elements of nabla operator
        REAL(q) PSCORE              ! is equal to V(q) + 4 Pi / q^2
        INTEGER,POINTER :: LPS(:)   ! L quantum number for each PP
        INTEGER,POINTER :: NLPRO(:) ! unused
        INTEGER,POINTER :: NDEP(:,:)! number of augmentation channels per ll'
      ! quantities defined on radial grid
        TYPE (rgrid)    :: R        ! radial grid
	REAL(q),POINTER :: POTAE(:) ! froze core AE potential on r-grid (valence only)
	REAL(q),POINTER :: POTPS(:) ! local PP on r-grid (valence only)
	REAL(q),POINTER :: POTPSC(:)! local PP on r-grid (core only)
	REAL(q),POINTER :: RHOAE(:) ! frozen core charge rho(r)r^2 on r-grid
	REAL(q),POINTER :: TAUAE(:) ! kinetic energy density of core
	REAL(q),POINTER :: RHOPS(:) ! frozen core pseudo charge rho(r)r^2 on r-grid
	REAL(q),POINTER :: WAE(:,:) ! ae valence wavefunction on r-grid
	REAL(q),POINTER :: WPS(:,:) ! pseudo valence wavefunction on r-grid
	REAL(q),POINTER :: AUG(:,:) ! L-dependent augmentation charge on r-grid
	REAL(q)         :: DEXCCORE ! exchange correlation energy of frozen core
	REAL(q)        :: DEXCCOREM ! exchange correlation energy of frozen core (MetaGGA)
      ! most integers are on the end (alignment ?)
        INTEGER      LDIM           ! dimension for l channels
        INTEGER      LMAX           ! l-channels for non local PP
        INTEGER      LMDIM          ! dimension for lm channels
        INTEGER      LMMAX          ! lm-channels for non local PP
        INTEGER      LDIM2          ! dimension for augmentation arryas
        LOGICAL      LREAL          ! real space optimized ?
	LOGICAL      LUNSCR         ! partial core has been unscreened
        INTEGER      LMAX_CALC      ! maximum L for onsite terms in PAW
	     ! might be overwritten by LMAXPAW line in the INCAR file
        REAL(q)      ZCORE          ! charge of core
        CHARACTER*40 SZNAMP         ! header
        CHARACTER*2  ELEMENT        ! Name of element
      END TYPE

