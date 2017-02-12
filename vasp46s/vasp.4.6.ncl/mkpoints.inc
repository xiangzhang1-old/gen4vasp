!
!  structure required for kpoints generation
!
      TYPE kpoints_struct
!only  KPOINTS
        INTEGER NKDIM                ! maximal number of k-points
        INTEGER NKPTS                ! actual number of k-points
        REAL(q),POINTER:: VKPT(:,:)     ! coordinate of k-point
        REAL(q),POINTER:: WTKPT(:)      ! symmetry weight-factor for each k-point
        REAL(q)        :: VOLWGT        ! volume weight for each tetrahedron
        INTEGER,POINTER:: IDTET(:,:) ! link list for tetrahedron
        INTEGER     :: NTET          ! number of tetrahedrons
        INTEGER,POINTER:: IKPT(:,:,:) ! temporary work array for tetrahedron method
        LOGICAL     :: LTET          ! use tetrahedron method ?
        INTEGER ISMEAR               ! type of smearing
        REAL(q)    SIGMA                ! type of smearing
        REAL(q)    EMIN                 ! minimal E for DOS
        REAL(q)    EMAX                 ! maximal E for DOS
        CHARACTER*40  SZNAMK         ! name of k-points file
      END TYPE
