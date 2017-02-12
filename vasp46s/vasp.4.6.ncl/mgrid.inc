!
! grid description input file
! only included if MODULES are not supported
! defined: GRID, GRIDC, GRIDC_UNB, GRIDB
      TYPE layout
        INTEGER NP                   ! total number of points
        INTEGER NALLOC               ! allocation required
        INTEGER NFAST                ! which index is fast (1-x, 2-y, 3-z)
        INTEGER NCOL                 ! number of columns in the grid
        INTEGER NROW                 ! number of elements in each column
        INTEGER,POINTER :: I2(:)     ! y/z/x-index of each column
        INTEGER,POINTER :: I3(:)     ! z/x/y-index of each column
        INTEGER,POINTER :: INDEX(:,:)! column index for each yz, zx or xy pair
      END TYPE

      TYPE grid_map
        INTEGER,POINTER :: N(:)      ! number of elements send by each node
        INTEGER,POINTER :: PTR(:)    ! sum_j=1,I N(j)
        INTEGER,POINTER :: RMT(:)    ! remote address for send  (shmem t3d)
        INTEGER,POINTER :: TBL(:)    ! address of each element
      ! inverse transformation (i.e. receiver information)
        INTEGER,POINTER :: NI(:)
        INTEGER,POINTER :: PTRI(:)
        INTEGER,POINTER :: TBLI(:)
        INTEGER,POINTER :: RMTI(:)
	LOGICAL LOCAL                ! all information is local
        LOGICAL LOCAL_COPY           ! no data redistribution required
      END TYPE

      TYPE grid_3d
!only  GRID
        INTEGER NGX,NGY,NGZ          ! number of grid points in x,y,z
	INTEGER NGX_rd,NGY_rd,NGZ_rd ! in the complex mode the _rd values
               ! are equal to NGX, NGY, NGZ
	       ! if a real to complex FFT is used only half of the data
               ! are stored in one direction and the corresponding NG?_rd is 
	       ! set to (NG?+1)/2
        INTEGER NPLWV                ! total number of grid points (full grid)
        INTEGER MPLWV                ! allocation in complex words
        INTEGER NGPTAR(3)            ! equivalent to /(NGX,NGY,NGZ/)
        INTEGER,POINTER :: LPCTX(:)  ! loop counters in x
        INTEGER,POINTER :: LPCTY(:)  ! loop counters in y
        INTEGER,POINTER :: LPCTZ(:)  ! loop counters in z
! loop counters, in which the unbalanced contribution is zeroed
        INTEGER,POINTER :: LPCTX_(:) ! loop counters in x
        INTEGER,POINTER :: LPCTY_(:) ! loop counters in y
        INTEGER,POINTER :: LPCTZ_(:) ! loop counters in z
! reciprocal space layout (x is always  fast index)
        TYPE(layout)    :: RC
! intermediate layout (y is always  fast index, used only in parallel version)
        TYPE(layout)    :: IN
! real space layout   (x or z is the fast index)
        TYPE(layout)    :: RL
! information only required for real space representation
        INTEGER NGZ_complex          ! number of grid points for z fast
! mapping for parallel version
        TYPE(grid_map)  :: RC_IN     ! recip -> intermediate
        TYPE(grid_map)  :: IN_RL     ! intermediate -> real space
        TYPE(communic), POINTER :: COMM
	LOGICAL         :: LREAL     ! are data stored as complex or real number 
                                     ! in the direct space ?
      END TYPE


!
! transition table used to go from a large to a small grid
! or vice versa
!
      TYPE transit
        INTEGER,POINTER :: IND1(:)   ! fast index transition table
        INTEGER,POINTER :: INDCOL(:) ! column to column transition table
      END TYPE

