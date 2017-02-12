!
! communication description include file
! this structure allows the use a several communicators
! it currently supports MPI and CRAY shmem calls
!
! Nbranch is the number of branches at each node in the gsum routines
! two should be always fine
!
      INTEGER Nbranch
      PARAMETER( Nbranch=2 )

     TYPE communic
!only COMM
        INTEGER MPI_COMM          ! MPI_Communicator
        INTEGER NODE_ME           ! node id starting from 1 ... NCPU
        INTEGER IONODE            ! node which has to do IO (set to 0 for no IO)
        INTEGER NCPU              ! total number of proc in this communicator
! following variables are used for reductions if MPI_allreduce is not used
        INTEGER parent            ! Id of our parent in the global tree
        INTEGER nchild            ! number of children in the global tree
        INTEGER childs(Nbranch)   ! Ids of our children in the global tree
! following variables are used only within MPI (private data)
        INTEGER,POINTER :: orders(:)  ! order of pair exchanges in the global exchange
        INTEGER hccrd             ! coordinate of this node in the hypercube: 1..hcnds (NODE_ME)
        INTEGER hcdim             ! dimension of the hypercube
        INTEGER hcnds             ! number of nodes in the hypercube (NCPU)
        INTEGER hcroot            ! hypercube coordinate of root node
        INTEGER,POINTER :: hid(:) ! mapping of ids (i.e. virtual node number -> real node id)
     END TYPE
