!======================================================================
! This file contains PARAMETER constants and COMMON variables
! describing the virtual hypercube structure of the parallel machine.
!
! Dependencies: NONE
!----------------------------------------------------------------------
! PARAMETER constant section

! NDIMEN is the number of dimensions in the hypercube. Change this to
!        on more or less numbers of processors.
!
! NODES  is the number of nodes in the hypercube. Thus is set to the
!        value of the expression 2**NDIMEN.

      INTEGER NDIMEN, NODES

      PARAMETER (NDIMEN = 5)
      PARAMETER (NODES = 2**NDIMEN)

!----------------------------------------------------------------------
! COMMON variable section

! Block /CUBE/ holds node specific cube data.
!
! INODE is the number of this node in the range 1:NODES. It should be
!       assigned in the main program before calling GENSET etc.
!
! IONODE is the number of the node that does I/O for the whole program
!

      INTEGER INODE, IONODE

      COMMON/CUBE/INODE, IONODE
      INTEGER nbandtot_com,band_index
      common/bands/nbandtot_com,band_index
      INTEGER nrplwv_node_com,nbands_com,lmdnio_node_com
      common/rearrange/nrplwv_node_com,nbands_com,lmdnio_node_com
