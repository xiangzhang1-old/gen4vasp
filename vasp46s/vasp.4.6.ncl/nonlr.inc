!
!  structure required to support non local projection operators in real space
!
      TYPE nonlr_proj
         REAL(q), POINTER :: PSPRNL(:,:,:)
         INTEGER, POINTER :: LPS(:)
      END TYPE nonlr_proj
      
      TYPE nonlr_struct
!only NONLR_S
         LOGICAL LREAL                 ! structure set up ?
         INTEGER NTYP                  ! number of types
         INTEGER NIONS                 ! number of ions
         INTEGER IRMAX                 ! maximum number points in sphere
         INTEGER IRALLOC               ! size for allocation =IRMAX*LMDIM*NIONS
         INTEGER NK                    ! kpoint for which CRREXP is set up
         INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
         INTEGER, POINTER :: LMAX(:)   ! max l-quantum number for each type
         INTEGER, POINTER :: LMMAX(:)  ! number lmn-quantum numbers for each type
         INTEGER, POINTER ::CHANNELS(:)! number of ln-quantum for each type
         REAL(q), POINTER :: PSRMAX(:) ! real space cutoff
         INTEGER, POINTER :: NLIMAX(:) ! maximum index for each ion
         INTEGER, POINTER :: NLI(:,:)  ! index for gridpoints
         REAL(q) ,POINTER :: RPROJ(:)  ! projectors on real space grid
         REAL(qs),POINTER :: RPROJS(:) ! projectors on real space grid (single precission)
         REAL(q), POINTER :: POSION(:,:)! positions (required for setup)
         TYPE(nonlr_proj), POINTER :: BETA(:) ! a set of structures containing pointers to 
                                       ! the parameters for the non local projection operators
!-MM- spin spiral stuff
         LOGICAL LSPIRAL               ! do we want to calculate spin spirals
! change in phaser array to accomodate spin spirals
! original statement
!        COMPLEX(q),POINTER:: CRREXP(:,:) ! phase factor exp (i k (R(ion)-r(grid)))
         COMPLEX(q),POINTER:: CRREXP(:,:,:) ! phase factor exp (i k (R(ion)-r(grid)))
!-MM- end of alteration
      END TYPE
