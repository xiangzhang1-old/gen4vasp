!
!  structure required to support non local projection operators in recip space
!
      TYPE nonl_struct
!only NONL_S
         LOGICAL LRECIP                ! structure set up ?
         INTEGER NTYP                  ! number of types
         INTEGER NIONS                 ! number of ions
         INTEGER NK                    ! kpoint for which CREXP is set up
         INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
         INTEGER, POINTER :: LMMAX(:)  ! max l-quantum number for each type
!-MM- spin spiral stuff
         LOGICAL LSPIRAL               ! do we want to calculate spin spirals
! changes to accomodate spin spirals
! original statement
!        REAL(q), POINTER ::QPROJ(:,:,:,:) ! projectors in reciprocal space for each k-point
         REAL(q), POINTER ::QPROJ(:,:,:,:,:) ! projectors in reciprocal space for each k-point
!-MM-
         COMPLEX(q),POINTER ::CREXP(:,:)  ! phase factor exp (i (G+k) R(ion))
         REAL(q),  POINTER  ::POSION(:,:) ! positions (required for setup)
         COMPLEX(q),POINTER ::CQFAK(:,:)  ! i^l
      END TYPE



