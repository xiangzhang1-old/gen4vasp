	parameter (maxpet=1000)
c
        real*8 mflop
        real*8 cputot,cpu0
        real*8 vptime,cptime
        integer numcal
        character*32  apetname
        real*4 etime_, tarray(2)
c
        common /SNI_time_1/ mflop(maxpet),cputot(maxpet),
     &                      cpu0(maxpet),numcal(maxpet)
        common /SNI_time_2/ apetname(maxpet)
