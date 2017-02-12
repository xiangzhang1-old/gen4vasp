define(ENTRY,
         global    $1
         global    $1_
$1:
$1_:)
define(RETURN,
         b         0(,$s32))
#
         text
#
#        itmtv()  - vector execution clock.
#
#        Exit:
#        ($s123) = 1/overflow flag bit,15/0, 48/vector execution clock.
#
         ENTRY(itmtv)
         smir    $s123,0
         RETURN
