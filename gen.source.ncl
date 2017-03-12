#!/bin/bash
#----------vars:nbands-----------------------
if [[ ${modnames[nonmag]} != 1 ]] ; then
 NBANDS=$(($NION/2+$NELECT*3/5))
else
 NBANDS=$(($NION/2+$NELECT/2))
fi
if [[ ${modnames[ncl]} == 1 ]] ; then
    NBANDS=$(($NBANDS*2))
fi
if [[ -z ${keywords[npar]} ]] ; then
 echo "error: gen.source.ncl require that npar be set. it is not. exiting."
 exit 1
fi
if [[ ${modnames[prehf]} == 1 ]] || [[ ${modnames[hse06_l]} == 1 ]] || [[ ${modnames[hse06_h]} == 1 ]] ; then
 NBANDS=$((($NBANDS+${keywords[npar]}-1)/${keywords[npar]}*${keywords[npar]}))
fi
fromvar=NBANDS
tovar=nbands
fromto_comment='nbands value are determined, as enforcing, using the VASP default algorithm read from main.F, in gen.source.ncl.'
source $SCRIPTDIR/gen.source.fromto
#--------vars: lmaxmix---------------------
TEMPVAR=`awk "/Description/,/Error/" POTCAR | awk '{print $1}' | egrep -x '[0-9]+' | sort -rk1 | head -1`
LMAXMIX=$((($TEMPVAR-1)*2))
if [[ $LMAXMIX == 0 ]] ; then
 LMAXMIX=2
fi
fromvar=LMAXMIX
tovar=lmaxmix
fromto_comment='lmaxmix value determined in gen.source.ncl.'
