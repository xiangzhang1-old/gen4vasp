#!/bin/bash
# ja412317s: kpoints 8 8 6
# ja412317s: encut: non-mae 450; mae 500; 400 for FM & AFM because otherwise memory usage would rise to 23GB
# ja412317s: ediff (ediffg): non-mae 1E-6 (-0.01); mae 1E-7
# use fast eth: eth_fast corespernode 24
# nelm 80
# nkred 1
# remove symprec
# pre-things
shopt -s extglob
gen_dir="$1"
if [[ ! -f "$1/gen" ]] ; then
  echo "error: $1/gen kpar 4 mempernode 64 not found. exiting."
  exit 64
fi
run_flag="$2"
${gen_dir}/gen kpar 4 npar 1 mempernode 64 generic sname opt_low opt istart 0 icharg 2 prec Low kpoints "8 8 6" isym 0 fm ismetal nsw 60 isif 3 lreal Auto eth_fast corespernode 24 totalnumbercores 24 ldau ediff 1E-4 ediffg 1E-3 
${gen_dir}/gen kpar 4 npar 1 mempernode 64 generic opt sname opt_normal istart 1 icharg 1 prec Normal kpoints "8 8 6" isym 0 fm ismetal nsw 60 isif 3 lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ldau ediff 1E-5 ediffg 1E-4 
${gen_dir}/gen kpar 4 mempernode 64 generic opt ldau sname opt_publevel istart 1 icharg 1 prec Accurate kpoints "8 8 6" isym 0 fm ismetal nsw 60 isif 3 lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-6 ediffg -0.01 npar 1
${gen_dir}/gen  mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Normal kpoints "8 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8  
${gen_dir}/gen mempernode 64 generic hse06_l nelm 80 sname static_hse06_l prec Normal kpoints "8 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-4 nkredx 2 nkredy 2 nkredz 3 bader
${gen_dir}/gen_vasp  "$run_flag"

echo -e "4 \n 1 \n 2 1 1 \n" | vaspkit > /dev/null
mv SC211.vasp POSCAR.2x1x1
mkdir FM AFM
cp POSCAR.2x1x1 FM/POSCAR
cp POSCAR.2x1x1 AFM/POSCAR

cd FM
${gen_dir}/gen kpar 4 mempernode 64 generic opt ldau sname opt_publevel istart 1 icharg 1 prec Accurate kpoints "4 8 6" isym 0 fm ismetal nsw 60 isif 3 lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-6 ediffg -0.01 npar 1
${gen_dir}/gen mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Accurate kpoints "4 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8 nkredx 2 nkredy 4 nkredz 3 
${gen_dir}/gen mempernode 64 generic hse06_l nelm 80 sname static_hse06_l prec Normal kpoints "4 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-6 nkredx 4 nkredy 4 nkredz 3 
${gen_dir}/gen_vasp  "$run_flag"
cd ../

cd AFM
${gen_dir}/gen kpar 4 mempernode 64 generic opt ldau sname opt_publevel istart 1 icharg 1 prec Accurate kpoints "4 8 6" isym 0 afm ismetal nsw 60 isif 3 lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-6 ediffg -0.01 npar 1
${gen_dir}/gen mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Accurate kpoints "4 8 6" isym 0 afm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8 nkredx 2 nkredy 4 nkredz 3 
${gen_dir}/gen mempernode 64 generic hse06_l nelm 80 sname static_hse06_l prec Normal kpoints "4 8 6" isym 0 afm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-6 nkredx 4 nkredy 4 nkredz 3 
${gen_dir}/gen_vasp  "$run_flag"
cd ..

mkdir 100 110 001
cp POSCAR 100
cp POSCAR 110
cp POSCAR 001

cd 100
${gen_dir}/gen mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Normal kpoints "8 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8 nkredx 2 nkredy 4 nkredz 3 
${gen_dir}/gen mempernode 64 generic hse06_h nelm 80 sname static_hse06_h_ncl prec Normal kpoints "8 8 6"  fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-7 ncl saxis "1 0 0" 
${gen_dir}/gen_vasp  "$run_flag"
cd ..

cd 110
${gen_dir}/gen mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Normal kpoints "8 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8 nkredx 2 nkredy 4 nkredz 3 
${gen_dir}/gen mempernode 64 generic hse06_h nelm 80 sname static_hse06_h_ncl prec Normal kpoints "8 8 6"  fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-7 ncl saxis "1 1 0"
${gen_dir}/gen_vasp  "$run_flag"
cd ..

cd 001
${gen_dir}/gen mempernode 64 generic sname static_prehf prehf istart 0 icharg 2 prec Normal kpoints "8 8 6" isym 0 fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-8 nkredx 2 nkredy 4 nkredz 3 
${gen_dir}/gen mempernode 64 generic hse06_h nelm 80 sname static_hse06_h_ncl prec Normal kpoints "8 8 6"  fm ismetal lreal .FALSE. eth_fast corespernode 24 totalnumbercores 24 ediff 1E-7 ncl saxis "0 0 1" 
${gen_dir}/gen_vasp  "$run_flag"
cd ..
