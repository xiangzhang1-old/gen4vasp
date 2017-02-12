==Idea==
The script is based on the notion of:
 - modules [ changeable_keyword, nonchangeable_keyword, customized_keyword, othercustomized_keyword ]
 - completeness: in each module, all publevel keywords are compulsorily set by hand, or defaulted to the mandatory value. if something goes unset, you'll definitely know
 - explicit self-consistency: no overwriting is ever allowed. nonchangeable keyword values are fixed, and changeable keyword values cannot be rewritten. if conflicts arrise, you'll definitely now
 - partial logical self-consistency: some logical errors (i.e. isnotmetal+ismetal,hse06+isym=2) is checked
 - portability: through gen.conf.modules,gen.conf.cust,gen.con.othercust,gen.conf.opt_man_modules,gen.conf.kwban,gen.conf.kwnullify, you get to establish, with ease and from scratch/newbie, a module-based command-line quantumchem-input-file-establishing system. Any of them!


==Data structure==
Associative arrays:
modnames[$modname=1/0] #flags   #mod "virtual" collects ckw's and nckw's from cust and othercust. note that module "virtual"'s kw's values are not used at all. only names are used.
$modname_ckw[$kwname=""]        #changeable_keyword                 #are necessary in module, but value changeable, therefore read or overwritten by nckw. still, empty at the start
$modname_nckw[$kwname=$kwval]   #nonchangeable_keyword and value    #are necessary in module, and value fixed
$keywords[$kwname=$kwval]
$nullifiedkw[$kwname=$kwval]    #effective after keyword nullifying takes place

Files: 
gen.conf.modules        =>  module infrastructure (& with cust,othercust...)    =>  $modname:$nckwname=$nckwval,...:$ckwname
gen.conf.cust           =>  customized keywords if modules $ifyes are present   =>  $ifyes:$ifnot:$nckwname=$nckwval,...:$ckwname:1/0   #1 if must, 0 if recommended. 1 comes first!
gen.conf.othercust      =>  script-customized kws if modules $ifyes are present =>  $ifyes:$ifnot:$nckwname:$ckwname[1]:$batchname,...  #$batchname is run if $ifyes !$ifnot
..opt_mand_modules      =>  customizable optional mandatory modules             =>  $sum_act_mod_number1,...:$modname1,...              #active mods summed equals $sum_act_1 or ...
----
[1] for othercust, nckw means 'keywords whose values must be dictated by othercust' and ckw means 'keywords whose values must survive as non-zero but can be changed by anyone'


==Additional feature==
gen.conf.kwban          =>  ban certain value of keywords (or kws) if $ifyes    =>  $ifyes:$ifnot:$kwban_kw=$kwban_kwval:$kwban_kw
gen.conf.kwnullify      =>  nullify certain kw if $ifyes !$ifnot                =>  $ifyes:$ifnot:$kwnullify
gen.conf.check          =>  perform non-complete syntax&semantic check          =>  $ifyes:$ifnot:$batchname,...                        #$batchname is run if $ifyes !$ifnot
gen.conf.postcheck      =>  check after KPOINTS/INCAR/POTCAR writing            =>  $ifyes:$ifnot:$batchname,...                        #$batchname is run if $ifyes !$ifnot

==Usage==
gen [modname] [kwname] [kwval]
