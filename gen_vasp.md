Purpose:    "gen" generates INCAR-opt_hse06_001, which is too complicated. Therefore, each generated file-group goes into gen.log, and a gen-vasp command does them all. 

Workflow:   create a NOW folder
            loop [read $filename from gen.log, remove tag_* and create tag_$filename file inside NOW, cp input file inside, calculate, cp -r to $filename]
            If any command fails, gen_vasp will quit.

Usage: gen_vasp $RUN_LOG_DIR $RUN_FLAG(test/ustccloud/ustcscc)
