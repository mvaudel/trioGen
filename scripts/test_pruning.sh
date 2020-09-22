#!/usr/bin/env bash

##
# This script makes a test run on the moba data.
##

java -Xmx5G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.cmd.ld_pruning.LdPruning \
             -res /mnt/work/marc/moba/run/bolt/bolt_output/childGeno_z_bmi3-stats-bgen.gz \
			 -l /mnt/work/marc/moba/run/triogen/ld \
			 -id SNP \
			 -pn P_BOLT_LMM \
			 -cn CHR \
			 -o /mnt/work/marc/moba/run/bolt/bolt_output/pruned/childGeno_z_bmi3.gz