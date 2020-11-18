#!/usr/bin/env bash

##
# This script makes a test run on the moba data.
##

java -Xmx32G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.cmd.association.LinearModel \
             -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/23.vcf.gz \
             -gf 1 \
             -vi /mnt/work/marc/moba/mobaRun/resources/triogen/targets/targets_h2020 \
             -f /mnt/work/marc/moba/run/triogen/pheno/trio \
             -p /mnt/work/marc/moba/run/triogen/pheno/trio_pheno \
             -pn "z_bmi0" \
             -cg "child_harvest,child_Rotterdam1,child_Rotterdam2,child_NormentMay16,child_NormentFeb18,child_Ted,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
             -cs /mnt/work/marc/moba/run/triogen/pheno/standardization_details.json \
             -maf 0.005 \
             -dos \
             -o /mnt/work/marc/moba/run/triogen/debug/test_X.lm_target.gz \
             -nv 12 \
             -vl
