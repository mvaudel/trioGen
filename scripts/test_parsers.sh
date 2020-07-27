#!/usr/bin/env bash

##
# This script checks that vcf file parsers return the same results.
##

java -cp bin/triogen-0.4.0-beta/triogen-0.4.0-beta.jar no.uib.triogen.test_scripts.TestVcfParsers \
             /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/1.vcf.gz \
             /mnt/work/marc/moba/run/triogen/pheno/trio
