#!/usr/bin/env bash

##
# This script extracts summary statistics on the number of Mendelian errors in phased trio genotypes.
##

java -cp bin/triogen-0.4.0-beta/triogen-0.4.0-beta.jar no.uib.triogen.cmd.mendelian_check.MendelianCheck \
             -g /mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf/22.vcf.gz \
             -gf 1 \
             -f /mnt/work/marc/moba/run/triogen/pheno/trio \
             -maf 0.005 \
             -o /mnt/work/marc/moba/run/triogen/m_err/chr22 \
             -nv 42
