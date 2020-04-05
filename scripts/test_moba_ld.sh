#!/usr/bin/env bash

##
# This script creates ld matrices on all MOBA.
##


# Paths
vcfFolder=~/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf
trioFile=~/mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919
outputFolder=~/mnt/work/marc/moba/test_TrioGen/ld

# Variables

# Commands
for chr in 22
do

    echo "Processing chromosome $chr"

    java -Xmx124G -cp bin/triogen-0.4.0-beta/triogen-0.4.0-beta.jar no.uib.triogen.cmd.ld.LdMatrix -g $vcfFolder/$chr.vcf.gz -gf 1 -f $trioFile -o $outputFolder/chr_$chr -nv 32

done
