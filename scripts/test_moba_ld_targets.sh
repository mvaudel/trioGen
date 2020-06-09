#!/usr/bin/env bash

##
# This script creates ld matrices on targets in MOBA.
##


# Paths
vcfFolder=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf
trioFile=/mnt/work2/helgeland/helgeland2020/sensitive/290520/full-triads-sentrixids-n23706-290520
outputFolder=/mnt/work/marc/moba/test_TrioGen/ld

# Variables

# Commands
for chr in 7 10 9 6 5 3 2 1
do

    echo "Processing chromosome $chr"

    java -Xmx100G -cp bin/triogen-0.4.0-beta/triogen-0.4.0-beta.jar no.uib.triogen.cmd.ld.LdMatrix -g $vcfFolder/$chr.vcf.gz -gf 1 -vi resources/targets_chr$chr -f $trioFile -o $outputFolder/dos_chr_$chr -nv 2 -d 2000000

done