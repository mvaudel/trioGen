#!/usr/bin/env bash

##
# This script runs linear association on targets in MOBA.
##


# Paths
vcfFolder=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf
trioFile=/mnt/work2/helgeland/helgeland2020/sensitive/290520/full-triads-sentrixids-n23706-290520
phenoFile=/mnt/work/marc/moba/test_TrioGen/pheno/phenos
outputFolder=/mnt/work/marc/moba/test_TrioGen/targets

# Variables

# Commands
for chr in 7 10 9 6 5 3 2 1
do

    echo "Processing targets in chromosome $chr"

    java -Xmx64G -cp bin/triogen-0.4.0-beta/triogen-0.4.0-beta.jar no.uib.triogen.cmd.association.LinearModel -g $vcfFolder/$chr.vcf.gz -gf 1 -vi resources/targets_chr$chr -maf 0.05 -f $trioFile -p $phenoFile -pn breastmilk_duration,formula_freq_6m,pregnancy_duration,z_umbilical_chord_length,z_placenta_weight,z_bmi0,z_bmi1,z_bmi2,z_bmi3,z_bmi4,z_bmi5,z_bmi6,z_bmi7,z_bmi8,z_bmi9,z_bmi10,z_bmi11,z_mother_height,z_father_bmi -cv harvest,rotterdam1,rotterdam2,normentMay16,normentMay18,ted,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 -o $outputFolder/$chr.lm_target.gz -nv 2 -d 1000000

done
