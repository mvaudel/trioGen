#!/usr/bin/env bash

##
# This script runs triogen on all MOBA.
##


# Paths
vcfFolder=/mnt/archive/MOBAGENETICS/genotypes-base/imputed/all/vcf
trioFile=/mnt/archive/helgeland/helgeland2020/sensitive/fulltriads-n24280-180919
phenoFile=/mnt/work/marc/moba/test_TrioGen/pheno/phenos
phenoFile=/mnt/work/marc/moba/test_TrioGen/pheno/phenos
outputFolder=/mnt/work/marc/moba/test_TrioGen/results

# Variables
phenoName=breastmilk_duration,formula_freq_6m,pregnancy_duration,z_umbilical_chord_length,z_placenta_weight,z_bmi0,z_bmi1,z_bmi2,z_bmi3,z_bmi4,z_bmi5,z_bmi6,z_bmi7,z_bmi8,z_bmi9,z_bmi10,z_bmi11,z_mother_height,z_father_bmi
covariates=harvest,rotterdam1,rotterdam2,normentMay16,normentMay18,ted,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

# Commands
for chr in {1..7}
do

    echo "Processing chromosome $chr"

    java -Xmx32G -cp bin/triogen-0.3.0-beta/triogen-0.3.0-beta.jar no.uib.triogen.cmd.association.LinearModel -g $vcfFolder/$chr.vcf.gz -maf 0.05 -gf 1 -f $trioFile -p $phenoFile -pn $phenoName -cv $covariates -o $outputFolder/chr_$chr.gz

done

echo "Processing chromosome X"

java -Xmx32G -cp bin/triogen-0.3.0-beta/triogen-0.3.0-beta.jar no.uib.triogen.cmd.association.LinearModel -g $vcfFolder/X.vcf.gz -maf 0.05 -gf 1 -f $trioFile -p $phenoFile -pn $phenoName -cv $covariates -o $outputFolder/chr_X.gz
