#!/usr/bin/env bash

##
# This script splits results into matrices and makes illustrations.
##


# Paths
resultFolder=/mnt/work/marc/moba/test_TrioGen/results
postProcessingFolder=/mnt/work/marc/moba/test_TrioGen/results/post_processed
markersInfoFolder=/mnt/archive/MOBAGENETICS/genotypes-base/aux/markerinfo

# Commands
for chr in {1..22} X
do

    echo "extracting p-values for chromosome $chr"

    # p-values
    echo java -Xmx4G -cp bin/triogen-0.3.0-beta/triogen-0.3.0-beta.jar no.uib.triogen.cmd.results.Extract -i $resultFolder/chr_$chr.gz -o $postProcessingFolder/p_chr_$chr -cat phenotype -val variantId,cmf_h_p,h_B1_p,h_B2_p,h_B3_p,h_B4_p,cmf_Bc_p,cmf_Bm_p,cmf_Bf_p,cmf_mt_Bmt_p,cmf_ft_Bft_p

done

for pheno in z_father_bmi z_mother_height z_placenta_weight z_umbilical_chord_length
do

    echo "MH and QQ for $pheno"

    # MH and QQ
    Rscript src/R/buildMHs.R $postProcessingFolder $markersInfoFolder $pheno docs/lm_test "~/R"

done
