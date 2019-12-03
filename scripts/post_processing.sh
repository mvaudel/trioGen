#!/usr/bin/env bash

##
# This script splits results into matrices and makes illustrations.
##


# Paths
resultFolder=~/mnt/work/marc/moba/test_TrioGen/results

# Commands
for chr in {8..9}
do

    echo "extracting p-values for chromosome $chr"

    # p-values
    java -Xmx32G -cp bin/triogen-0.3.0-beta/triogen-0.3.0-beta.jar no.uib.triogen.cmd.results.Extract -i $resultFolder/chr_$chr.gz -o $resultFolder/p_chr_$chr -cat phenotype -val cmf_h_p,h_B1_p,h_B2_p,h_B3_p,h_B4_p,cmf_Bc_p,cmf_Bm_p,cmf_Bf_p,cmf_mt_Bmt_p,cmf_ft_Bft_p


    echo "Building MHs and QQs for $chr"



done
