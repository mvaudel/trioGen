#!/usr/bin/env bash

##
# This script splits results into matrices that can be fed to R.
##


# Paths
resultFolder=~/mnt/work/marc/moba/test_TrioGen/results

# Commands
for chr in {8}
do

    echo "Splitting chromosome $chr"

    # p-values
    java -Xmx32G -cp bin/triogen-0.3.0-beta/triogen-0.3.0-beta.jar no.uib.triogen.cmd.results.Extract -i $resultFolder/chr_$chr.gz -o $resultFolder/h_cmf_p_chr_$chr -cat phenotype -val cmf_h_p,h_B1_p,h_B2_p,h_B3_p,h_B4_p,cmf_Bc_p,cmf_Bm_p,cmf_Bf_p,cmf_mt_Bmt_p,cmf_ft_Bft_p

done
