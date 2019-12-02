# Chromosome 7 - z_BMI


Chromosome 7 run with TrioGen v.0.3.0-beta on 27,451 full trios, ADHD cases, related, and ethnic outliers excluded. 10 PCs and genotyping batch as covariates. Reference panel AF > 5%, info score >= 0.7, singularities excluded

z_BMI “à la Chris”, see details [here](../pheno/plots.md):

- not the same as the pheno tables used so far!

- adjusted for pregnancy duration!


4 Models:

- h: `y = ß1 h1 + ß2 h2 + ß3 h3 + ß4 h4  + e`

- cmf: `y = ßm (h1 + h2) + ßc (h1 + h3) + ßf (h3 + h4) + e`

- cmf_mt: `y = ßm (h1 + h2) + ßc (h1 + h3) + ßf (h3 + h4) + ßmt h1 + e`

- cmf_ft: `y = ßm (h1 + h2) + ßc (h1 + h3) + ßf (h3 + h4) + ßft h1 + e`


### Birth - h vs. cmf (F-test p-value)

![](z_bmi0_cmf_h_p_MH.png)

![](z_bmi0_cmf_h_p_QQ.png)


### Birth - h betas (Prob(|t| > 0))

- B1

![](z_bmi0_h_B1_p_MH.png)

![](z_bmi0_h_B1_p_QQ.png)

- B2

![](z_bmi0_h_B2_p_MH.png)

![](z_bmi0_h_B2_p_QQ.png)

- B3

![](z_bmi0_h_B3_p_MH.png)

![](z_bmi0_h_B3_p_QQ.png)

- B4

![](z_bmi0_h_B4_p_MH.png)

![](z_bmi0_h_B4_p_QQ.png)


### Birth - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi0_cmf_Bc_p_MH.png)

![](z_bmi0_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi0_cmf_Bm_p_MH.png)

![](z_bmi0_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi0_cmf_Bf_p_MH.png)

![](z_bmi0_cmf_Bf_p_QQ.png)


### Birth - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi0_cmf_mt_Bmt_p_MH.png)

![](z_bmi0_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi0_cmf_ft_Bft_p_MH.png)

![](z_bmi0_cmf_ft_Bft_p_QQ.png)


### 6 w - h vs. cmf (F-test p-value)

![](z_bmi1_cmf_h_p_MH.png)

![](z_bmi1_cmf_h_p_QQ.png)


### 6 w - h betas (Prob(|t| > 0))

- B1

![](z_bmi1_h_B1_p_MH.png)

![](z_bmi1_h_B1_p_QQ.png)

- B2

![](z_bmi1_h_B2_p_MH.png)

![](z_bmi1_h_B2_p_QQ.png)

- B3

![](z_bmi1_h_B3_p_MH.png)

![](z_bmi1_h_B3_p_QQ.png)

- B4

![](z_bmi1_h_B4_p_MH.png)

![](z_bmi1_h_B4_p_QQ.png)


### 6 w - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi1_cmf_Bc_p_MH.png)

![](z_bmi1_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi1_cmf_Bm_p_MH.png)

![](z_bmi1_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi1_cmf_Bf_p_MH.png)

![](z_bmi1_cmf_Bf_p_QQ.png)


### 6 w - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi1_cmf_mt_Bmt_p_MH.png)

![](z_bmi1_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi1_cmf_ft_Bft_p_MH.png)

![](z_bmi1_cmf_ft_Bft_p_QQ.png)


### 3 m - h vs. cmf (F-test p-value)

![](z_bmi2_cmf_h_p_MH.png)

![](z_bmi2_cmf_h_p_QQ.png)


### 3 m - h betas (Prob(|t| > 0))

- B1

![](z_bmi2_h_B1_p_MH.png)

![](z_bmi2_h_B1_p_QQ.png)

- B2

![](z_bmi2_h_B2_p_MH.png)

![](z_bmi2_h_B2_p_QQ.png)

- B3

![](z_bmi2_h_B3_p_MH.png)

![](z_bmi2_h_B3_p_QQ.png)

- B4

![](z_bmi2_h_B4_p_MH.png)

![](z_bmi2_h_B4_p_QQ.png)


### 3 m - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi2_cmf_Bc_p_MH.png)

![](z_bmi2_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi2_cmf_Bm_p_MH.png)

![](z_bmi2_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi2_cmf_Bf_p_MH.png)

![](z_bmi2_cmf_Bf_p_QQ.png)


### 3 m - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi2_cmf_mt_Bmt_p_MH.png)

![](z_bmi2_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi2_cmf_ft_Bft_p_MH.png)

![](z_bmi2_cmf_ft_Bft_p_QQ.png)


### 6 m - h vs. cmf (F-test p-value)

![](z_bmi3_cmf_h_p_MH.png)

![](z_bmi3_cmf_h_p_QQ.png)


### 6 m - h betas (Prob(|t| > 0))

- B1

![](z_bmi3_h_B1_p_MH.png)

![](z_bmi3_h_B1_p_QQ.png)

- B2

![](z_bmi3_h_B2_p_MH.png)

![](z_bmi3_h_B2_p_QQ.png)

- B3

![](z_bmi3_h_B3_p_MH.png)

![](z_bmi3_h_B3_p_QQ.png)

- B4

![](z_bmi3_h_B4_p_MH.png)

![](z_bmi3_h_B4_p_QQ.png)


### 6 m - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi3_cmf_Bc_p_MH.png)

![](z_bmi3_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi3_cmf_Bm_p_MH.png)

![](z_bmi3_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi3_cmf_Bf_p_MH.png)

![](z_bmi3_cmf_Bf_p_QQ.png)


### 6 m - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi3_cmf_mt_Bmt_p_MH.png)

![](z_bmi3_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi3_cmf_ft_Bft_p_MH.png)

![](z_bmi3_cmf_ft_Bft_p_QQ.png)


### 8 m - h vs. cmf (F-test p-value)

![](z_bmi4_cmf_h_p_MH.png)

![](z_bmi4_cmf_h_p_QQ.png)


### 8 m - h betas (Prob(|t| > 0))

- B1

![](z_bmi4_h_B1_p_MH.png)

![](z_bmi4_h_B1_p_QQ.png)

- B2

![](z_bmi4_h_B2_p_MH.png)

![](z_bmi4_h_B2_p_QQ.png)

- B3

![](z_bmi4_h_B3_p_MH.png)

![](z_bmi4_h_B3_p_QQ.png)

- B4

![](z_bmi4_h_B4_p_MH.png)

![](z_bmi4_h_B4_p_QQ.png)


### 8 m - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi4_cmf_Bc_p_MH.png)

![](z_bmi4_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi4_cmf_Bm_p_MH.png)

![](z_bmi4_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi4_cmf_Bf_p_MH.png)

![](z_bmi4_cmf_Bf_p_QQ.png)


### 8 m - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi4_cmf_mt_Bmt_p_MH.png)

![](z_bmi4_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi4_cmf_ft_Bft_p_MH.png)

![](z_bmi4_cmf_ft_Bft_p_QQ.png)


### 1 y - h vs. cmf (F-test p-value)

![](z_bmi5_cmf_h_p_MH.png)

![](z_bmi5_cmf_h_p_QQ.png)


### 1 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi5_h_B1_p_MH.png)

![](z_bmi5_h_B1_p_QQ.png)

- B2

![](z_bmi5_h_B2_p_MH.png)

![](z_bmi5_h_B2_p_QQ.png)

- B3

![](z_bmi5_h_B3_p_MH.png)

![](z_bmi5_h_B3_p_QQ.png)

- B4

![](z_bmi5_h_B4_p_MH.png)

![](z_bmi5_h_B4_p_QQ.png)


### 1 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi5_cmf_Bc_p_MH.png)

![](z_bmi5_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi5_cmf_Bm_p_MH.png)

![](z_bmi5_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi5_cmf_Bf_p_MH.png)

![](z_bmi5_cmf_Bf_p_QQ.png)


### 1 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi5_cmf_mt_Bmt_p_MH.png)

![](z_bmi5_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi5_cmf_ft_Bft_p_MH.png)

![](z_bmi5_cmf_ft_Bft_p_QQ.png)


### 1.5 y - h vs. cmf (F-test p-value)

![](z_bmi6_cmf_h_p_MH.png)

![](z_bmi6_cmf_h_p_QQ.png)


### 1.5 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi6_h_B1_p_MH.png)

![](z_bmi6_h_B1_p_QQ.png)

- B2

![](z_bmi6_h_B2_p_MH.png)

![](z_bmi6_h_B2_p_QQ.png)

- B3

![](z_bmi6_h_B3_p_MH.png)

![](z_bmi6_h_B3_p_QQ.png)

- B4

![](z_bmi6_h_B4_p_MH.png)

![](z_bmi6_h_B4_p_QQ.png)


### 1.5 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi6_cmf_Bc_p_MH.png)

![](z_bmi6_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi6_cmf_Bm_p_MH.png)

![](z_bmi6_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi6_cmf_Bf_p_MH.png)

![](z_bmi6_cmf_Bf_p_QQ.png)


### 1.5 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi6_cmf_mt_Bmt_p_MH.png)

![](z_bmi6_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi6_cmf_ft_Bft_p_MH.png)

![](z_bmi6_cmf_ft_Bft_p_QQ.png)


### 2 y - h vs. cmf (F-test p-value)

![](z_bmi7_cmf_h_p_MH.png)

![](z_bmi7_cmf_h_p_QQ.png)


### 2 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi7_h_B1_p_MH.png)

![](z_bmi7_h_B1_p_QQ.png)

- B2

![](z_bmi7_h_B2_p_MH.png)

![](z_bmi7_h_B2_p_QQ.png)

- B3

![](z_bmi7_h_B3_p_MH.png)

![](z_bmi7_h_B3_p_QQ.png)

- B4

![](z_bmi7_h_B4_p_MH.png)

![](z_bmi7_h_B4_p_QQ.png)


### 2 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi7_cmf_Bc_p_MH.png)

![](z_bmi7_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi7_cmf_Bm_p_MH.png)

![](z_bmi7_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi7_cmf_Bf_p_MH.png)

![](z_bmi7_cmf_Bf_p_QQ.png)


### 2 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi7_cmf_mt_Bmt_p_MH.png)

![](z_bmi7_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi7_cmf_ft_Bft_p_MH.png)

![](z_bmi7_cmf_ft_Bft_p_QQ.png)


### 3 y - h vs. cmf (F-test p-value)

![](z_bmi8_cmf_h_p_MH.png)

![](z_bmi8_cmf_h_p_QQ.png)


### 3 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi8_h_B1_p_MH.png)

![](z_bmi8_h_B1_p_QQ.png)

- B2

![](z_bmi8_h_B2_p_MH.png)

![](z_bmi8_h_B2_p_QQ.png)

- B3

![](z_bmi8_h_B3_p_MH.png)

![](z_bmi8_h_B3_p_QQ.png)

- B4

![](z_bmi8_h_B4_p_MH.png)

![](z_bmi8_h_B4_p_QQ.png)


### 3 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi8_cmf_Bc_p_MH.png)

![](z_bmi8_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi8_cmf_Bm_p_MH.png)

![](z_bmi8_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi8_cmf_Bf_p_MH.png)

![](z_bmi8_cmf_Bf_p_QQ.png)


### 3 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi8_cmf_mt_Bmt_p_MH.png)

![](z_bmi8_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi8_cmf_ft_Bft_p_MH.png)

![](z_bmi8_cmf_ft_Bft_p_QQ.png)


### 5 y - h vs. cmf (F-test p-value)

![](z_bmi9_cmf_h_p_MH.png)

![](z_bmi9_cmf_h_p_QQ.png)


### 5 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi9_h_B1_p_MH.png)

![](z_bmi9_h_B1_p_QQ.png)

- B2

![](z_bmi9_h_B2_p_MH.png)

![](z_bmi9_h_B2_p_QQ.png)

- B3

![](z_bmi9_h_B3_p_MH.png)

![](z_bmi9_h_B3_p_QQ.png)

- B4

![](z_bmi9_h_B4_p_MH.png)

![](z_bmi9_h_B4_p_QQ.png)


### 5 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi9_cmf_Bc_p_MH.png)

![](z_bmi9_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi9_cmf_Bm_p_MH.png)

![](z_bmi9_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi9_cmf_Bf_p_MH.png)

![](z_bmi9_cmf_Bf_p_QQ.png)


### 5 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi9_cmf_mt_Bmt_p_MH.png)

![](z_bmi9_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi9_cmf_ft_Bft_p_MH.png)

![](z_bmi9_cmf_ft_Bft_p_QQ.png)


### 7 y - h vs. cmf (F-test p-value)

![](z_bmi10_cmf_h_p_MH.png)

![](z_bmi10_cmf_h_p_QQ.png)


### 7 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi10_h_B1_p_MH.png)

![](z_bmi10_h_B1_p_QQ.png)

- B2

![](z_bmi10_h_B2_p_MH.png)

![](z_bmi10_h_B2_p_QQ.png)

- B3

![](z_bmi10_h_B3_p_MH.png)

![](z_bmi10_h_B3_p_QQ.png)

- B4

![](z_bmi10_h_B4_p_MH.png)

![](z_bmi10_h_B4_p_QQ.png)


### 7 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi10_cmf_Bc_p_MH.png)

![](z_bmi10_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi10_cmf_Bm_p_MH.png)

![](z_bmi10_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi10_cmf_Bf_p_MH.png)

![](z_bmi10_cmf_Bf_p_QQ.png)


### 7 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi10_cmf_mt_Bmt_p_MH.png)

![](z_bmi10_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi10_cmf_ft_Bft_p_MH.png)

![](z_bmi10_cmf_ft_Bft_p_QQ.png)


### 8 y - h vs. cmf (F-test p-value)

![](z_bmi11_cmf_h_p_MH.png)

![](z_bmi11_cmf_h_p_QQ.png)


### 8 y - h betas (Prob(|t| > 0))

- B1

![](z_bmi11_h_B1_p_MH.png)

![](z_bmi11_h_B1_p_QQ.png)

- B2

![](z_bmi11_h_B2_p_MH.png)

![](z_bmi11_h_B2_p_QQ.png)

- B3

![](z_bmi11_h_B3_p_MH.png)

![](z_bmi11_h_B3_p_QQ.png)

- B4

![](z_bmi11_h_B4_p_MH.png)

![](z_bmi11_h_B4_p_QQ.png)


### 8 y - cmf betas (Prob(|t| > 0))

- Bc

![](z_bmi11_cmf_Bc_p_MH.png)

![](z_bmi11_cmf_Bc_p_QQ.png)

- Bm

![](z_bmi11_cmf_Bm_p_MH.png)

![](z_bmi11_cmf_Bm_p_QQ.png)

- Bf

![](z_bmi11_cmf_Bf_p_MH.png)

![](z_bmi11_cmf_Bf_p_QQ.png)


### 8 y - cmf_mt and cmf_ft (Prob(|t| > 0))

- Bmt

![](z_bmi11_cmf_mt_Bmt_p_MH.png)

![](z_bmi11_cmf_mt_Bmt_p_QQ.png)

- Bft

![](z_bmi11_cmf_ft_Bft_p_MH.png)

![](z_bmi11_cmf_ft_Bft_p_QQ.png)

