##
#
# This script extracts the hits with best p-value. 
# 
# Results are expected to be split by chromosome under the naming scheme "p_chrX.pheno.gz"
# 
# SNP coordinates tables are expected to be under the naming scheme "8-markerinfo"
# Marker infos are exected to contain the following columns: "chrom", "pos", "variantId", "ref", "alt", "typed", "info", "refPanelAF"
# 
# This script excepts the following arguments:
# 1- Folder containing Triogen results.
# 2- Path to folder containing SNP coordinates table
# 3- Folder where to write the snp lists
# 4- Path where the libraries are installed
#
##


# Command line arguments

args <- commandArgs(TRUE)

outputFolder <- args[3]


# Libraries

lib = args[4]

library(scales, lib.loc = lib)
library(backports, lib = lib)
library(vctrs, lib = lib)
library(crayon, lib = lib)
library(tidyr, lib = lib)
library(dplyr, lib = lib)
library(withr, lib.loc = lib)
library(labeling, lib.loc = lib)
library(digest, lib.loc = lib)
library(reshape2, lib.loc = lib)
library(ggplot2, lib.loc = lib)
library(grid, lib.loc = lib)
library(scico, lib.loc = lib)
library(gtable, lib.loc = lib)
library(conflicted, lib.loc = lib)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# Parameters

distance <- 10E3
hPvalueThreshold <- 1E-6
cmfPvalueThreshold <- 1E-6
poePvalueThreshold <- 1E-2


# Functions

#' Returns the top hits.
#' 
#' @param associationDF the association data frame
#' @param pColumn the name of the column containing the p-values
#' @param threshold the p-value threshold
#' 
#' @return a ggplot object with the MH
getTopHits <- function(
    associationDF, 
    pColumn, 
    threshold
) {
    
    # Sanity checks
    
    if (! pColumn %in% names(associationDF)) {
        
        stop(paste0(pColumn), " not found in data frame.")
        
    }
    if (! "snp" %in% names(associationDF)) {
        
        stop("snp column not found in data frame.")
        
    }
    if (! "chrom" %in% names(associationDF)) {
        
        stop("chrom column not found in data frame.")
        
    }
    if (! "pos" %in% names(associationDF)) {
        
        stop("pos column not found in data frame.")
        
    }
    if (length(unique(associationDF$chrom)) > 1) {
        
        stop("Only one chromosome at a time should be provided to top-hits pruning.")
        
    }
    
    
    # Filter p-values
    
    pruneDF <- associationDF %>%
        select(
            pos, !!sym(pColumn)
        ) %>%
        rename(
            p = !!sym(pColumn)
        ) %>%
        filter(
            p <= threshold
        ) %>%
        arrange(
            p
        )
    
    
    # Pick best hit and remove neighbors iteratively
    
    results <- c()
    
    while(nrow(pruneDF) > 0) {
        
        results[lenth(results) + 1] <- pruneDF$snp[1]
        
        pos <- pruneDF$pos[1]
        
        pruneDF %>%
            filter(
                pos <= pos - distance | pos >= pos - distance
            ) -> pruneDF
        
    }
    
    return(results)
    
}


# Plot settings

theme_set(theme_bw(base_size = 24))

## MH
mhColor1 <- "grey20"
mhColor2 <- "grey40"


# Chromosome lengths in GRCh37.p13 (hg19) from Ensembl
chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
chromosomeLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 	155270560)
genomeLength <- sum(chromosomeLength)
chromosomeStart <- cumsum(chromosomeLength) - chromosomeLength
chromosomeMiddle <- chromosomeStart + chromosomeLength / 2


# Load TrioGen results and SNP coordinates

resultsPath <- args[1]
coordinatesPath <- args[2]

pDFs <- list()

for (chr in c(1:22, "X")) {
    
    trioGenFile <- file.path(resultsPath, paste0("p_chr_", chr, ".", pheno, ".gz"))
    markerInfoFile <- file.path(coordinatesPath, paste0(chr, "-markerinfo"))
    
    if (file.exists(trioGenFile) && file.exists(markerInfoFile)) {
        
        trioGenDF <- read.table(
            file = trioGenFile,
            header = T,
            quote = "",
            stringsAsFactors = F
        )
        
        markerInfoDF <- read.table(
            file = markerInfoFile,
            header = F,
            quote = "",
            stringsAsFactors = F
        )
        names(markerInfoDF) <- c("chrom", "pos", "variantId", "ref", "alt", "typed", "info", "refPanelAF")
        
        trioGenDF %>%
            left_join(
                markerInfoDF,
                by = "variantId"
            ) %>%
            filter(
                !is.na(info) & info > 0.7
            ) -> trioGenDF
        
        pDFs[[length(pDFs) + 1]] <- trioGenDF
        
    }
}

pDF <- do.call("rbind", pDFs)


# Free memory

pDFs <- NULL
trioGenDF <- NULL
markerInfoDF <- NULL


# Sanity check

if (length(pDF) == 0 || nrow(pDF) == 0) {
    
    stop(paste0("No data for phenotype ", pheno))
    
}
if (! "chrom" %in% names(pDF)) {
    
    stop("Missing crom column")
    
}
if (! "pos" %in% names(pDF)) {
    
    stop("Missing pos column")
    
}


# Order data frame

pDF %>% 
    arrange(
        chrom, pos
    ) -> pDF


# F-test

plot <- getMh(
    associationDF = pDF,
    pColumn = "cmf_h_p"
)

png(file.path(resultsPath, paste0(pheno, "_", "cmf_h_p", "_MH.png")), width = 900, height = 600)
plot(plot)
dummy <- dev.off()

plot <- getQQ(
    associationDF = pDF,
    pColumn = "cmf_h_p"
)

png(file.path(resultsPath, paste0(pheno, "_", "cmf_h_p", "_QQ.png")), width = 900, height = 600)
plot(plot)
dummy <- dev.off()

maxY <- max(
    pDF$h_B1_p,
    pDF$h_B2_p,
    pDF$h_B3_p,
    pDF$h_B4_p,
    pDF$cmf_Bc_p,
    pDF$cmf_Bm_p,
    pDF$cmf_Bf_p,
    pDF$cmf_mt_Bmt_p,
    pDF$cmf_ft_Bft_p
)

for (colname in c("h_B1_p", "h_B2_p", "h_B3_p", "h_B4_p", "cmf_Bc_p", "cmf_Bm_p", "cmf_Bf_p", "cmf_mt_Bmt_p", "cmf_ft_Bft_p")) {
    
    plot <- getMh(
        associationDF = pDF,
        pColumn = colname,
        maxY = maxY
    )
    
    png(file.path(resultsPath, paste0(pheno, "_", colname, "_MH.png")), width = 900, height = 600)
    plot(plot)
    dummy <- dev.off()
    
    plot <- getQQ(
        associationDF = pDF,
        pColumn = colname
    )
    
    png(file.path(resultsPath, paste0(pheno, "_", colname, "_QQ.png")), width = 900, height = 600)
    plot(plot)
    dummy <- dev.off()
    
}


# Write doc

docFile <- file.path(resultsPath, paste0(pheno, ".md"))
write(x = paste0("# ", pheno, "\n\n"), file = docFile, append = F)

write(x = paste0("\n### h vs. cmf (F-test p-value)\n"), file = docFile, append = T)
write(x = paste0("![](", pheno, "_", "cmf_h_p", "_MH.png)\n"), file = docFile, append = T)
write(x = paste0("![](", pheno, "_", "cmf_h_p", "_QQ.png)\n"), file = docFile, append = T)

for (i in 1:4) {
    
    write(x = paste0("\n### h B", i, " (Prob(|t| > 0))\n"), file = docFile, append = T)
    
    write(x = paste0("![](", pheno, "_h_B", i, "_p_MH.png)\n"), file = docFile, append = T)
    write(x = paste0("![](", pheno, "_h_B", i, "_p_QQ.png)\n"), file = docFile, append = T)
    
    
}

for (i in c("c", "m", "f")) {
    
    write(x = paste0("\n### cmf B", i, " (Prob(|t| > 0))\n"), file = docFile, append = T)
    
    write(x = paste0("![](", pheno, "_cmf_B", i, "_p_MH.png)\n"), file = docFile, append = T)
    write(x = paste0("![](", pheno, "_cmf_B", i, "_p_QQ.png)\n"), file = docFile, append = T)
    
    
}

write(x = paste0("\n### cmf_mt Bmt (Prob(|t| > 0))\n"), file = docFile, append = T)

write(x = paste0("![](", pheno, "_cmf_mt_Bmt_p_MH.png)\n"), file = docFile, append = T)
write(x = paste0("![](", pheno, "_cmf_mt_Bmt_p_QQ.png)\n"), file = docFile, append = T)


write(x = paste0("\n### cmf_ft Bft (Prob(|t| > 0))\n"), file = docFile, append = T)

write(x = paste0("![](", pheno, "_cmf_ft_Bft_p_MH.png)\n"), file = docFile, append = T)
write(x = paste0("![](", pheno, "_cmf_ft_Bft_p_QQ.png)\n"), file = docFile, append = T)


