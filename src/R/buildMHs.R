##
#
# This script loads results from Triogen. 
# 
# Results are expected to be split by chromosome under the naming scheme "p_chrX.pheno.gz"
# 
# SNP coordinates tables are expected to be under the naming scheme "8-markerinfo"
# Marker infos are exected to contain the following columns: "chrom", "pos", "variantId", "ref", "alt", "typed", "info", "refPanelAF"
# 
# This script excepts the following arguments:
# 1- Folder containing Triogen results.
# 2- Path to folder containing SNP coordinates table
# 3- Phenotype
# 4- folder where to write the plots
# 5- Path where the libraries are installed
#
##


# Command line arguments

args <- commandArgs(TRUE)

pheno <- args[3]
outputFolder <- args[4]


# Libraries

lib = args[5]

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


# Functions

#' Returns a ggplot object with the MH for the given association data frame.
#' 
#' @param associationDF the association data frame
#' @param pColumn the name of the column containing the p-values
#' @param maxY the max value to use for the y axis
#' 
#' @return a ggplot object with the MH
getMh <- function(
    associationDF, 
    pColumn, 
    maxY = NA
) {
    
    # Sanity checks
    
    if (! pColumn %in% names(associationDF)) {
        
        stop(paste0(pColumn), " not found in data frame")
        
    }
    if (! "chrom" %in% names(associationDF)) {
        
        stop("chrom column not found in data frame")
        
    }
    if (! "pos" %in% names(associationDF)) {
        
        stop("pos column not found in data frame")
        
    }
    
    
    # make data frame with only the data needed to plot
    
    plotDF <- associationDF %>%
        select(
            chrom, pos, !!sym(pColumn)
        ) %>%
        rename(
            p = !!sym(pColumn)
        ) %>%
        filter(
            p > 0
        ) %>%
        mutate(
            logP = -log10(p),
            x = chromosomeStart[chrom] + pos,
            color = factor(chrom %% 2, levels = c(0, 1))
        )
    
    
    # Chromosome labels
    
    xLabels <- 1:22
    xLabels[xLabels %% 2 == 0 & xLabels > 17] <- ""
    xLabels <- c(xLabels, "X")
    
    # y axis
    
    if (is.na(maxY)) {
        
        maxY <- 10 * ceiling(max(plotDF$logP / 10))
        
    }
    
    maxY <- max(maxY, 10)
    
    yBreaks <- c(0, 5, -log10(5e-8))
    yLabels <- c("", 5, round(-log10(5e-8), digits = 1))
    
    lastBreak <- floor(max(maxY / 10))
    
    if (lastBreak > 0) {
        
        newBreaks <- 10*(1:lastBreak)
        
        while(length(newBreaks) > 3) {
            
            newBreaks <- newBreaks[c(T, F)]
            
        }
        
        yBreaks <- c(yBreaks, newBreaks)
        yLabels <- c(yLabels, round(newBreaks, digits = 1))
        
    }
    
    
    # Make plot
    mhPlot <- ggplot(
        data = plotDF
    ) + 
        geom_hline(
            yintercept = -log10(5e-8), 
            col = "green4", 
            size = 0.3
        ) + 
        geom_point(
            aes(x = x, 
                y = logP, 
                col = color
            ), 
            size = 2
        ) + 
        scale_y_continuous(
            name = paste0(pColumn, " [-log10]"), 
            breaks = yBreaks, 
            labels = yLabels, 
            expand = expand_scale(
                mult = c(0, 0.05)
            ), 
            limits = c(0, maxY)
        ) + 
        scale_x_continuous(
            name = "Chromosome", 
            breaks = chromosomeMiddle, 
            labels = xLabels, 
            limits = c(0, genomeLength), 
            expand = expand_scale(
                mult = 0.01
            )
        ) + 
        scale_color_manual(
            values = c(mhColor1, mhColor2)
        ) + 
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(size = 0.3),
            strip.background = element_rect(
                fill = "grey99"
            )
        )
    
    return(mhPlot)
    
}

#' Returns a ggplot object with the QQ for the given association data frame.
#' 
#' @param associationDF the association data frame
#' @param pColumn the name of the column containing the p-values
#' @param maxAxis the max value to use for the axes
#' 
#' @return a ggplot object with the MH
getQQ <- function(
    associationDF, 
    pColumn
) {
    
    # Sanity checks
    
    if (! pColumn %in% names(associationDF)) {
        
        stop(paste0(pColumn), " not found in data frame")
        
    }
    if (! "chrom" %in% names(associationDF)) {
        
        stop("chrom column not found in data frame")
        
    }
    if (! "pos" %in% names(associationDF)) {
        
        stop("pos column not found in data frame")
        
    }
    
    
    # make data frame with only the data needed to plot
    
    plotDF <- associationDF %>%
        select(
            chrom, pos, !!sym(pColumn)
        ) %>%
        rename(
            p = !!sym(pColumn)
        ) %>%
        filter(
            p > 0
        ) %>%
        mutate(
            logP = -log10(p)
        ) %>%
        arrange(
            p
        ) %>%
        mutate(
            expectedP = sort(-log10(runif(n = n()))),
            x = chromosomeStart[chrom] + pos,
            color = factor(chrom %% 2, levels = c(0, 1))
        )
    
    
    # Axis
    
    maxValue <- max(plotDF$logP, plotDF$expectedP)
    
    
    # Make plot
    qqPlot <- ggplot(
        data = plotDF
    ) + 
        geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dotted"
        ) +
        geom_point(
            mapping = aes(
                x = expectedP,
                y = logP
            ),
            size = 2
        ) +
        scale_y_continuous(
            name = paste0(pColumn, " [-log10]"), 
            limits = c(0, maxValue),
            expand = expand_scale(
                mult = 0.02
            )
        ) +
        scale_x_continuous(
            name = "Expected p-value [-log10]",
            limits = c(0, maxValue),
            expand = expand_scale(
                mult = 0.02
            )
        )
    
    return(qqPlot)
    
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

if (nrow(pDF) == 0) {
    
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

write(x = "TrioGen v.0.3.0-beta on 27,451 full trios, ADHD cases, related, and ethnic outliers excluded. 10 PCs and genotyping batch as covariates. Reference panel AF > 5%, info score >= 0.7, singularities excluded\n", file = docFile, append = T)
write(x = "z_BMI “à la Chris”, see details [here](../pheno/plots.md):\n", file = docFile, append = T)
write(x = "- not the same as the pheno tables used so far!\n", file = docFile, append = T)
write(x = "- BMI adjusted for pregnancy duration!\n\n", file = docFile, append = T)

write(x = "4 Models:\n", file = docFile, append = T)
write(x = "- h: `y = β1 h1 + β2 h2 + β3 h3 + β4 h4  + ε`\n", file = docFile, append = T)
write(x = "- cmf: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + ε`\n", file = docFile, append = T)
write(x = "- cmf_mt: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + βmt h1 + ε`\n", file = docFile, append = T)
write(x = "- cmf_ft: `y = βm (h1 + h2) + βc (h1 + h3) + βf (h3 + h4) + βft h1 + ε`\n", file = docFile, append = T)


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


