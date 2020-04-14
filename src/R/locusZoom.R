##
#
# This script loads results from Triogen and plots LocusZoom plots.
# 
# This script excepts the following arguments:
# 1- Triogen results as obtained from the LocusZoom command.
# 2- Target SNP.
# 3- Phenotype.
# 4- Output stem
# 5- Path where the libraries are installed.
#
##

# Command line arguments

# args <- commandArgs(TRUE)

args <- c(
    "docs/tmp/rs287621_zBMI3_locusZoomData.gz",
    "rs287621",
    "zbmi_3",
    "docs/tmp/rs287621_zBMI3_locusZoom",
    NULL
)

resultsFile <- args[1]
variant <- args[2]
pheno <- args[3]
outputStem <- args[4]


# Libraries

lib <- args[5]

library(scales, lib.loc = lib)
library(backports, lib.loc = lib)
library(vctrs, lib.loc = lib)
library(crayon, lib.loc = lib)
library(tidyr, lib.loc = lib)
library(dplyr, lib.loc = lib)
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

#' Returns a ggplot object with the LocusZoom for the given model.
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


# Plot settings

theme_set(theme_bw(base_size = 24))


# Load TrioGen results and recombination rates

resultsDF <- read.table(
    file = resultsFile,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

contigs <- unique(
    resultsDF$contig
)

recombinationRatesDFs <- list()

for (i in 1:length(contigs)) {
    
    contig <- contigs[i]
    
    recombinationRateFilePath <- paste0("resources/recombination_rates/genetic_map_GRCh37_chr", contig, ".txt.gz")
    
}