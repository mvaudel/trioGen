##
#
# This script loads results from Triogen and plots LocusZoom plots.
# 
# This script excepts the following arguments:
# 1- Triogen results as obtained from the LocusZoom command.
# 2- Ensembl gene mapping as obtained from the LocusZoom command.
# 3- Target SNP.
# 4- Phenotype.
# 5- Output stem
# 6- Path where the libraries are installed.
#
##

# Command line arguments

# args <- commandArgs(TRUE)

args <- c(
    "docs/lz/data/rs287621_z_bmi3_locusZoomData.gz",
    "docs/lz/data/rs287621_z_bmi3_locusZoomGenes.gz",
    "rs287621",
    "zbmi_3",
    "docs/lz/plots/rs287621_z_bmi3_locusZoom",
    NULL
)

resultsFile <- args[1]
genesFile <- args[2]
variant <- args[3]
pheno <- args[4]
outputStem <- args[5]


# Libraries

lib <- NULL

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

#' Returns a grob with the LocusZoom for the given model.
#' 
#' @param resultsPlotDF The results data frame to plot.
#' @param genesDF The genes data frame to plot.
#' @param recombinationRatesDF The recombination rates data frame to plot.
#' @param minBp The minimum of the bp window.
#' @param maxBp The maximum of the bp window.
#' 
#' @return a ggplot object with the MH
getLocusZoom <- function(
    resultsPlotDF, 
    genesPlotDF,
    recombinationRatesPlotDF,
    minBp,
    maxBp
) {
    
    genes <- unique(genesPlotDF$name)
    biotypes <- levels(genesPlotDF$biotype)
    
    yBreaks <- 5 * (0:(floor(max(resultsPlotDF$logP, -log10(5e-8))/5)))
    yLabels <- as.character(yBreaks)
    yLabels[1] <- ""
    
    
    # Plots
    
    locusPlot <- ggplot() +
        geom_line(
            data = recombinationRatesPlotDF,
            mapping = aes(
                x = position,
                y = scaledY
            ),
            col = "blue3"
        ) + 
        geom_hline(
            yintercept = -log10(5e-8), 
            col = "green4", 
            size = 0.3
        ) +
        geom_point(
            data = resultsPlotDF %>%
                filter(typed == 0),
            mapping = aes(
                x = position,
                y = logP,
                col = ldBinned
            ),
            size = 4
        ) +
        geom_point(
            data = resultsPlotDF %>%
                filter(typed == 1),
            mapping = aes(
                x = position,
                y = logP,
                fill = ldBinned
            ),
            shape = 21,
            col = "black",
            size = 4
        ) +
        geom_point(
            data = targetDF,
            mapping = aes(
                x = position,
                y = logP
            ),
            shape = 18,
            col = "red3",
            size = 4
        ) +
        scale_x_continuous(
            name = paste0("Chromosome ", contig, " [Mbp]"),
            limits = c(minBp, maxBp)
        ) + 
        scale_y_continuous(
            name = "p-value [-log10]",
            limits = c(0, 1.1 * max(resultsPlotDF$logP, -log10(5e-8))),
            expand = expansion(
                mult = 0.02
            ),
            breaks = yBreaks,
            labels = yLabels
        ) +
        scale_color_manual(
            name = expression(r^2),
            values = c(scico(
                n = 4,
                palette = locusPalette,
                direction = -1,
                end = 0.9
            ), "grey"),
            breaks = rev(ldLevels),
            labels = rev(ldLevels),
            drop = F
        ) +
        scale_fill_manual(
            name = expression(r^2),
            values = c(scico(
                n = 4,
                palette = locusPalette,
                direction = -1,
                end = 0.9
            ), "grey"),
            breaks = rev(ldLevels),
            labels = rev(ldLevels),
            drop = F
        ) +
        facet_grid(
            variable ~ .
        ) +
        theme(
            panel.spacing = unit(0, "lines")
        ) + 
        guides(
            fill = guide_legend(
                override.aes = list(size = 8)
            )
        )
    
    locusZoomGrob <- ggplotGrob(locusPlot)
    
    if (nrow(genesPlotDF) > 0) {
        
        genesPlot <- ggplot() +
            geom_rect(
                data = genesPlotDF,
                mapping = aes(
                    xmin = start,
                    xmax = end,
                    ymin = y - 0.2,
                    ymax = y + 0.2,
                    fill = biotype,
                    col = biotype
                ),
                alpha = 0.2
            ) +
            scale_x_continuous(
                limits = c(minBp, maxBp)
            ) +
            scale_y_continuous(
                breaks = 1:length(genes),
                labels = genes,
                limits = c(0, length(genes) + 1)
            ) +
            scale_color_manual(
                name = NULL,
                values = scico(
                    n = length(biotypes),
                    palette = "batlow",
                    direction = -1,
                    begin = 0,
                    end = 0.8
                )
            ) +
            scale_fill_manual(
                name = NULL,
                values = scico(
                    n = length(biotypes),
                    palette = "batlow",
                    direction = -1,
                    begin = 0,
                    end = 0.8
                )
            ) +
            theme(
                axis.text.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank()
            )
        
        genesGrob <- ggplotGrob(genesPlot)
        
        genesGrob <- gtable_add_cols(
            x = genesGrob, 
            widths = unit(1, "null"),
            pos = 7
        )
        
        resultGrob <- rbind(locusZoomGrob[1:(nrow(locusZoomGrob) - 5), ], genesGrob[7, ], locusZoomGrob[(nrow(locusZoomGrob) - 4):nrow(locusZoomGrob), ])
        
        
        nVariables <- length(unique(resultsPlotDF$variable))
        nGenes <- nrow(genesPlotDF)
        
        relativeSize <- nGenes / (nGenes + nVariables * panelRefHeight / geneRefHeight)
        
        resultGrob$heights[14] <- unit(relativeSize, "null")
        
        return(resultGrob)
        
    } else {
        
        return(locusZoomGrob)
        
    }
}


# Plot settings

theme_set(theme_bw(base_size = 24))

ldThreshold <- 0.05
xMargin <- 0.4

panelRefHeight <- 300
geneRefHeight <- 120

locusPalette <- "hawaii"


# Load TrioGen results and recombination rates

resultsDF <- read.table(
    file = resultsFile,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

genesDF <- read.table(
    file = genesFile,
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

contigs <- unique(
    resultsDF$contig
)

if (length(contigs) > 1) {
    stop("Only one contig expected.")
}

contig <- contigs[1]

recombinationRateFilePath <- paste0("resources/recombination_rates/genetic_map_GRCh37_chr", contig, ".txt.gz")

if (file.exists(recombinationRateFilePath)) {
    
    recombinationRatesDF <- read.table(
        file = recombinationRateFilePath,
        header = T,
        sep = "\t",
        stringsAsFactors = F
    )
    
    names(recombinationRatesDF) <- c("contig", "position", "rate", "map")
    
}


# Make one plot per model

models <- unique(resultsDF$model)

for (targetModel in models) {
    
    # Limits in genomics coordinates
    
    maxLdX <- max(resultsDF$position[resultsDF$ld >= ldThreshold])
    minLdX <- min(resultsDF$position[resultsDF$ld >= ldThreshold])
    
    minBp <- minLdX - xMargin * (maxLdX - minLdX)
    maxBp <- maxLdX + xMargin * (maxLdX - minLdX)
    
    ldLevels <- c("[0.0, 0.2]", "[0.2, 0.4]", "[0.4, 0.6]", "[0.6, 0.8]", "[0.8, 1.0]")
    
    resultsPlotDF <- resultsDF %>%
        filter(
            model == targetModel & position >= minBp & position <= maxBp & !is.na(p) & p > 0 & variantId != variant
        ) %>%
        mutate(
            position = position / 1e6,
            logP = -log10(p),
            ldBinned = ifelse(ld < 0.2, "[0.0, 0.2]", "[0.2, 0.4]"),
            ldBinned = ifelse(ld >= 0.4, "[0.4, 0.6]", ldBinned),
            ldBinned = ifelse(ld >= 0.6, "[0.6, 0.8]", ldBinned),
            ldBinned = ifelse(ld >= 0.8, "[0.8, 1.0]", ldBinned),
            ldBinned = factor(ldBinned, ldLevels),
            variable = factor(variable)
        ) %>%
        arrange(
            typed,
            ld
        )
    
    if (nrow(resultsPlotDF) > 0) {
        
        targetDF <- resultsDF %>%
            filter(
                model == targetModel & variantId == variant
            ) %>%
            mutate(
                position = position / 1e6,
                logP = -log10(p),
                variable = factor(variable)
            )
        
        
        recombinationRatesPlotDF <- recombinationRatesDF %>%
            filter(
                position >= minBp & position <= maxBp
            ) %>%
            mutate(
                position = position / 1e6,
                scaledY = rate / 100 * 1.1 * max(resultsPlotDF$logP, -log10(5e-8))
            )
        
        genesPlotDF <- genesDF %>%
            filter(
                start < maxBp & end > minBp
            ) %>%
            mutate(
                start = ifelse(start < minBp, minBp, start),
                end = ifelse(end > maxBp, maxBp, end),
                start = start / 1e6,
                end = end / 1e6,
                biotype = factor(biotype)
            ) %>%
            arrange(
                biotype,
                start,
                end
            ) %>%
            mutate(
                y = row_number()
            )
        
        # Make plot
        
        plotGrob <- getLocusZoom(
            resultsPlotDF = resultsPlotDF, 
            genesPlotDF = genesPlotDF,
            recombinationRatesPlotDF = recombinationRatesPlotDF,
            minBp = minBp / 1e6,
            maxBp = maxBp / 1e6
        )
        
        nVariables <- length(unique(resultsPlotDF$variable))
        nGenes <- nrow(genesPlotDF)
        
        png(
            filename = paste0(outputStem, ".", targetModel, ".png"),
            height = nVariables * panelRefHeight + nGenes * geneRefHeight,
            width = 1000
        )
        grid.draw(plotGrob)
        dummy <- dev.off()
        
    }
}
