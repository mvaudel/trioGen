##
#
# This script loads results from Triogen.
#
##

# Libraries
lib = NULL
lib = "~/R"

library(scales, lib.loc = lib)
library(backports, lib = lib)
library(vctrs, lib = lib)
library(crayon, lib = lib)
library(tidyr, lib = lib)
library(dplyr, lib = lib)
library(gamlss.data, lib = lib)
library(gamlss.dist, lib = lib)
library(gamlss, lib = lib)
library(withr, lib.loc = lib)
library(labeling, lib.loc = lib)
library(digest, lib.loc = lib)
library(reshape2, lib.loc = lib)
library(ggplot2, lib.loc = lib)
library(grid, lib.loc = lib)
library(scico, lib.loc = lib)
library(gtable, lib.loc = lib)
library(conflicted, lib.loc = lib)

theme_set(theme_bw(base_size = 24))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# pheno names

phenos <- paste0("z_bmi", 0:11)

# variant coordinates

variantDF <- read.table(
    file = "docs/lm_test/7-markerinfo.gz",
    header = F,
    quote = "",
    stringsAsFactors = F
)
names(variantDF) <- c("chrom", "pos", "variantId", "ref", "alt", "typed", "info", "refPanelAF")

# Iterate all phenos

for (pheno in phenos) {
    
    # Load results
    
    print(paste0(pheno, " - Loading"))
    
    phenoDF <- read.table(
        file = paste0("docs/lm_test/chr_7_", pheno, ".gz"),
        header = T,
        quote = "",
        stringsAsFactors = F
    )
    
    phenoDF %>%
        left_join(
            variantDF,
            by = "variantId"
        ) -> phenoDF
    
    
    # Export QQ and MHs
    
    phenoDF %>%
        filter(
            refPanelAF > 0.05 & info >= 0.7
        ) %>%
        arrange(
            pos
        ) -> phenoDF
    
    pColumns <- c("cmf_h_p", "h_B1_p", "h_B2_p", "h_B3_p", "h_B4_p", "cmf_Bc_p", "cmf_Bm_p", "cmf_Bf_p", "cmf_mt_Bmt_p", "cmf_ft_Bft_p")
    colors <- c("black", scico(n = 4, palette = "batlow", end = 0.8), scico(n = 3, palette = "hawaii", end = 0.8), scico(n = 2, palette = "cork", begin = 0.2, end = 0.8))
    labels <- c("cmf vs. h\nF-test p-value [-log10]", "h\nβ1 p-value [-log10]", "h\nβ2 p-value [-log10]", "h\nβ3 p-value [-log10]", "h\nβ4 p-value [-log10]", "cmf\nβc p-value [-log10]", "cmf\nβm p-value [-log10]", "cmf\nβf p-value [-log10]", "cmf_mt\nβmt p-value [-log10]", "cmf_fr\nβft p-value [-log10]")
    
    for (i in 1:length(pColumns)) {
        
        pColumn <- pColumns[i]
        axisLabel <- labels[i]
        
        plotDF <- phenoDF %>%
            filter(
                !is.na(!!sym(pColumn))
            )
        
        
        # MH
        
        print(paste0(pheno, " - MH ", pColumn))
        
        plotDF$logValues <- -log10(plotDF[[pColumn]])
        
        yMax <- ceiling(max(plotDF$logValues))
        
        pValuePlot <- ggplot(
            data = plotDF
        ) +
            geom_point(
                mapping = aes(
                    x = pos,
                    y = logValues
                ),
                col = colors[i],
                size = 2
            ) +
            scale_y_continuous(
                name = axisLabel,
                limits = c(0, yMax),
                expand = expand_scale(
                    mult = c(0, 0.05)
                )
            ) +
            scale_x_continuous(
                name = "Chromosome 7",
                expand = expand_scale(
                    mult = 0.02
                )
            ) +
            theme(
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank()
            )
        
        png(paste0("docs/lm_test/", pheno, "_", pColumn, "_MH.png"), width = 900, height = 600)
        plot(pValuePlot)
        dummy <- dev.off()
        
        
        # QQ
        
        print(paste0(pheno, " - QQ ", pColumn))
        
        plotDF %>% 
            arrange(
                logValues
            ) %>%
            mutate(
                expectedP = sort(-log10(runif(n = nrow(plotDF))))
            ) -> plotDF
        
        yMax = max(yMax, ceiling(max(plotDF$expectedP)))
        
        pValuePlot <- ggplot(
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
                    y = logValues
                ),
                col = colors[i],
                size = 2
            ) +
            scale_y_continuous(
                name = axisLabel,
                limits = c(0, yMax),
                expand = expand_scale(
                    mult = 0.02
                )
            ) +
            scale_x_continuous(
                name = "Expected p-value [-log10]",
                limits = c(0, yMax),
                expand = expand_scale(
                    mult = 0.02
                )
            )
        
        png(paste0("docs/lm_test/", pheno, "_", pColumn, "_QQ.png"), width = 900, height = 600)
        plot(pValuePlot)
        dummy <- dev.off()
        
    }
}
