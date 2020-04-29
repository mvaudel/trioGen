##
#
# This script compares the cmf and cm models.
# 
##


# Libraries

library("dplyr")
library("ggplot2")
library("scico")

# Plot settings

theme_set(theme_bw(base_size = 24))

timePoints <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y")


# Load data

resultsDF <- read.table(
    file = "docs/tmp/targets/lm_targets_7.gz",
    header = T,
    sep = "\t",
    stringsAsFactors = F
)


# Select BMI

resultsDF %>%
    filter(
        phenotype == "z_bmi3"
    ) -> resultsDF


# Plot h

hPlotDF <- data.frame(
    beta = c(resultsDF$h.B1, resultsDF$h.B2, resultsDF$h.B3, resultsDF$h.B4),
    seMin =  c(resultsDF$h.B1 - resultsDF$h.B1.se, resultsDF$h.B2 - resultsDF$h.B2.se, resultsDF$h.B3 - resultsDF$h.B3.se, resultsDF$h.B4 - resultsDF$h.B4.se),
    seMax =  c(resultsDF$h.B1 + resultsDF$h.B1.se, resultsDF$h.B2 + resultsDF$h.B2.se, resultsDF$h.B3 + resultsDF$h.B3.se, resultsDF$h.B4 + resultsDF$h.B4.se),
    logP = round(-log10(c(resultsDF$h.B1.p, resultsDF$h.B2.p, resultsDF$h.B3.p, resultsDF$h.B4.p)), digits = 2),
    variable = c("B1", "B2", "B3", "B4")
)

hPlot <- ggplot() +
    geom_segment(
        data = hPlotDF,
        mapping = aes(
            x = variable,
            xend = variable,
            y = seMin,
            yend = seMax
        ),
        size = 1,
        alpha = 0.5,
        col = "darkblue"
    ) +
    geom_label(
        data = hPlotDF,
        mapping = aes(
            x = variable,
            y = beta,
            label = logP
        ),
        col = "darkblue"
    ) + 
    scale_y_continuous(
        name = "Beta [Z-score]"
    ) +
    theme(
        axis.title.x = element_blank()
    )

png(
    filename = "docs/tmp/cmf_cm.h.png",
    height = 400,
    width = 300
)
plot(hPlot)
dummy <- dev.off()


# Plot cmf_mt

cmfPlotDF <- data.frame(
    beta = c(resultsDF$cmf_mt.Bc, resultsDF$cmf_mt.Bm, resultsDF$cmf_mt.Bf, resultsDF$cmf_mt.Bmt),
    seMin =  c(resultsDF$cmf_mt.Bc - resultsDF$cmf_mt.Bc.se, resultsDF$cmf_mt.Bm - resultsDF$cmf_mt.Bm.se, resultsDF$cmf_mt.Bf - resultsDF$cmf_mt.Bf.se, resultsDF$cmf_mt.Bmt - resultsDF$cmf_mt.Bmt.se),
    seMax =  c(resultsDF$cmf_mt.Bc + resultsDF$cmf_mt.Bc.se, resultsDF$cmf_mt.Bm + resultsDF$cmf_mt.Bm.se, resultsDF$cmf_mt.Bf + resultsDF$cmf_mt.Bf.se, resultsDF$cmf_mt.Bmt + resultsDF$cmf_mt.Bmt.se),
    logP = round(-log10(c(resultsDF$cmf_mt.Bc.p, resultsDF$cmf_mt.Bm.p, resultsDF$cmf_mt.Bf.p, resultsDF$cmf_mt.Bmt.p)), digits = 2),
    variable = factor(c("Bc", "Bm", "Bf", "Bmt"), levels = c("Bc", "Bm", "Bf", "Bmt"))
)

cmfPlot <- ggplot() +
    geom_segment(
        data = cmfPlotDF,
        mapping = aes(
            x = variable,
            xend = variable,
            y = seMin,
            yend = seMax
        ),
        size = 1,
        alpha = 0.5,
        col = "darkred"
    ) +
    geom_label(
        data = cmfPlotDF,
        mapping = aes(
            x = variable,
            y = beta,
            label = logP
        ),
        col = "darkred"
    ) + 
    scale_y_continuous(
        name = "Beta [Z-score]"
    ) +
    theme(
        axis.title.x = element_blank()
    )

png(
    filename = "docs/tmp/cmf_cm.cmf.png",
    height = 400,
    width = 300
)
plot(cmfPlot)
dummy <- dev.off()


# Plot cm_mt

cmfPlotDF <- data.frame(
    beta = c(resultsDF$cm_mt.Bc, resultsDF$cm_mt.Bm, resultsDF$cm_mt.Bmt),
    seMin =  c(resultsDF$cm_mt.Bc - resultsDF$cm_mt.Bc.se, resultsDF$cm_mt.Bm - resultsDF$cm_mt.Bm.se, resultsDF$cm_mt.Bmt - resultsDF$cm_mt.Bmt.se),
    seMax =  c(resultsDF$cm_mt.Bc + resultsDF$cm_mt.Bc.se, resultsDF$cm_mt.Bm + resultsDF$cm_mt.Bm.se, resultsDF$cm_mt.Bmt + resultsDF$cm_mt.Bmt.se),
    logP = round(-log10(c(resultsDF$cm_mt.Bc.p, resultsDF$cm_mt.Bm.p, resultsDF$cm_mt.Bmt.p)), digits = 2),
    variable = factor(c("Bc", "Bm", "Bmt"), levels = c("Bc", "Bm", "Bmt"))
)

cmPlot <- ggplot() +
    geom_segment(
        data = cmfPlotDF,
        mapping = aes(
            x = variable,
            xend = variable,
            y = seMin,
            yend = seMax
        ),
        size = 1,
        alpha = 0.5,
        col = "darkgreen"
    ) +
    geom_label(
        data = cmfPlotDF,
        mapping = aes(
            x = variable,
            y = beta,
            label = logP
        ),
        col = "darkgreen"
    ) + 
    scale_y_continuous(
        name = "Beta [Z-score]"
    ) +
    theme(
        axis.title.x = element_blank()
    )

png(
    filename = "docs/tmp/cmf_cm.cm.png",
    height = 400,
    width = 300
)
plot(cmPlot)
dummy <- dev.off()
