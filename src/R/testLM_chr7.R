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

theme_set(theme_bw(base_size = 13))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# Parameters

timePoints <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y")


# Load results from chr 3 targets

print("Loading")

targets3DF <- read.table(
    file = "docs/lm_test/lm_targets_3",
    header = T,
    quote = "",
    stringsAsFactors = F
)

targets3DF %>%
    filter(
        startsWith(x = phenotype, prefix = "z_bmi") & variantId == "rs11708067"
    ) %>%
    mutate(
        phenotype = factor(phenotype, levels = paste0("z_bmi", 0:11))
    ) -> targets3DF


# Export p-value profiles

yMax = ceiling(max(-log10(targets3DF$cmf_h_p)))

pValuePlot <- ggplot(
    data = targets3DF
) +
    geom_point(
        mapping = aes(
            x = phenotype,
            y = -log10(cmf_h_p)
        ),
        size = 2,
        alpha = 0.8
    ) +
    geom_line(
        mapping = aes(
            x = phenotype,
            y = -log10(cmf_h_p),
            group = 1
        ),
        size = 1.2,
        alpha = 0.5
    ) +
    scale_y_continuous(
        name = "cmf vs. h\nF-test p-value [-log10]",
        limits = c(0, yMax)
    ) +
    scale_x_discrete(
        labels = timePoints
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1,
            vjust = 1
            )
    )

png("docs/lm_test/lm_targets_3_cmf_h_p.png", width = 900, height = 300)
pValuePlot
dummy <- dev.off()


# Export cmf p-values

targets3DF %>%
    select(
        phenotype, cmf_Bc_p, cmf_Bm_p, cmf_Bf_p
    ) %>%
    gather(
        cmf_Bc_p, cmf_Bm_p, cmf_Bf_p,
        key = "individual",
        value = "p"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF

levels(betaDF$individual) <- c("Child", "Father", "Mother")

yMax = ceiling(max(-log10(betaDF$p)))

pPlot <- ggplot(
    data = betaDF
) +
    geom_point(
        mapping = aes(
            x = phenotype,
            y = -log10(p),
            col = individual
        ),
        size = 2,
        alpha = 0.8
    ) +
    geom_line(
        mapping = aes(
            x = phenotype,
            y = -log10(p),
            group = individual,
            col = individual
        ),
        size = 1,
        alpha = 0.5
    ) +
    scale_y_continuous(
        name = "cmf\np-value [-log10]",
        limits = c(0, yMax)
    ) +
    scale_x_discrete(
        labels = timePoints
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1,
            vjust = 1
        ),
        legend.title = element_blank()
    )

png("docs/lm_test/lm_targets_3_cmf_p.png", width = 900, height = 300)
pPlot
dummy <- dev.off()


# Export beta profiles

targets3DF %>%
    select(
        phenotype, cmf_Bc, cmf_Bm, cmf_Bf
    ) %>%
    gather(
        cmf_Bc, cmf_Bm, cmf_Bf,
        key = "individual",
        value = "beta"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF
levels(betaDF$individual) <- c("Child", "Father", "Mother")

targets3DF %>%
    select(
        phenotype, cmf_Bc_se, cmf_Bm_se, cmf_Bf_se
    ) %>%
    gather(
        cmf_Bc_se, cmf_Bm_se, cmf_Bf_se,
        key = "individual",
        value = "se"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> seDF
levels(seDF$individual) <- c("Child", "Father", "Mother")

betaSeDF <- merge(betaDF, seDF, by = c("phenotype", "individual"), all = T)

betaPlot <- ggplot(
    data = betaSeDF
) +
    geom_ribbon(
        mapping = aes(
            x = as.numeric(phenotype),
            ymin = beta - se,
            ymax = beta + se,
            fill = individual
        ),
        alpha = 0.2
    ) +
    geom_point(
        mapping = aes(
            x = as.numeric(phenotype),
            y = beta,
            col = individual
        ),
        size = 2,
        alpha = 0.8
    ) +
    geom_line(
        mapping = aes(
            x = as.numeric(phenotype),
            y = beta,
            group = individual,
            col = individual
        ),
        size = 1,
        alpha = 0.5
    ) +
    scale_y_continuous(
        name = "cmf\nBeta ± se"
    ) +
    scale_x_continuous(
        breaks = 1:12,
        labels = timePoints
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1,
            vjust = 1
        ),
        legend.title = element_blank()
    )

png("docs/lm_test/lm_targets_3_cmf_beta.png", width = 900, height = 300)
betaPlot
dummy <- dev.off()
