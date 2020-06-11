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


# Load results from chr 7 targets

print("Loading")

targetsDF <- read.table(
    file = "docs/lm_test/target/7.lm_target.gz",
    header = T,
    quote = "",
    stringsAsFactors = F
)

targetsDF %>%
    filter(
        startsWith(x = phenotype, prefix = "z_bmi") & variantId == "rs287621"
    ) %>%
    mutate(
        phenotype = factor(phenotype, levels = paste0("z_bmi", 0:11))
    ) -> targetsDF


# Export p-value profiles

yMax = ceiling(max(-log10(targetsDF$cmf_h_p)))

pValuePlot <- ggplot(
    data = targetsDF
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

png("docs/lm_test/lm_targets_7_cmf_h_p.png", width = 900, height = 300)
pValuePlot
dummy <- dev.off()


# Export beta profiles

targetsDF %>%
    select(
        phenotype, h_B1_p, h_B2_p, h_B3_p, h_B4_p
    ) %>%
    gather(
        h_B1_p, h_B2_p, h_B3_p, h_B4_p,
        key = "individual",
        value = "p"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF

levels(betaDF$individual) <- c("h1", "h2", "h3", "h4")

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
        name = "cmf mt\np-value [-log10]",
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

png("docs/lm_test/lm_targets_7_h_p.png", width = 900, height = 300)
pPlot
dummy <- dev.off()

targetsDF %>%
    select(
        phenotype, h_B1, h_B2, h_B3, h_B4
    ) %>%
    gather(
        h_B1, h_B2, h_B3, h_B4,
        key = "individual",
        value = "beta"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF
levels(betaDF$individual) <- c("h1", "h2", "h3", "h4")

targetsDF %>%
    select(
        phenotype, h_B1_se, h_B2_se, h_B3_se, h_B4_se
    ) %>%
    gather(
        h_B1_se, h_B2_se, h_B3_se, h_B4_se,
        key = "individual",
        value = "se"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> seDF
levels(seDF$individual) <- c("h1", "h2", "h3", "h4")

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
        name = "h\nBeta ± se"
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

png("docs/lm_test/lm_targets_7_h_beta.png", width = 900, height = 300)
betaPlot
dummy <- dev.off()


# Export cmf_mt p-values

targetsDF %>%
    select(
        phenotype, cmf_mt.Bc.p, cmf_mt.Bm.p, cmf_mt.Bf.p, cmf_mt.Bmt.p
    ) %>%
    gather(
        cmf_mt.Bc.p, cmf_mt.Bm.p, cmf_mt.Bf.p, cmf_mt.Bmt.p,
        key = "individual",
        value = "p"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF

levels(betaDF$individual) <- c("Child", "Father", "Mother", "Mother T")

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
        name = "cmf mt\np-value [-log10]",
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

png("docs/lm_test/lm_targets_7_cmf_mt_p.png", width = 450, height = 300)
pPlot
dummy <- dev.off()


# Export beta profiles

targetsDF %>%
    select(
        phenotype, cmf_mt.Bc, cmf_mt.Bm, cmf_mt.Bf, cmf_mt.Bmt
    ) %>%
    gather(
        cmf_mt.Bc, cmf_mt.Bm, cmf_mt.Bf, cmf_mt.Bmt,
        key = "individual",
        value = "beta"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> betaDF
levels(betaDF$individual) <- c("Child", "Mother", "Father", "Mother T")

targetsDF %>%
    select(
        phenotype, cmf_mt.Bc.se, cmf_mt.Bm.se, cmf_mt.Bf.se, cmf_mt.Bmt.se
    ) %>%
    gather(
        cmf_mt.Bc.se, cmf_mt.Bm.se, cmf_mt.Bf.se, cmf_mt.Bmt.se,
        key = "individual",
        value = "se"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> seDF
levels(seDF$individual) <- c("Child", "Mother", "Father", "Mother T")

betaSeDF <- merge(betaDF, seDF, by = c("phenotype", "individual"), all = T)

betaPlot <- ggplot(
    data = betaSeDF
) +
    theme_bw(
        base_size = 24
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
        name = "Beta [Z-score]"
    ) +
    scale_x_continuous(
        breaks = 1:12,
        labels = timePoints
    ) +
    scale_color_manual(
        name = "rs287621",
        values = c(
            scico(
                n = 3, 
                  palette = "hawai",
                  begin = 0.1,
                  end = 0.8
                ),
            "black"
        )
    ) +
    scale_fill_manual(
        name = "rs287621",
        values = c(
            scico(
                n = 3, 
                  palette = "hawai",
                  begin = 0.1,
                  end = 0.8
                ),
            "black"
        )
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1,
            vjust = 1
        ),
        panel.grid.minor.x = element_blank()
    )

png("docs/lm_test/poe.png", width = 900, height = 600)
betaPlot
dummy <- dev.off()


# Export cmf_ft p-values

targetsDF %>%
    select(
        phenotype, cmf_ft_Bc_p, cmf_ft_Bm_p, cmf_ft_Bf_p, cmf_ft_Bft_p
    ) %>%
    gather(
        cmf_ft_Bc_p, cmf_ft_Bm_p, cmf_ft_Bf_p, cmf_ft_Bft_p,
        key = "individual",
        value = "p"
    ) %>%
    mutate(
        individual = factor(individual, levels = c("cmf_ft_Bc_p", "cmf_ft_Bm_p", "cmf_ft_Bf_p", "cmf_ft_Bft_p"))
    ) -> betaDF

levels(betaDF$individual) <- c("Child", "Father", "Mother", "Father T")

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
        name = "cmf ft\np-value [-log10]",
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

png("docs/lm_test/lm_targets_7_cmf_ft_p.png", width = 450, height = 300)
pPlot
dummy <- dev.off()


# Export beta profiles

targetsDF %>%
    select(
        phenotype, cmf_ft_Bc, cmf_ft_Bm, cmf_ft_Bf, cmf_ft_Bft
    ) %>%
    gather(
        cmf_ft_Bc, cmf_ft_Bm, cmf_ft_Bf, cmf_ft_Bft,
        key = "individual",
        value = "beta"
    ) %>%
    mutate(
        individual = factor(individual, levels = c("cmf_ft_Bc", "cmf_ft_Bm", "cmf_ft_Bf", "cmf_ft_Bft"))
    ) -> betaDF
levels(betaDF$individual) <- c("Child", "Father", "Mother", "Father T")

targetsDF %>%
    select(
        phenotype, cmf_ft_Bc_se, cmf_ft_Bm_se, cmf_ft_Bf_se, cmf_ft_Bft_se
    ) %>%
    gather(
        cmf_ft_Bc_se, cmf_ft_Bm_se, cmf_ft_Bf_se, cmf_ft_Bft_se,
        key = "individual",
        value = "se"
    ) %>%
    mutate(
        individual = factor(individual)
    ) -> seDF
levels(seDF$individual) <- c("Child", "Father", "Mother", "Father T")

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
        name = "cmf ft\nBeta ± se"
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

png("docs/lm_test/lm_targets_7_cmf_ft_beta.png", width = 450, height = 300)
betaPlot
dummy <- dev.off()
