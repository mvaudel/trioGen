
library(ggplot2)
library(scico)


# Plot settings

theme_set(theme_bw(base_size = 24))

cmfColors = scico(
    n = 4,
    palette = "batlow",
    begin = 0.1,
    end = 0.7,
    direction = -1
)

hColors = scico(
    n = 4,
    palette = "cork",
    begin = 0.2,
    end = 0.8,
    direction = -1
)


# h

hBetaDF <- data.frame(
    x = factor(c("h2", "h1", "h3", "h4"), levels = c("h2", "h1", "h3", "h4")),
    y = c(2, 6, 4, 1),
    fill = factor(c("h2", "h1", "h3", "h4"), levels = c("h2", "h1", "h3", "h4")),
    stringsAsFactors = F
)

betaPlot <- ggplot(
    data = hBetaDF
) +
    geom_col(
        mapping = aes(
            x = x,
            y = y,
            fill = fill
        )
    ) +
    geom_col(
        mapping = aes(
            x = x,
            y = y
        ),
        fill = NA,
        col = "black",
        size = 1
    ) +
    scale_y_continuous(
        name = "β per parental allele",
        expand = expand_scale(
            mult = c(0, 0.05)
        )
    ) +
    scale_fill_manual(
        values = hColors
    ) +
    theme(
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
    )

png(
    filename = "docs/illustrations/h_beta.png",
    height = 400,
    width = 300
)
plot(betaPlot)
dummy <- dev.off()


# cmf

cmfBetaDF <- data.frame(
    x = factor(c("h2", "h1", "h3", "h4", "h1", "h3"), levels = c("h2", "h1", "h3", "h4")),
    y = c(2, 2, 1, 1, 3, 3),
    fill = factor(c("m", "m", "f", "f", "c", "c"), levels = c("c", "m", "f")),
    stringsAsFactors = F
)

betaPlot <- ggplot() +
    geom_col(
        data = cmfBetaDF,
        mapping = aes(
            x = x,
            y = y,
            fill = fill
        )
    ) +
    geom_col(
        data = hBetaDF,
        mapping = aes(
            x = x,
            y = y
        ),
        fill = NA,
        col = "black",
        size = 1
    ) +
    scale_y_continuous(
        expand = expand_scale(
            mult = c(0, 0.05)
        )
    ) +
    scale_fill_manual(
        values = cmfColors[c(1, 3, 2)]
    ) +
    theme(
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
    )

png(
    filename = "docs/illustrations/cmf_beta.png",
    height = 400,
    width = 300
)
plot(betaPlot)
dummy <- dev.off()


# cmf_mt

cmfBetaDF <- data.frame(
    x = factor(c("h2", "h1", "h3", "h4", "h1", "h3", "h1"), levels = c("h2", "h1", "h3", "h4")),
    y = c(2, 2, 1, 1, 3, 3, 1),
    fill = factor(c("m", "m", "f", "f", "c", "c", "mt"), levels = c("c", "m", "f", "mt")),
    stringsAsFactors = F
)

betaPlot <- ggplot() +
    geom_col(
        data = cmfBetaDF,
        mapping = aes(
            x = x,
            y = y,
            fill = fill
        )
    ) +
    geom_col(
        data = hBetaDF,
        mapping = aes(
            x = x,
            y = y
        ),
        fill = NA,
        col = "black",
        size = 1
    ) +
    scale_y_continuous(
        expand = expand_scale(
            mult = c(0, 0.05)
        )
    ) +
    scale_fill_manual(
        values = cmfColors[c(1, 3, 2, 4)],
        breaks = c("c", "mt", "m", "f")
    ) +
    theme(
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
    )

png(
    filename = "docs/illustrations/cmf_beta_mt.png",
    height = 400,
    width = 300
)
plot(betaPlot)
dummy <- dev.off()
