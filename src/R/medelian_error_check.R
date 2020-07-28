##
#
# This script makes plots from the Mendelian error check cli.
#
##

# Libraries

library("tidyr")
library("dplyr")
library("ggplot2")
library("scico")


# Load data

mendelianErrorDF <- read.table(
    file = "tmp/chr22.gz",
    header = T,
    stringsAsFactors = F,
    sep = "\t",
    comment.char = '#'
)

mendelianErrorDF %>% gather(
    key = "error_type",
    value = "p",
    h2_1_p,
    h2_2_p,
    h4_1_p,
    h4_2_p
) -> mendelianErrorDF2

mendelianErrorDF %>% gather(
    key = "error_type",
    value = "p",
    prevalence_before_check,
    prevalence_after_check
) -> mendelianErrorDF3

# Plot density

ggplot() +
    theme_bw() +
    geom_point(
        data = mendelianErrorDF,
        mapping = aes(
            x = h4_2_obs
        )
    )

ggplot() +
    theme_bw() +
    geom_point(
        data = mendelianErrorDF2 %>% filter(p < 2),
        mapping = aes(
            x = prevalence_before_check,
            y = prevalence_after_check,
            col = error_type
        )
    )

ggplot() +
    theme_bw() +
    geom_density(
        data = mendelianErrorDF3 %>% filter(p < 2),
        mapping = aes(
            x = p,
            col = error_type
        )
    )

table(mendelianErrorDF$typed)
