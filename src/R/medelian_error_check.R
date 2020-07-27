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
    file = "tmp/chr22/chr22",
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

# Plot density

ggplot() +
    theme_bw() +
    geom_density(
        data = mendelianErrorDF %>% filter(typed == 0),
        mapping = aes(
            x = p,
            col = error_type
        )
    )

ggplot() +
    theme_bw() +
    geom_density(
        data = mendelianErrorDF2 %>% filter(typed == 0 & p < 2),
        mapping = aes(
            x = p,
            col = error_type
        )
    )

table(mendelianErrorDF$typed)
