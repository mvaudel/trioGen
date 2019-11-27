##
#
# This script loads results from Triogen.
#
##


# Load results from chr 7

print("Loading")

chr7DF <- read.table(
    file = "tmp/chr_7.gz",
    header = T,
    quote = "",
    stringsAsFactors = F
)

for (pheno in unique(chr7DF$phenotype)) {
    
    chr7PhenoDF <- chr7DF %>%
        filter(
            phenotype == pheno
        )
    
    fileName <- paste0("tmp/chr_7", pheno, ".gz")
    
    write.table(
        x = chr7PhenoDF,
        file = fileName,
        sep = "\t",
        row.names = F,
        col.names = T
    )
    
}
