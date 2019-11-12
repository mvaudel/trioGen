
##
#
# This file builds test files to evaluate the transmission command.
# Warning: any change to this script must be implemented in the command test.
# 
##

# Libraries

library(tidyr)
library(dplyr)
library(ggplot2)

theme_set(theme_bw(base_size = 11))

# Parameters

set.seed(03102019)
separators <- c("/", "|")

possibleGenotypes <- list()

for (i in 1:length(separators)) {
    
    separator <- separators[i]
    
    alleleCombinations <- c()
    
    for (allele1 in 0:1) {
        
        for (allele2 in 0:1) {
            
            genotype <- paste0(allele1, separator, allele2)
            
            alleleCombinations[length(alleleCombinations) + 1] <- genotype
            
        }
    }
    
    possibleGenotypes[[i]] <- alleleCombinations
    
}


# Function

getGenotypeCode <- function(genotype) {
    
    allele1 <- substr(genotype, 1, 1)
    allele2 <- substr(genotype, 3, 3)
    
    return(
        ifelse(
            allele1 == "0" & allele2 == "0",
            0,
            ifelse(
                allele1 == "0" & allele2 == "1" | allele1 == "1" & allele2 == "0",
                1,
                2
            )
        )
    )
}

getAlleleCode <- function(genotype) {
    
    allele1 <- substr(genotype, 1, 1)
    allele2 <- substr(genotype, 3, 3)
    
    return(
        ifelse(
            allele1 == "0" & allele2 == "0",
            0,
            ifelse(
                allele1 == "0" & allele2 == "1",
                1,
                ifelse(
                    allele1 == "1" & allele2 == "0",
                    2,
                    3
                )
            )
        )
    )
}


# Make a data frame with all alleles combinations and expected results

print("Iterating allele combinations")

trioDFs <- list()

for (variantI in 1:100) {
    
    variantId <- paste0("rs", variantI)
    
    separator <- sample(
        x = separators,
        size = 1
    )
    
    separatorI <- (variantI %% length(separators)) + 1
    childGenotypes <- possibleGenotypes[[separatorI]]
    childGenotypes <- sample(childGenotypes, length(childGenotypes))
    
    separatorI <- ((variantI + 1) %% length(separators)) + 1
    motherGenotypes <- possibleGenotypes[[separatorI]]
    motherGenotypes <- sample(motherGenotypes, length(motherGenotypes))
    
    separatorI <- ((variantI + 2) %% length(separators)) + 1
    fatherGenotypes <- possibleGenotypes[[separatorI]]
    fatherGenotypes <- sample(fatherGenotypes, length(fatherGenotypes))
    
    triadI <- 1
    
    for (rep in 1:10) {
        
        for (childGenotype in childGenotypes) {
            
            for (motherGenotype in motherGenotypes) {
                
                for (fatherGenotype in fatherGenotypes) {
                    
                    newLine <- data.frame(
                        variantI = variantI,
                        triadI = triadI,
                        variant = variantId,
                        childId = paste0("CHILD", triadI),
                        motherId = paste0("MOTHER", triadI),
                        fatherId = paste0("FATHER", triadI),
                        childGenotype = childGenotype,
                        motherGenotype = motherGenotype,
                        fatherGenotype = fatherGenotype,
                        stringsAsFactors = F
                    )
                    
                    trioDFs[[length(trioDFs) + 1]] <- newLine
                    
                    triadI <- triadI + 1
                    
                }
            }
        }
    }
    
    if (variantI %% 10 == 0) {
        
        print(paste(variantI, "/ 100"))
        
    }
}


print("Merging")

trioDF <- do.call("rbind", trioDFs)


# Ground truth results for extraction

print("Computing h")

trioDF %>%
    arrange(
        variantI,
        triadI
    ) %>%
    mutate(
        motherGenotypeCode = getGenotypeCode(motherGenotype),
        fatherGenotypeCode = getGenotypeCode(fatherGenotype),
        childAlleleCode = getAlleleCode(childGenotype),
        h1 = ifelse(childAlleleCode == 0 | childAlleleCode == 2, 0, 1),
        h2 = motherGenotypeCode - h1,
        h3 = ifelse(childAlleleCode < 2, 0, 1),
        h4 = fatherGenotypeCode - h3
    ) -> trioDF

write.table(
    x = trioDF,
    file = "src/main/resources/transmission/ground_truth_extraction.txt",
    row.names = F,
    col.names = T,
    quote = F,
    sep = "\t"
)


# Export extraction test files

print("Exporting")

trioLines <- c(
    "child father mother",
    unique(
        paste(
            trioDF$childId, trioDF$fatherId, trioDF$motherId, 
            sep = " "
        )
    )
)
writeLines(trioLines, "src/main/resources/transmission/test_trio")

snpIds <- unique(trioDF$variant)
motherIds <- unique(trioDF$motherId)
fatherIds <- unique(trioDF$fatherId)
childIds <- unique(trioDF$childId)

write(
    "## Transmission test file",
    file = "src/main/resources/transmission/test_transmission.vcf",
    append = F
)
write(
    paste0("## ", Sys.time()),
    file = "src/main/resources/transmission/test_transmission.vcf",
    append = T
)
write(
    "## ",
    file = "src/main/resources/transmission/test_transmission.vcf",
    append = T
)
write(
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", childIds, motherIds, fatherIds), collapse = "\t"),
    file = "src/main/resources/transmission/test_transmission.vcf",
    append = T
)

for (i in 1:length(snpIds)) {
    
    pass <- ifelse(i %% 2 == 0, "PASS", "")
    
    header <- c(
        "1", i, snpIds[i], "A", "B", "OK", pass, "bla"
    )
    
    variantDF <- trioDF %>% 
        filter(
            variant == snpIds[i]
        )
    
    childGenotype <- paste0(
        variantDF$childGenotype,
        ":0.05,0.05:0.1:0.9025,0.095,0.0025"
    )
    
    motherGenotype <- paste0(
        variantDF$motherGenotype,
        ":0.123,0.123:0.12:0.5,0.0956,0.002"
    )
    
    fatherGenotype <- paste0(
        variantDF$fatherGenotype,
        ":0.1,0.12:0.145:0.1234,0.8,0.1234"
    )
    
    write(
        paste(c(header, childGenotype, motherGenotype, fatherGenotype), collapse = "\t"),
        file = "src/main/resources/transmission/test_transmission.vcf",
        append = T
    )
    
    if (i %% 10 == 0) {
        
        print(paste(i, "/ 100"))
        
    }
}


# Ground truth results for linear model

print("Computing phenos")

phenoDF <- trioDF %>%
    filter(
        variantI == 1
    ) %>%
    select(
        childId, h1, h2, h3, h4
    )

lmDF <- data.frame(
    h = character(4*4),
    pheno = character(4*4),
    beta = numeric(4*4),
    se = numeric(4*4),
    t = numeric(4*4),
    p = numeric(4*4),
    stringsAsFactors = F
)

for (i in 1:4) {
    
    hColumn <- paste0("h", i)
    
    noise = rnorm(
        n = nrow(phenoDF),
        mean = 0,
        sd = 0.1
    )
    
    phenoColumn <- paste0("pheno", i)
    
    phenoDF[[phenoColumn]] <- phenoDF[[hColumn]] + noise
    
}

for (i in 1:4) {
    for (j in 1:4) {
        
        hColumn <- paste0("h", i)
        phenoColumn <- paste0("pheno", j)
        
        lmResults <- lm(
            formula = as.formula(paste(phenoColumn, "~", hColumn)),
            data = phenoDF
        )
        
        lmSummary <- summary(lmResults)
        
        plot <- ggplot(
            data = phenoDF
        ) +
            geom_point(
                mapping = aes(
                    x = !!sym(hColumn),
                    y = !!sym(phenoColumn)
                ),
                alpha = 0.2
            ) +
            geom_smooth(
                mapping = aes(
                    x = !!sym(hColumn),
                    y = !!sym(phenoColumn)
                ),
                method = "lm",
                col = "blue3",
                fill = "blue3",
                linetype = "dotted"
            ) +
            geom_abline(
                intercept = lmSummary$coefficients[1, 1],
                slope = lmSummary$coefficients[2, 1],
                col = "red",
                linetype = "dashed"
            )
        
        png(
            filename = paste0("src/main/resources/transmission/", hColumn, "_", phenoColumn, ".png"),
            width = 800,
            height = 600
        )
        plot(plot)
        dummy <- dev.off()
        
        k <- (i-1)*4 + j
        lmDF$h[k] <- hColumn
        lmDF$pheno[k] <- phenoColumn
        lmDF$beta[k] <- lmSummary$coefficients[2, 1]
        lmDF$se[k] <- lmSummary$coefficients[2, 2]
        lmDF$t[k] <- lmSummary$coefficients[2, 3]
        lmDF$p[k] <- lmSummary$coefficients[2, 4]
        
    }
}

phenoDF <- phenoDF %>%
    select(
        childId, starts_with("pheno")
    )

write.table(
    x = phenoDF,
    file = "src/main/resources/transmission/phenos_linear_model.txt",
    row.names = F,
    col.names = T,
    quote = F,
    sep = "\t"
)

write.table(
    x = lmDF,
    file = "src/main/resources/transmission/ground_truth_linear_model.txt",
    row.names = F,
    col.names = T,
    quote = F,
    sep = "\t"
)






