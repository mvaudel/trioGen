
##
#
# This file builds test files to evaluate the transmission command.
# Warning: any change to this script must be implemented in the command test.
# 
##

# Libraries

library(tidyr)
library(dplyr)

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


# Ground truth results

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
    file = "src/main/resources/transmission/ground_truth.txt",
    row.names = F,
    col.names = T,
    quote = F
)


# Export test files

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

