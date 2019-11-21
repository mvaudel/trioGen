##
#
# This script creates pheno files to test TrioGen.
#
##

# Libraries

library(scales, lib.loc = "~/R")
library(backports, lib = "~/R")
library(vctrs, lib = "~/R")
library(crayon, lib = "~/R")
library(dplyr, lib = "~/R")
library(gamlss.data, lib = "~/R")
library(gamlss.dist, lib = "~/R")
library(gamlss, lib = "~/R")
library(withr, lib.loc = "~/R")
library(labeling, lib.loc = "~/R")
library(digest, lib.loc = "~/R")
library(reshape2, lib.loc = "~/R")
library(ggplot2, lib.loc = "~/R")
library(grid, lib.loc = "~/R")
library(scico, lib.loc = "~/R")
library(gtable, lib.loc = "~/R")
library(conflicted, lib.loc = "~/R")

theme_set(theme_bw(base_size = 13))

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# Functions

#' Builds a plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
plotPhenos <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY
) {
    
    if (!pheno1 %in% names(df)) {
        stop(paste0("Pheno ", pheno1, " not found in phenoDF: ", paste(names(df))))
    }
    if (!pheno2 %in% names(df)) {
        stop(paste0("Pheno ", pheno2, " not found in phenoDF: ", paste(names(df))))
    }
    
    plotDF <- data.frame(
        x = df[[pheno1]],
        y = df[[pheno2]],
        stringsAsFactors = F
    ) %>%
        filter(
            !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
        ) %>%
        arrange(
            rev(abs(y - x))
        )
    
    # Build the scatter plot
    
    scatterPlot <- ggplot(
        data = plotDF
    ) +
        geom_point(
            mapping = aes(
                x = x,
                y = y
            ),
            col = "black",
            alpha = 0.1
        ) + 
        geom_density_2d(
            mapping = aes(
                x = x,
                y = y
            ),
            col = "white",
        ) +
        scale_x_continuous(
            name = labelX
        ) +
        scale_y_continuous(
            name = labelY
        ) +
        theme(
            legend.position = "none"
        )
    
    
    # Build the density plots
    
    xDensityPlot <- ggplot(
        data = plotDF
    ) + theme_minimal() + 
        geom_density(
            mapping = aes(
                x = x
            ),
            fill = "black",
            alpha = 0.1
        ) +
        scale_x_continuous(
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            expand = c(0, 0)
        ) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank()
        )
    
    yDensityPlot <- ggplot(
        data = plotDF
    ) + theme_minimal() + 
        geom_density(
            mapping = aes(
                x = y
            ),
            fill = "black",
            alpha = 0.1
        ) +
        scale_x_continuous(
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            expand = c(0, 0)
        ) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank()
        ) + 
        coord_flip()
    
    
    # Make grobs from plots
    
    scatterGrob <- ggplotGrob(scatterPlot)
    xDensityGrob <- ggplotGrob(xDensityPlot)
    yDensityGrob <- ggplotGrob(yDensityPlot)
    
    
    # Insert the densities as new row and column in the scatter grob
    
    mergedGrob <- rbind(scatterGrob[1:6, ], xDensityGrob[7, ], scatterGrob[7:nrow(scatterGrob), ], size = "last")
    mergedGrob$heights[7] <- unit(0.15, "null")
    
    yDensityGrob <- gtable_add_rows(
        x = yDensityGrob, 
        heights = unit(rep(0, nrow(mergedGrob) - nrow(yDensityGrob)), "null"), 
        pos = 0
    )
    
    mergedGrob <- cbind(mergedGrob[, 1:5], yDensityGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
    mergedGrob$widths[6] <- unit(0.15, "null")
    
    
    # Plot
    
    png(
        filename = file.path(docsFolder, paste0(pheno1, "-", pheno2, ".png")),
        width = 800,
        height = 600
    )
    grid.draw(mergedGrob)
    dummy <- dev.off()
    
}


# Parameters

timePoints <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y")

theme_set(theme_bw(base_size = 13))


# Paths

phenoFolder <- "/mnt/archive/moba/pheno/v10/V10_1.0.0-190506"
adhdCasesFile <- "/mnt/archive/TED/ted-aux/connections/PDB1382_toDECODE_Iceland_adhd_cases.csv"
adhdBridgeFile <- "/mnt/archive/TED/ted-aux/connections/321-FinalDelJan2017-5410PNs-Sample_Map.txt"
pcaFile <- "/mnt/archive/MOBAGENETICS/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc"

docsFolder <- "docs/pheno"
docsFile <- file.path(docsFolder, "plots.md")


# Load data

print(paste0(Sys.time(), "    Loading data"))

idDF <- read.table(
    file = file.path(phenoFolder, "id.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>%
    filter(
        !is.na(child_SentrixID)
    ) %>%
    mutate(
        sex_number = as.numeric(factor(sex, levels = c("Boy", "Girl"))),
        child_genotyping_batch = ifelse(!is.na(child_SentrixID) & !is.na(child_Harvest_SentrixID) & child_SentrixID == child_Harvest_SentrixID, "Harvest",
                                        ifelse(!is.na(child_SentrixID) & !is.na(child_Rotterdam1_SentrixID) & child_SentrixID == child_Rotterdam1_SentrixID, "Rotterdam1",
                                               ifelse(!is.na(child_SentrixID) & !is.na(child_Rotterdam2_SentrixID) & child_SentrixID == child_Rotterdam2_SentrixID, "Rotterdam2",
                                                      ifelse(!is.na(child_SentrixID) & !is.na(child_NormentMay16_SentrixID) & child_SentrixID == child_NormentMay16_SentrixID, "NormentMay16",
                                                             ifelse(!is.na(child_SentrixID) & !is.na(child_NormentFeb18_SentrixID) & child_SentrixID == child_NormentFeb18_SentrixID, "NormentFeb18",
                                                                    ifelse(!is.na(child_SentrixID) & !is.na(child_Ted_SentrixID) & child_SentrixID == child_Ted_SentrixID, "Ted",
                                                                           NA
                                                                    )  
                                                             )    
                                                      )    
                                               )   
                                        )
        ),
        child_genotyping_batch_number = as.numeric(
            factor(
                child_genotyping_batch
            )
        )
    )

parentDF <- read.table(
    file = file.path(phenoFolder, "parents.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>% 
    select(
        child_SentrixID, 
        mother_age, 
        father_age, 
        father_height,
        father_weight,
        mother_height
    ) %>% mutate(
        father_bmi = ifelse(!is.na(father_weight) & !is.na(father_height) & father_height > 0, 10000 * father_weight / (father_height)^2, NA)
    )

lwDF <- read.table(
    file = file.path(phenoFolder, "length_weight_bmi.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>%
    filter(
        !is.na(child_SentrixID)
    ) %>%
    select(
        child_SentrixID,
        starts_with("length"),
        starts_with("weight")
    )

pregnancyDF <- read.table(
    file = file.path(phenoFolder, "pregnancy.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>%
    select(
        child_SentrixID,
        pregnancy_duration
    )

deliveryDF <- read.table(
    file = file.path(phenoFolder, "delivery.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>%
    select(
        child_SentrixID,
        placenta_weight,
        umbilical_chord_length
    )

childNutritionDF <- read.table(
    file = file.path(phenoFolder, "child_nutrition.gz"),
    header = T,
    sep = "\t",
    stringsAsFactors = F
) %>%
    filter(
        !is.na(child_SentrixID)
    ) %>% 
    select(
        child_SentrixID,
        breastmilk_first_week,
        sugarwater_first_week,
        formula_first_week,
        breastmilk_0m,
        breastmilk_1m,
        breastmilk_2m,
        breastmilk_3m,
        breastmilk_4m,
        breastmilk_5m,
        breastmilk_6m,
        breastmilk_6_8m,
        breastmilk_9_11m,
        breastmilk_12_14m,
        breastmilk_15_18m,
        formula_freq_6m,
        breastmilk_freq_18m
    ) %>%
    mutate(
        sugarwater_first_week = ifelse(!is.na(sugarwater_first_week), 1, 0),
        formula_first_week = ifelse(!is.na(formula_first_week), 1, 0),
        breastmilk_0m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_0m), 1, 0), NA),
        breastmilk_1m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_1m), 1, 0), NA),
        breastmilk_2m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_2m), 1, 0), NA),
        breastmilk_3m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_3m), 1, 0), NA),
        breastmilk_4m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_4m), 1, 0), NA),
        breastmilk_5m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_5m), 1, 0), NA),
        breastmilk_6m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_6m), 1, 0), NA),
        breastmilk_6_8m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_6_8m), 1, 0), NA),
        breastmilk_9_11m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_9_11m), 1, 0), NA),
        breastmilk_12_14m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_12_14m), 1, 0), NA),
        breastmilk_15_18m = ifelse(!is.na(breastmilk_first_week), ifelse(!is.na(breastmilk_15_18m), 1, 0), NA),
        breastmilk_duration = ifelse(!is.na(breastmilk_first_week), ifelse(breastmilk_15_18m == 1, 16.5,
                                                                           ifelse(breastmilk_12_14m == 1, 13,
                                                                                  ifelse(breastmilk_9_11m == 1, 10,
                                                                                         ifelse(breastmilk_6_8m == 1, 7,
                                                                                                ifelse(breastmilk_5m == 1, 6, 
                                                                                                       NA))))), NA),
        formula_freq_6m = factor(formula_freq_6m, levels = c("Never / seldom", "1-3 times a week", "4-6 times a week", "At least once a day"))
    ) %>% 
    select(
        child_SentrixID,
        breastmilk_duration,
        formula_freq_6m
    )

levels(childNutritionDF$formula_freq_6m) <- c(0, 2, 5, 7)
childNutritionDF$formula_freq_6m <- as.numeric(childNutritionDF$formula_freq_6m)


# Filter ethnic outliers and ADHD cases

print(paste0(Sys.time(), "    Excluding cases"))

adhdStatusDF <- read.table(
    file = adhdCasesFile, 
    header = T, 
    stringsAsFactors = F,
    sep = ","
)

adhdSamplemap <- read.table(
    file = adhdBridgeFile, 
    header = T, 
    stringsAsFactors = F, 
    sep='\t'
)

adhdDF <- merge(
    x = adhdSamplemap, 
    y = adhdStatusDF, 
    by.x="ID", 
    by.y="RetrievalDetail_ID"
)

nStart <- nrow(idDF)

idDF %>%
    filter(
        !child_Ted_SentrixID %in% adhdDF$SentrixPosition
    ) -> idDF

nEnd <- nrow(idDF)

print(paste0("Number of ADHD cases excluded: ", nStart - nEnd))

nStart <- nrow(idDF)

idDF %>%
    filter(
        child_core == 1
    ) -> idDF

nEnd <- nrow(idDF)

print(paste0("Number of ethnic outliers excluded: ", nStart - nEnd))
print(paste0("Number of children in pheno file: ", nrow(idDF)))


# Merge

print(paste0(Sys.time(), "    Merging"))

phenoDF <- idDF %>% 
    select(
        child_SentrixID, sex_number, child_genotyping_batch_number
    ) %>% 
    left_join(
        y = childNutritionDF,
        by = "child_SentrixID"
    ) %>%
    left_join(
        y = lwDF,
        by = "child_SentrixID"
    ) %>%
    left_join(
        y = pregnancyDF,
        by = "child_SentrixID"
    ) %>%
    left_join(
        y = deliveryDF,
        by = "child_SentrixID"
    ) %>%
    left_join(
        y = parentDF,
        by = "child_SentrixID"
    )


# Compute child BMI

print(paste0(Sys.time(), "    Computing child BMI"))

for (ageI in 0:11) {
    
    lengthColumn <- paste0("length", ageI)
    weightColumn <- paste0("weight", ageI)
    bmiColumn <- paste0("bmi", ageI)
    
    phenoDF[[bmiColumn]] <- ifelse(!is.na(phenoDF[[lengthColumn]]) & !is.na(phenoDF[[weightColumn]]), 10000 * phenoDF[[weightColumn]] / (phenoDF[[lengthColumn]] * phenoDF[[lengthColumn]]), NA)
    
}

phenoDF %>%
    select(
        -starts_with("length"), -starts_with("weight")
    ) -> phenoDF


# Standardize umbilical chord length

print(paste0(Sys.time(), "    Standardizing umbilical chord length"))

maleDF <- phenoDF %>%
    filter(
        sex_number == 1
    ) %>%
    select(
        child_SentrixID, pregnancy_duration, umbilical_chord_length
    ) %>%
    filter(
        !is.na(pregnancy_duration) & !is.na(umbilical_chord_length)
    )
femaleDF <- phenoDF %>%
    filter(
        sex_number == 2
    ) %>%
    select(
        child_SentrixID, pregnancy_duration, umbilical_chord_length
    ) %>%
    filter(
        !is.na(pregnancy_duration) & !is.na(umbilical_chord_length)
    )

maleModel <- gamlss(
    formula = umbilical_chord_length ~ fp(pregnancy_duration),
    sigma.formula = ~fp(pregnancy_duration),
    family = BCT,
    data = maleDF
)
femaleModel <- gamlss(
    formula = umbilical_chord_length ~ fp(pregnancy_duration),
    sigma.formula = ~fp(pregnancy_duration),
    family = BCT,
    data = femaleDF
)

maleDF$z_umbilical_chord_length = centiles.pred(
    obj = maleModel, 
    xname = "pregnancy_duration", 
    xvalues = maleDF$pregnancy_duration, 
    yval = maleDF$umbilical_chord_length, 
    type = "z-scores"
)

femaleDF$z_umbilical_chord_length = centiles.pred(
    obj = femaleModel, 
    xname = "pregnancy_duration", 
    xvalues = femaleDF$pregnancy_duration, 
    yval = femaleDF$umbilical_chord_length, 
    type = "z-scores"
)

maleDF %>% 
    select(
        child_SentrixID, z_umbilical_chord_length
    ) -> maleDF
femaleDF %>% 
    select(
        child_SentrixID, z_umbilical_chord_length
    ) -> femaleDF
zDF <- rbind(maleDF, femaleDF)

phenoDF %>%
    left_join(
        zDF,
        by = "child_SentrixID"
    ) -> phenoDF


# Standardize placenta weight

print(paste0(Sys.time(), "    Standardizing placenta weight"))

maleDF <- phenoDF %>%
    filter(
        sex_number == 1
    ) %>%
    select(
        child_SentrixID, pregnancy_duration, placenta_weight
    ) %>%
    filter(
        !is.na(pregnancy_duration) & !is.na(placenta_weight)
    )
femaleDF <- phenoDF %>%
    filter(
        sex_number == 2
    ) %>%
    select(
        child_SentrixID, pregnancy_duration, placenta_weight
    ) %>%
    filter(
        !is.na(pregnancy_duration) & !is.na(placenta_weight)
    )

maleModel <- gamlss(
    formula = placenta_weight ~ fp(pregnancy_duration),
    sigma.formula = ~fp(pregnancy_duration),
    family = BCT,
    data = maleDF
)
femaleModel <- gamlss(
    formula = placenta_weight ~ fp(pregnancy_duration),
    sigma.formula = ~fp(pregnancy_duration),
    family = BCT,
    data = femaleDF
)

maleDF$z_placenta_weight = centiles.pred(
    obj = maleModel, 
    xname = "pregnancy_duration", 
    xvalues = maleDF$pregnancy_duration, 
    yval = maleDF$placenta_weight, 
    type = "z-scores"
)

femaleDF$z_placenta_weight = centiles.pred(
    obj = femaleModel, 
    xname = "pregnancy_duration", 
    xvalues = femaleDF$pregnancy_duration, 
    yval = femaleDF$placenta_weight, 
    type = "z-scores"
)

maleDF %>% 
    select(
        child_SentrixID, z_placenta_weight
    ) -> maleDF
femaleDF %>% 
    select(
        child_SentrixID, z_placenta_weight
    ) -> femaleDF
zDF <- rbind(maleDF, femaleDF)

phenoDF %>%
    left_join(
        zDF,
        by = "child_SentrixID"
    ) -> phenoDF


# Standardize BMI

for (ageI in 0:11) {
    
    bmiColumn <- paste0("bmi", ageI)
    zBmiColumn <- paste0("z_bmi", ageI)
    
    print(paste0(Sys.time(), "    Standardizing ", bmiColumn))
    
    maleDF <- phenoDF %>%
        filter(
            sex_number == 1
        ) %>%
        select(
            child_SentrixID, pregnancy_duration, !!sym(bmiColumn)
        ) %>%
        filter(
            !is.na(pregnancy_duration) & !is.na(!!sym(bmiColumn))
        )
    femaleDF <- phenoDF %>%
        filter(
            sex_number == 2
        ) %>%
        select(
            child_SentrixID, pregnancy_duration, !!sym(bmiColumn)
        ) %>%
        filter(
            !is.na(pregnancy_duration) & !is.na(!!sym(bmiColumn))
        )
    
    if (ageI < 3) {
        sigmaFormula <- " ~ fp(pregnancy_duration)"
    } else {
        sigmaFormula <- " ~ pregnancy_duration"
    }
        formula <- paste0(bmiColumn, sigmaFormula)
    
    maleModel <- gamlss(
        formula = as.formula(formula),
        sigma.formula = as.formula(sigmaFormula),
        family = LOGNO,
        data = maleDF
    )
    femaleModel <- gamlss(
        formula = as.formula(formula),
        sigma.formula = as.formula(sigmaFormula),
        family = LOGNO,
        data = femaleDF
    )
    
    maleDF[[zBmiColumn]] = centiles.pred(
        obj = maleModel, 
        xname = "pregnancy_duration", 
        xvalues = maleDF$pregnancy_duration, 
        yval = maleDF[[bmiColumn]], 
        type = "z-scores"
    )
    
    femaleDF[[zBmiColumn]] = centiles.pred(
        obj = femaleModel, 
        xname = "pregnancy_duration", 
        xvalues = femaleDF$pregnancy_duration, 
        yval = femaleDF[[bmiColumn]], 
        type = "z-scores"
    )
    
    maleDF %>% 
        select(
            child_SentrixID, !!sym(zBmiColumn)
        ) -> maleDF
    femaleDF %>% 
        select(
            child_SentrixID, !!sym(zBmiColumn)
        ) -> femaleDF
    zDF <- rbind(maleDF, femaleDF)
    
    phenoDF %>%
        left_join(
            zDF,
            by = "child_SentrixID"
        ) -> phenoDF
    
}


# Standardize mother height

print(paste0(Sys.time(), "    Standardizing mother height"))

motherDF <- phenoDF %>%
    select(
        child_SentrixID, mother_height, mother_age
    ) %>%
    filter(
        !is.na(mother_age) & !is.na(mother_height)
    )

motherModel <- gamlss(
    formula = mother_height ~ fp(mother_age),
    sigma.formula = ~fp(mother_age),
    family = NO,
    data = motherDF
)

motherDF$z_mother_height = centiles.pred(
    obj = motherModel, 
    xname = "mother_age", 
    xvalues = motherDF$mother_age, 
    yval = motherDF$mother_height, 
    type = "z-scores"
)

motherDF %>% 
    select(
        child_SentrixID, z_mother_height
    ) -> motherDF

phenoDF %>%
    left_join(
        motherDF,
        by = "child_SentrixID"
    ) -> phenoDF


# Standardize father BMI

print(paste0(Sys.time(), "    Standardizing father BMI"))

fatherDF <- phenoDF %>%
    select(
        child_SentrixID, father_bmi, father_age
    ) %>%
    filter(
        !is.na(father_age) & !is.na(father_bmi)
    )

fatherModel <- gamlss(
    formula = father_bmi ~ fp(father_age),
    sigma.formula = ~fp(father_age),
    family = LOGNO,
    data = fatherDF
)

fatherDF$z_father_bmi = centiles.pred(
    obj = fatherModel, 
    xname = "father_age", 
    xvalues = fatherDF$father_age, 
    yval = fatherDF$father_bmi, 
    type = "z-scores"
)

fatherDF %>% 
    select(
        child_SentrixID, z_father_bmi
    ) -> fatherDF

phenoDF %>%
    left_join(
        fatherDF,
        by = "child_SentrixID"
    ) -> phenoDF


# Export DF

print(paste0(Sys.time(), "    Exporting"))

write.table(
    x = phenoDF,
    file = "tmp/phenos",
    quote = F,
    sep = "\t",
    row.names = F
)


# Export docs

print(paste0(Sys.time(), "    Writing docs"))

write(x = "# TrioGen test phenotypes\n", file = docsFile, append = F)
write(x = paste0("Genotyped samples only, ADHD cases and ethnic outliers removed (N = ", nrow(phenoDF), ")\n\n"), file = docsFile, append = T)

write(x = paste0("Phenotypes version V10_1.0.0-190506, standardization using [GAMLSS](https://www.gamlss.com/).\n\n"), file = docsFile, append = T)
write(x = paste0("| Name | variable | Formula | Distribution | Normalization | n |\n"), file = docsFile, append = T)
write(x = paste0("| --------- | ------- | ------------ | ------------- | - |\n"), file = docsFile, append = T)
write(x = paste0("| Standardized Mother height | mother_height | `mother_height ~ fp(mother_age)` | `NO` | `centiles.pred` Z-scores | ", sum(!is.na(phenoDF$z_mother_height)), " |\n"), file = docsFile, append = T)
write(x = paste0("| Standardized Father BMI | z_father_bmi | `father_bmi ~ fp(father_age)` | `LOGNO` | `centiles.pred` Z-scores | ", sum(!is.na(phenoDF$z_father_bmi)), " |\n"), file = docsFile, append = T)
write(x = paste0("| Pregnancy Duration | pregnancy_duration |  |  |  | ", sum(!is.na(phenoDF$pregnancy_duration)), " |\n"), file = docsFile, append = T)
write(x = paste0("| Standardized Placenta Weight | z_placenta_weight | `placenta_weight ~ fp(pregnancy_duration)` per child sex | `BCT` | `centiles.pred` Z-scores | ", sum(!is.na(phenoDF$z_placenta_weight)), " |\n"), file = docsFile, append = T)
write(x = paste0("| Standardized Umbilical Cord Length | z_umbilical_chord_length | `umbilical_chord_length ~ fp(pregnancy_duration)` per child sex | `BCT` | `centiles.pred` Z-scores | ", sum(!is.na(phenoDF$z_umbilical_chord_length)), " |\n"), file = docsFile, append = T)

for (ageI in 0:11) {
    
    variable <- paste0("z_bmi", ageI)
    phenoName <- paste0("BMI at ", timePoints[ageI + 1])
    
    if (ageI < 3) {
        sigmaFormula <- " ~ fp(pregnancy_duration)"
    } else {
        sigmaFormula <- " ~ pregnancy_duration"
    }
    formula <- paste0(variable, sigmaFormula)
    
    write(x = paste0("| ", phenoName, " | ", variable, " | `", formula, "` per child sex | `LOGNO` | `centiles.pred` Z-scores | ", sum(!is.na(phenoDF[[variable]])), " |\n"), file = docsFile, append = T)
    
}
write(x = paste0("| Breastmilk Duration | breastmilk_duration |  |  |  | ", sum(!is.na(phenoDF$breastmilk_duration)), " |\n"), file = docsFile, append = T)
write(x = paste0("| Formula Frequency at 6m | formula_freq_6m |  |  |  | ", sum(!is.na(phenoDF$formula_freq_6m)), " |\n\n"), file = docsFile, append = T)


pheno1 <- "mother_age"
pheno2 <- "mother_height"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Mother Age [Years]",
    labelY = "Mother Height [cm]"
)

write(x = paste0("### Mother Height vs. Mother Age\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "mother_age"
pheno2 <- "z_mother_height"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Mother Age [Years]",
    labelY = "Mother Height [Z-score]"
)

write(x = paste0("### Standardized Mother height vs. Mother Age\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "mother_height"
pheno2 <- "z_mother_height"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Mother Height [cm]",
    labelY = "Mother Height [Z-score]"
)

write(x = paste0("### Standardized Mother Height vs. Mother Height\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "father_age"
pheno2 <- "father_bmi"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Father Age [Years]",
    labelY = "Father BMI [kg/m2]"
)

write(x = paste0("### Father BMI vs. Father Age\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "father_age"
pheno2 <- "z_father_bmi"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Father Age [Z-score]",
    labelY = "Father BMI [kg/m2]"
)

write(x = paste0("### Standardized Father BMI vs. Father Age\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "father_bmi"
pheno2 <- "z_father_bmi"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Father BMI [kg/m2]",
    labelY = "Father BMI [Z-score]"
)

write(x = paste0("### Standardized Father BMI vs. Father BMI\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "pregnancy_duration"
pheno2 <- "umbilical_chord_length"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Pregnancy Duration [Days]",
    labelY = "Umbilical Cord Length [cm]"
)

write(x = paste0("### Umbilical Cord Length vs. Pregnancy Duration\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "pregnancy_duration"
pheno2 <- "placenta_weight"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Pregnancy Duration [Days]",
    labelY = "Placenta Weight [kg]"
)

write(x = paste0("### Placenta Weight vs. Pregnancy Duration\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "umbilical_chord_length"
pheno2 <- "placenta_weight"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Umbilical Cord Length [cm]",
    labelY = "Placenta Weight [kg]"
)

write(x = paste0("### Placenta Weight vs. Umbilical Cord Length\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "umbilical_chord_length"
pheno2 <- "z_umbilical_chord_length"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Umbilical Cord Length [cm]",
    labelY = "Umbilical Cord Length [Z-score]"
)

write(x = paste0("### Standardized Umbilical Cord Length vs. Umbilical Cord Length\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "pregnancy_duration"
pheno2 <- "z_umbilical_chord_length"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Pregnancy Duration [Days]",
    labelY = "Umbilical Cord Length [Z-score]"
)

write(x = paste0("### Standardized Umbilical Cord Length vs. Pregnancy Duration\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "pregnancy_duration"
pheno2 <- "z_placenta_weight"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Pregnancy Duration [Days]",
    labelY = "Placenta Weight [Z-score]"
)

write(x = paste0("### Standardized Placenta Weight vs. Pregnancy Duration\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "placenta_weight"
pheno2 <- "z_placenta_weight"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Placenta Weight [kg]",
    labelY = "Placenta Weight [Z-score]"
)

write(x = paste0("### Standardized Placenta Weight vs. Placenta Weight\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

pheno1 <- "z_umbilical_chord_length"
pheno2 <- "z_placenta_weight"

plotPhenos(
    df = phenoDF,
    pheno1 = pheno1,
    pheno2 = pheno2,
    labelX = "Umbilical Cord Length [Z-score]",
    labelY = "Placenta Weight [Z-score]"
)

write(x = paste0("### Standardized Placenta Weight vs. Standardized Umbilical Cord Length\n"), file = docsFile, append = T)
write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)

for (ageI in 0:11) {
    
    bmiColumn <- paste0("bmi", ageI)
    zBmiColumn <- paste0("z_bmi", ageI)
    
    pheno1 <- "pregnancy_duration"
    pheno2 <- bmiColumn
    
    plotPhenos(
        df = phenoDF,
        pheno1 = pheno1,
        pheno2 = pheno2,
        labelX = "Pregnancy Duration [Days]",
        labelY = paste0("BMI ", timePoints[ageI + 1], " [kg/m2]")
    )
    
    write(x = paste0("### BMI at ", timePoints[ageI + 1], " vs. Pregnancy Duration\n"), file = docsFile, append = T)
    write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)
    
    pheno1 <- "pregnancy_duration"
    pheno2 <- zBmiColumn
    
    plotPhenos(
        df = phenoDF,
        pheno1 = pheno1,
        pheno2 = pheno2,
        labelX = "Pregnancy Duration [Days]",
        labelY = paste0("BMI ", timePoints[ageI + 1], " [Z-score]")
    )
    
    write(x = paste0("### Standardized BMI at ", timePoints[ageI + 1], " vs. Pregnancy Duration\n"), file = docsFile, append = T)
    write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)
    
    pheno1 <- bmiColumn
    pheno2 <- zBmiColumn
    
    plotPhenos(
        df = phenoDF,
        pheno1 = pheno1,
        pheno2 = pheno2,
        labelX = paste0("BMI ", timePoints[ageI + 1], " [kg/m2]"),
        labelY = paste0("BMI ", timePoints[ageI + 1], " [Z-score]")
    )
    
    write(x = paste0("### Standardized BMI at ", timePoints[ageI + 1], " vs. BMI at ", timePoints[ageI + 1], "\n"), file = docsFile, append = T)
    write(x = paste0("![](", pheno1, "-", pheno2, ".png", ")\n\n"), file = docsFile, append = T)
    
}



