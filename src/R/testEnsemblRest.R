
# Script to test Ensembl API queries

# Libraries

library(httr)
library(jsonlite)
library(xml2)


# Parameters

server <- "https://rest.ensembl.org"


# Functions

queryEnsembl <- function(
    getQuery,
    nAttempts = 100
    ) {
    
    
    for (try in 1:nAttempts) {
        
        tryCatch(
            {
                r <- GET(paste(server, getQuery, sep = ""), content_type("application/json"))
                
                responseCode <- status_code(r)
                
                if (responseCode == 200) {
                    
                    return(fromJSON(toJSON(content(r))))
                    
                } else {
                    
                    headers <- headers(r)
                    
                    remaining <- headers[["x-ratelimit-remaining"]]
                    reset <- headers[["x-ratelimit-reset"]]
                    
                    if (!is.null(remaining) && !is.null(reset) && remaining <= 1) {
                        
                        Sys.sleep(reset + 1)
                        
                    } else {
                        
                        Sys.sleep(try)
                        
                    }
                }
            },
            error = function(e) {
                
                print(e)
                return(NULL)
                
            }
        )
    }
    
}


# Get variant coordinates

variantQuery <- "/variant_recoder/homo_sapiens/rs5926278?fields=vcf_string"

variantVcf <- queryEnsembl(getQuery = variantQuery)


# List populations available for LD

populationsQuery <- "/info/variation/populations/homo_sapiens?filter=LD"

populationsJson <- queryEnsembl(getQuery = populationsQuery)


# Check LD for a variant

proxyQuery <- "/ld/homo_sapiens/rs62447179/1000GENOMES:phase_3:GBR?r2=0.2;attribs=T"

proxyJson <- queryEnsembl(getQuery = proxyQuery)







