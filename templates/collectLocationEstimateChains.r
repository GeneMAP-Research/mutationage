#!/usr/bin/env Rscript

#   COLLECT LOCATION ESTIMATES
#   ==========================
#
#   Author:
#           kevin esoh (eshkev001@myuct.ac.za)
#
#   Date:
#           July 7, 2021
#
#   This template takes as input the mutation location output files
#   and joins them into a single file containing all the chains
#   run during the estimate.
#
###################################################################


getPacklocations <- function()
{
    library(plyr)
}

getJointlocationEstimate <- function() 
{
    iterations <- data.frame(
        seq(
            from=0, 
            to=(${params.mainIterations} - 1), 
            by=1)
        )

    colnames(iterations) <- "ITER"
    jointlocationEstimate <- iterations

    for (run_estimate in seq(from=1, to=${params.numberOfSimultaneousRuns}, by=1)) 
    {
        print(
            paste0(
                "Reading file: ", 
                "${params.outputDir}${params.variantName}-ageEstimate.", 
                run_estimate, 
                ".params.output.mutloc", 
                " ..."
            ), 
            quote=F
        )
        
        locationEstimate <- read.table(
            paste0(
                "${params.outputDir}${params.variantName}-ageEstimate.", 
                run_estimate, 
                ".params.output.mutloc")
            )
    
        locationEstimate_col <- c("ITER")
        locationEstimate_ncol <- ncol(locationEstimate)

        for (clmn in seq(from=1, to=(locationEstimate_ncol - 1))) 
        {
           locationEstimate_col[clmn+1] <- paste0("CHAIN", clmn)
        }
    
        colnames(locationEstimate) <- locationEstimate_col
       
        jointlocationEstimate <- join(jointlocationEstimate, 
                      locationEstimate,
                      type = "inner", 
                      by="ITER")
    }

    return(jointlocationEstimate)
}

formatJointlocationEstimateColumns <- function()
{
    jointlocationEstimate_ncol <- ncol(jointlocationEstimate)
    jointlocationEstimate_col <- c("ITER")

    for (clmn in seq(from=1, to=(jointlocationEstimate_ncol - 1))) {
        jointlocationEstimate_col[clmn+1] <- paste0("CHAIN", clmn)
    }

    colnames(jointlocationEstimate) <- jointlocationEstimate_col

    return(jointlocationEstimate)
}

writeJointlocationFile <- function(j=jointlocationEstimate)
{
    write.table(
        j, 
        "${params.variantName}-jointlocationEstimate.txt", 
        col.names=T, 
        row.names=F, 
        quote=F, 
        sep=" "
    )
}

# suppress warning
defaultW <- getOption("warn");
options(warn = -1);

getPacklocations();
jointlocationEstimate <- getJointlocationEstimate();
jointlocationEstimate <- formatJointlocationEstimateColumns();
writeJointlocationFile(jointlocationEstimate);

print(paste0("Done!"), quote=F);

# unsuppress warning
options(warn = defaultW);