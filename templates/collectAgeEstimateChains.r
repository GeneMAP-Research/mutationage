#!/usr/bin/env Rscript

#   COLLECT AGE ESTIMATES
#   =====================
#
#   Author:
#           kevin esoh (eshkev001@myuct.ac.za)
#
#   Date:
#           July 7, 2021
#
#   This template takes as input the mutation age output files
#   and joins them into a single file containing all the chains
#   run during the estimate.
#
###############################################################


getPackages <- function()
{
    library(plyr)
}

getJointAgeEstimate <- function() 
{
    iterations <- data.frame(
        seq(
            from=0, 
            to=(${params.mainIterations} - 1), 
            by=1)
        )

    colnames(iterations) <- "ITER"
    jointAgeEstimate <- iterations

    for (run_estimate in seq(from=1, to=${params.numberOfSimultaneousRuns}, by=1)) 
    {
        print(
            paste0(
                "Reading file: ", 
                "${params.outputDir}${params.variantName}-ageEstimate.", 
                run_estimate, 
                ".params.output.mutage", 
                " ..."
            ), 
            quote=F
        )
        
        ageEstimate <- read.table(
            paste0(
                "${params.outputDir}${params.variantName}-ageEstimate.", 
                run_estimate, 
                ".params.output.mutage")
            )
    
        ageEstimate_col <- c("ITER")
        ageEstimate_ncol <- ncol(ageEstimate)

        for (clmn in seq(from=1, to=(ageEstimate_ncol - 1))) 
        {
    	   ageEstimate_col[clmn+1] <- paste0("CHAIN", clmn)
        }
    
        colnames(ageEstimate) <- ageEstimate_col
       
        jointAgeEstimate <- join(jointAgeEstimate, 
    			      ageEstimate,
                      type = "inner", 
    			      by="ITER")
    }

    return(jointAgeEstimate)
}

formatJointAgeEstimateColumns <- function()
{
    jointAgeEstimate_ncol <- ncol(jointAgeEstimate)
    jointAgeEstimate_col <- c("ITER")

    for (clmn in seq(from=1, to=(jointAgeEstimate_ncol - 1))) {
        jointAgeEstimate_col[clmn+1] <- paste0("CHAIN", clmn)
    }

    colnames(jointAgeEstimate) <- jointAgeEstimate_col

    return(jointAgeEstimate)
}

writeJointAgeFile <- function(j=jointAgeEstimate)
{
    write.table(
        j, 
        "${params.variantName}-jointAgeEstimate.txt", 
        col.names=T, 
        row.names=F, 
        quote=F, 
        sep=" "
    )
}

# suppress warning
defaultW <- getOption("warn");
options(warn = -1);

getPackages();
jointAgeEstimate <- getJointAgeEstimate();
jointAgeEstimate <- formatJointAgeEstimateColumns();
writeJointAgeFile(jointAgeEstimate);

print(paste0("Done!"), quote=F);

# unsuppress warning
options(warn = defaultW);