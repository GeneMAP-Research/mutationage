#!/usr/bin/env nextfow

nextflow.enable.dsl = 2

include {
	getAgeEstimateInputFile;
	getVariantAgeEstimate;
} from "${projectDir}/modules/mutationAgeEstimate.nf"

workflow {
	println "\nNow estimating ${params.variantId} age...\n"

	inputFile = getAgeEstimateInputFile()
	mutationAge = getVariantAgeEstimate( inputFile )

}

workflow.onComplete { println "Done estimating ${params.variantId} age...!\n" }
