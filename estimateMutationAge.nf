#!/usr/bin/env nextfow

nextflow.enable.dsl = 2

include {
	getAgeEstimateInputFile;
	getVariantAgeEstimate;
	collectAgeEstimateChains;
	collectLocationEstimateChains;
} from "${projectDir}/modules/mutationAgeEstimate.nf"

workflow {
	println "\nNow estimating ${params.variantId} age...\n"

	//inputFile = getAgeEstimateInputFile()

	nchains = channel.of(1..params.numberOfSimultaneousRuns)
	nchains
		.map { chain -> tuple("${params.variantName}-ageEstimate.${chain}", chain) } 		// create a group key similar to the one generated in paramsFiles bellow 
		.set { key_chains }

	paramsFiles = getAgeEstimateInputFile()
	paramsFiles
		.map { group_key, file -> tuple( group_key, file.first() ) }
		.set { key_params }

	key_params
		.join( key_chains )
		.set { age_estimate_input }


	mutationAge = getVariantAgeEstimate( age_estimate_input )

	collectAgeEstimateChains(mutationAge)
	collectLocationEstimateChains(mutationAge)

	//log.info "This is a test log info"

}

workflow.onComplete { println "Done estimating ${params.variantId} age...!\n" }
