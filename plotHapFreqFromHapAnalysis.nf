#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
	getHapFrequencies;
	getHapIds;
	plotHaplotypeFreqs;
} from "${projectDir}/modules/haplotypeAnalysis.nf"

workflow {
	println "\nPlotting starts here\n"

	hapFreqs = getHapFrequencies().collect()
	hapIds = getHapIds()
	hapPlot = plotHaplotypeFreqs(hapFreqs, hapIds)

}

workflow.onComplete { println "\nWorkflow completed successfully!\n" }