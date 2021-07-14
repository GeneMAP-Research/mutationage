#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
	println "\nPlotting starts here\n"

	hapFreqs = getHapFrequencies().collect()
	hapIds = getHapIds()

	hapPlot = plotHaplotypeFreqs(hapFreqs, hapIds)

}

def getHapFrequencies() {
	return channel.fromPath( params.outputDir + "*.freq" )
}

def getHapIds() {
	return channel.fromPath( params.outputDir + "${params.outPrefix}.haps.ids" )
}

process plotHaplotypeFreqs() {
	input:
		path hapFreqs
		path hapIds
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple \
			path("${params.outPrefix}-hap.pdf"), \
			path("${params.outPrefix}-transpoded-hap-matrix.csv"), \
			path("${params.outPrefix}-hap-matrix.csv")
	script:
		template "plotHapFreqsFromHapAnalysis.r"
}