#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
	getAffectedReformattedHapFreqs;
	getUnaffectedReformattedHapFreqs;
	getDmleInputFiles;
} from "${projectDir}/modules/mutationAgeEstimate.nf"

workflow {

	println "\nworkflow started. Making DMLE input file...\n"

	affected = getAffectedReformattedHapFreqs()
	unaffected = getUnaffectedReformattedHapFreqs()

	dmleInput = getDmleInputFiles( affected, unaffected )

}

workflow.onComplete { println "Done making DMLE input file successfully!\n" }