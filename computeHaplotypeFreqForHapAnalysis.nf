#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
	getSamples;
	getSnpList;
	getVcf;
	getHaplotype;
	transposeHapFiles;
	pasteHapFiles;
	collectHapIds;
	getHaplotypeCounts;
	getHaplotypeFreqsFromHapCounts;
} from "${projectDir}/modules/haplotypeAnalysis.nf"

workflow {
	println "\nWorkflow starts here\n"

	pop = getSamples()
	vcfFile = getVcf()
	pop.combine(vcfFile)
	   .set { get_hap_input }

	hapCsvs = getHaplotype(get_hap_input)
	transposed = transposeHapFiles(hapCsvs)
	haps = pasteHapFiles(transposed)

    haps
        .collect()
        .set { collect_hapids_input }    
	hapIds = collectHapIds(collect_hapids_input)

	haps
	    .map { hap, txt -> tuple(hap.baseName, hap, txt) }
	    .set { hap_count_input }

	hapCount = getHaplotypeCounts(hap_count_input)

	hapFreqs = getHaplotypeFreqsFromHapCounts(hapCount)
}

workflow.onComplete { println "\nWorkflow completed successfully!\n" }