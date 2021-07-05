#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
	getInputVcf;
	getAffectedSamples;
	getUnaffectedSamples;
	getVariantIdAndPositions;
	getVariantIdFile;
	getTagVariants;
	getListOfPositionsFromTagVariants;
	getHaplotypes;
	getTransposedHaplotypes;
	getHaplotypeFrequencies;
	reformatHaplotypeFreqFiles;
} from "${projectDir}/modules/mutationAgeEstimate.nf" 

workflow {

	println "\nworkflow started. Getting haplotype frequencies...\n"

	vcfFile = getInputVcf()
	variantIdFile = getVariantIdFile( params.variantId )
	rsidChrPosFile = getVariantIdAndPositions( vcfFile )
	tagVariantsFile = getTagVariants( vcfFile, variantIdFile, rsidChrPosFile )

	rsidChrPosFile
		.combine( tagVariantsFile )
		.set { positions_input }

	tagVariantPositionsFile = getListOfPositionsFromTagVariants( positions_input )

	// combine these channels because order matters
	tagVariantPositionsFile
		.combine( vcfFile )
		.set { tagVariantsAndVcfFile }

	affectedSampleIds = getAffectedSamples()
	unaffectedSampleIds = getUnaffectedSamples()

	/*  mix sample channels because order does not. But the mixed channel combine with
		tag_vcf channel since order matters when tag_vcf channels are involved */

	affectedSampleIds
		.mix(unaffectedSampleIds)
		.combine ( tagVariantsAndVcfFile )
		.set { sampleTagVariantsAndVcfFile }
	
	haplotypes = getHaplotypes( sampleTagVariantsAndVcfFile )
	transposedHaplotypes = getTransposedHaplotypes( haplotypes ).flatten()
	haplotypeFrequencies = getHaplotypeFrequencies( transposedHaplotypes )

	reformattedHapFreqFiles = reformatHaplotypeFreqFiles( haplotypeFrequencies )

}

workflow.onComplete { println "\nHaplotype frequencies written to '${params.outputDir}' successfully!\n" }