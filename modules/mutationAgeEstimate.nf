#!/usr/bin/env nextflow

nextflow.enable.dsl = 2




def getInputVcf() {
	return channel.fromPath( params.inputDir + params.vcfFile )
}

def getAffectedSamples() {
	return channel.fromPath( params.inputDir + params.sampleIdsOfAffectedPopulation )
}

def getUnaffectedSamples() {
	return channel.fromPath( params.inputDir + params.sampleIdsOfUnaffectedPopulation )
}

def getAffectedReformattedHapFreqs() {
	return channel.fromPath( params.outputDir + '/*.affected.haps.reform.freq' )
}

def getUnaffectedReformattedHapFreqs() {
	return channel.fromPath( params.outputDir + '/*.unaffected.haps.reform.freq' )
}

def getAgeEstimateInputFile() {
	return channel.fromFilePairs( params.outputDir + params.variantName + "*.params" , size: 1)
}




process getVariantIdAndPositions() {
	input:
		path vcfFile
	output:
		path "${params.variantName}-rsid-chr-pos.txt"
	script:
		"""
		bcftools \
			query \
			-f '%ID\t%CHROM:%POS\n' \
			${vcfFile} > \
			${params.variantName}-rsid-chr-pos.txt
		"""
}

process getVariantIdFile() {
	tag "Variant ID supplied: ${variantId}"
	input:
		val variantId
	output:
		path "${params.variantName}.rsid"
	script:
		"""
		echo ${params.variantId} > ${params.variantName}.rsid
		"""
}

process getTagVariants() {
	input:
		path vcfFile
		path variantIdFile
		path rsidChrPosFile
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		path("${params.variantName}.{tags,tags.list}")
	script:
		"""
		halfInterval=\$(( ${params.mutationRegionSize}/2 ))

		mutposition=\$(grep -w -f ${variantIdFile} ${rsidChrPosFile} | awk '{print \$2}' | cut -f2 -d':')

		if [ \$mutposition -le \$halfInterval ]; then
			downstream=\$(( \$mutposition - 500000 ))
			upstream=\$(( \$mutposition + 500000 ))
		else
			downstream=\$(( \$mutposition - \$halfInterval ))
			upstream=\$(( \$mutposition + \$halfInterval ))
		fi

		plink \
			--vcf ${vcfFile} \
			--chr ${params.chromosomeNumber} \
			--from-bp \${downstream} \
			--to-bp \${upstream} \
			--show-tags ${variantIdFile} \
			--tag-r2 ${params.leastLDbetweenTags} \
			--tag-kb 2000 \
			--list-all \
			--threads ${task.cpus} \
			--double-id \
			--keep-allele-order \
			--out ${params.variantName}
		"""
}

process getListOfPositionsFromTagVariants() {
	input:
		tuple path(rsidChrPosFile), path(tagVariantsFile), path(tagVariantslist)
	output:
		path "${params.variantName}-snps.list"
	script:
		"""
		grep \
			-f ${tagVariantsFile} ${rsidChrPosFile} | \
			awk '{print \$2}' | \
			tr '\\n' ',' | \
			sed 's/,\$/\\n/g' > "${params.variantName}-snps.list"
		"""
}

process getHaplotypes() {
	tag "${sampleIds}"
	input:
		tuple path(sampleIds), path(tagVariantPositionsFile), path(vcfFile)
	output:
		path "*.haps"
	script:
		"""
		bcftools \
			index \
			--threads ${task.cpus} \
			-ft \
			${vcfFile}

		bcftools \
			view \
			-v snps \
			--threads ${task.cpus} \
			-k \
			-m2 \
			-M2 \
			-S ${sampleIds} \
			-r \$(cat ${tagVariantPositionsFile}) \
			${vcfFile} | \
		bcftools \
			query \
			-H \
			-f '%POS[ %GT]\n' | \
		sed 's/1|1/2 2/g' | \
		sed 's/1|0/2 1/g' | \
		sed 's/0|1/1 2/g' | \
		sed 's/0|0/1 1/g' | \
		sed '1d' > ${sampleIds.baseName}.haps
		"""
}

process getTransposedHaplotypes() {
	tag "${haplotypeFiles}"
	input:
		path haplotypeFiles
	output:
		path "${haplotypeFiles.baseName}.haps.transposed"
	script:
		template 'gethaps.r'

}

process getHaplotypeFrequencies() {
	tag "${transposedHaplotypes}"
	input:
		path transposedHaplotypes
	output:
		path "*.freq"
	script:

		freq_out = transposedHaplotypes.baseName
	
		"""
		sort ${transposedHaplotypes} | \
			uniq -c | \
			sort -gr -k1 > "${freq_out}.freq"
		"""

}

process reformatHaplotypeFreqFiles() {
	tag "${haplotypeFrequencies}"
	input:
		path haplotypeFrequencies
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		path "*.reform.freq"
	script:
		template 'fix_freq_files.r'
}

process getDmleInputFiles() {
	echo true
	input:
		path affected
		path unaffected
	output:
		publishDir path: "${params.outputDir}"
		path "${params.variantName}-ageEstimate.*.params"
	script:
		template 'make_input_params.sh'
}

process getVariantAgeEstimate() {
	tag "${inputFile}"
	input:
		tuple val(unused_grp_key), path(inputFile), val(unused_chain_number)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		path "${inputFile}.{mat,hap,tre,sig,hpf,dat,log,output*}"
	script:
		"""
		DMLE+2.2 ${inputFile}

		if [ \$? -eq 0  ]; then
			cut -f1,2,4 "${inputFile}.dat" -d' ' | \
				sed '1 i ITER MUTLOC MUTAGE' > "${inputFile}.output"
		fi
		
		awk '{print \$1,\$2,\$11,\$20,\$29,\$38,\$47,\$56,\$65,\$74,\$83}' "${inputFile}.dat" > "${inputFile}.output.mutloc"
		awk '{print \$1,\$4,\$13,\$22,\$31,\$40,\$49,\$58,\$67,\$76,\$85}' "${inputFile}.dat" > "${inputFile}.output.mutage"
		"""
}