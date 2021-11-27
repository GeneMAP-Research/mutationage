#!/usr/bin/env nextflow

nextflow.enable.dsl = 2




def getInputVcf() {
	return channel.fromPath( params.inputDir + params.vcf )
}

def getAffectedSamples() {
	return channel.fromPath( params.inputDir + params.sampleIdsOfAffectedPopulation )
}

def getUnaffectedSamples() {
	return channel.fromPath( params.inputDir + params.sampleIdsOfUnaffectedPopulation )
}

def getAffectedReformattedHapFreqs() {
	return channel.fromPath( params.outputDir + params.variantName + "/hapfreqs/*.affected.haps.reform.freq" )
}

def getUnaffectedReformattedHapFreqs() {
	return channel.fromPath( params.outputDir + params.variantName + "/hapfreqs/*.unaffected.haps.reform.freq" )
				  .ifEmpty { error "Frequency files are empty! Please check and run the workflow again..." }
}

def getAgeEstimateInputFile() {
	return channel.fromFilePairs( params.outputDir + params.variantName + "/params/*.params" , size: 1)
}




process getVariantIdAndPositions() {
	tag "VCF file supplied: ${vcfFile}"
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
		publishDir path: "${params.outputDir}/${params.variantName}/", mode: 'copy'
		path("${params.variantName}.{ld,tags}")
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

		tagkb=\$(( (\$upstream - \$downstream)/1000 ))

		plink \
			--vcf ${vcfFile} \
			--make-bed \
			--keep-allele-order \
			--double-id \
			--chr ${params.chromosomeNumber} \
			--biallelic-only \
			--out temp

		awk '{print \$1,\$1":"\$4,\$3,\$4,\$5,\$6}' temp.bim > temp2.bim

		mv temp2.bim temp.bim

		plink \
			--bfile temp \
			--chr ${params.chromosomeNumber} \
			--from-bp \${downstream} \
			--to-bp \${upstream} \
			--r2 \
			--ld-snp ${params.variantId} \
			--ld-window-r2 ${params.leastLDbetweenTags} \
			--threads ${task.cpus} \
			--double-id \
			--keep-allele-order \
			--out ${params.variantName}

		sed '1d' ${params.variantName}.ld | awk '{print \$6}' > ${params.variantName}.tags
		"""
}

process getListOfPositionsFromTagVariants() {
	input:
		tuple path(rsidChrPosFile), path(ld_file), path(tagVariantsFile)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		path "${params.variantName}-snps.list"
	script:
		"""
		numberOfTags="\$(wc -l ${tagVariantsFile} | awk '{print \$1}')"

		if [ "\$numberOfTags" -gt 25 ]; then
			echo ${params.variantId} > tags.txt
			grep -wv "${params.variantId}" ${tagVariantsFile} | \
				shuf -n 24 >> tags.txt
			tags="tags.txt"
		else
			tags="${tagVariantsFile}"
		fi

		grep \
			-f \${tags} ${rsidChrPosFile} | \
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
			--force-samples \
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
			-f '%POS[ %GT]\\n' | \
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
		template 'getHapsForMutationAge.r'

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
		publishDir path: "${params.outputDir}/${params.variantName}/hapfreqs/", mode: 'copy'
		path "*.reform.freq"
	script:
		template 'fixFreqFiles.r'
}

process getDmleInputFiles() {
	echo true
	input:
		path affected
		path unaffected
	output:
		publishDir path: "${params.outputDir}/${params.variantName}/params/"
		path "${params.variantName}-ageEstimate*.params"
	script:
		template 'makeInputParamsMultipleChainsSingleJob.sh'
}

process getVariantAgeEstimate() {
	tag "${inputFile}"
        label 'ageEstimate'
        label 'dmle'
        //label 'rcran'
	input:
		tuple val(unused_grp_key), path(inputFile), val(unused_chain_number)
	output:
		publishDir path: "${params.outputDir}/${params.variantName}/", mode: 'copy'
		tuple 	path("${inputFile}.output.mutage"), \
				path("${inputFile}.output.mutloc"), \
				path("${inputFile}.mat"), \
				path("${inputFile}.hap"), \
				path("${inputFile}.tre"), \
				path("${inputFile}.sig"), \
				path("${inputFile}.hpf"), \
				path("${inputFile}.dat"), \
				path("${inputFile}.log"), \
				path("${inputFile}.output")
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

process collectAgeEstimateChains() {
	echo true
	input:
		path mutationAge
	output:
		publishDir path: "${params.outputDir}/${params.variantName}/", mode: 'copy'
		path "*.txt"
	script:
		mutAge = mutationAge[0]
		template "collectAgeEstimateChains.r"
}

process collectLocationEstimateChains() {
	echo true
	input:
		path mutationLoc
	output:
		publishDir path: "${params.outputDir}/${params.variantName}/", mode: 'copy'
		path "*.txt"
	script:
		mutLoc = mutationLoc[1]
		template "collectLocationEstimateChains.r"
}
