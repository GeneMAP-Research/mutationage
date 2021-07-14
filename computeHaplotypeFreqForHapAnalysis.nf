#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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


def getSamples() {
	return channel
	            .fromFilePairs( params.hapSamplesDir + "*.txt", size: 1 )
	            .map { popFile -> tuple( popFile ) }
}

def getSnpList() {
	return channel.fromPath( params.snpList )
}

def getVcf() {
	return channel.fromPath( params.inputDir + params.vcf )
	              .ifEmpty { error "\nNo VCF file found! " + \
	                        "Please make sure you have provided " + \
	                        "a VCF file in the config file and make sure to " + \
	                        "provide the input directory\n" }
}


process getHaplotype() {
	tag "processing ${popName}"
	input:
	    tuple val(popName), path(popFile), path(vcfFile)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple val(popName), path("${popName}.haps.csv")
	script:
		"""
    	bcftools \
    		view \
    		-S ${popFile} \
    		--force-samples \
    		-t \$(cat ${params.snpList}) \
    		${vcfFile} | \
	    		bcftools \
	    		query \
	    		-f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' | \
	    			cut -f3-4,6- | \
	    				sed 's/|/,/g' | \
	    					sed 's/\\//,/g' | \
	    						sed 's/\t/,/g' > \
	    						${popName}.haps.csv		
		"""
}

process transposeHapFiles() {
	tag "processing ${popName}"
	input:
		tuple val(popName), path(hapCsv)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple val(popName), path("*.transposed.tsv"), path("*.transposed.txt")
	script:
		template "getHapsForHapAnalysis.r"
}

process pasteHapFiles() {
	tag "processing ${popName}"
	input:
		tuple val(popName), path(tsvFile), path(txtFile)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple path("${popName}.haps"), path("${popName}.haps.txt")
	script:
		"""
		paste ${tsvFile} ${txtFile} > ${popName}.haps
    	sed '1,2d' ${popName}.haps | \
	    	cut -f5 | \
	    		sort | \
	    			uniq > ${popName}.haps.txt
		"""
}


process collectHapIds() {
	input:
		path input
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		path "${params.outPrefix}.haps.ids"
	script:
		"""
		for hapFile in ${input}; do
			if [ \${hapFile##*.} == "txt" ]; then
				cat \${hapFile}
			fi
		done | \
			sort | \
				uniq | \
					awk '{print "\\""\$1"\\""}' | \
						sed '1 i Hap' > \
						${params.outPrefix}.haps.ids
		"""
}

process getHaplotypeCounts() {
	tag "processing ${popName}"
	input:
		tuple val(popName), path(hap), path(txt)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple val(popName), path("${popName}.haps.frq"), path("${popName}.haps.count")
	script:
		"""
		head -2 ${hap} | \
	    	awk '{print "Count",\$0}' > \
	    	${popName}.haps.frq
    
    	sed '1,2d' ${hap} | \
	    	sort | \
	    		uniq -c | \
	    			sort \
	    			-g \
	    			-k6 >> \
	    			${popName}.haps.frq
    
    	sed '2d' ${popName}.haps.frq | \
	    	awk '{print \$1,\$2,\$3,\$4,\$5,"\\""\$6"\\""}' > ${popName}.haps.count
		"""
}

process getHaplotypeFreqsFromHapCounts() {
	tag "processing ${popName}"
	input:
		tuple val(popName), path(hapFrq), path(hapCount)
	output:
		publishDir path: "${params.outputDir}", mode: 'copy'
		tuple val(popName), path("${popName}.haps.freq")
	script:
		template "getHapsFreqsForHapAnalysis.r"
}