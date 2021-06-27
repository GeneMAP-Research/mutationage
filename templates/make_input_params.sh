#!/usr/bin/env bash


diseaseHaps=\$(awk '{sum += \$1} END {print sum}' "${affected}")
nhaps=\$(( \$diseaseHaps - 1 ))
echo "\nThere are \$nhaps haplotypes in "${affected}""

number_of_chrom=\${nhaps}
number_of_loci=\$( head -1 ${affected} | cut -f2- -d' ' | tr ' ' '\n' | wc -l)
pop_growth_rate=${params.popGrowthRate}
variant_allele_freq=${params.variantAlleleFreq}
burn_iterations=${params.burnIterations}
main_iterations=${params.mainIterations}
nchains=${params.numberOfSimultaneousRuns}


echo """Data as genotypes? Yes = 1, No = 0
0
Genetic model: Dominant=0,Recessive=1
1
Read old file?: (0=no, 1=yes):
0
Use fixed random seed?:(0 = no (=random), negative integer = yes (=fixed), and is the seed):
0
# chromosomes (N):
\${number_of_chrom}
# loci per chromosome (L):
\${number_of_loci}
Numbers of haplotypes in the normal(base) pop.:""" > "${params.variantName}-ageEstimate.params"



sed '\$d' ${unaffected} >> "${params.variantName}-ageEstimate.params"



echo "Map distances:" >> "${params.variantName}-ageEstimate.params"



tail -1 ${affected} | cut -f2- -d' ' >> "${params.variantName}-ageEstimate.params"



echo """Run simulation?:
0
Mutation location (only used for simulated data):
-0.0008559
Mutation's low and high boundaries
0 1
# simultaneous runs:
\${nchains}
Starting value(s) for recdist. for each simul. run (-99 for random):
-99
Population growth rate:
\${pop_growth_rate}
Proportion of population sampled:
\${variant_allele_freq}
Iterate ancestral states,mutation age,mutation location, allele freq. (0=no, 1=yes):
1 1 1 1
Flip (potentially) all loci? (0=no, 1=yes):
0
Adjustment level for tree, recdist, ancestral, and internal states,alleles:
1.0 3.0 0.005 0.5 0.5 0.1 1
Burn-in iterations:
\${burn_iterations}
Iterations:
\${main_iterations}
Screen update and file update intervals:
100 100 0 0
Number of histogram bars:
200
Alpha level for recdist histogram:
0.05
Mutation age (-99 for random):
-99
Mutation age boundaries:
0 50000
Star genealogy (0=no, 1=yes):
0
Loci for the root  (1xL) (-99):
\$( for i in \$(seq 1 \${number_of_loci}); do echo 1; done | tr '\\n' ' ' | sed 's/ \$/\\n/g' )
Frequency, and loci for the tip chromosomes (?x(L+1)):""" >> "${params.variantName}-ageEstimate.params"





sed '\$d' ${affected} >> "${params.variantName}-ageEstimate.params"




echo """Use sequence weights?
0
Weights for exons,introns,non-genes
1 0.17 0.02
Input file:
ncbi.txt
301500
301700
303234
311111
312456""" >> "${params.variantName}-ageEstimate.params"