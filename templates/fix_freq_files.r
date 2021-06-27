#!/usr/bin/env Rscript

freq <- read.table("${haplotypeFrequencies}",h=F)
write.table(freq, file="${haplotypeFrequencies.baseName}.reform.freq", col.names=F, row.names=F, quote=F, sep=" ")
