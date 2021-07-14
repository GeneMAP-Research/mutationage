#!/usr/bin/env Rscript

freq_out <- gsub(".count", 
				 ".freq", 
				 "${hapCount}")

f <- read.csv2(
				"${hapCount}",
				h=T,
				quote="",
				sep=" ")

n <- ncol(f)
f\$Hap <- f[,n]
f[,n] <- NULL
f\$Freq <- ((f\$Count)/(sum(f\$Count)))
hap_frq <- data.frame(f\$Hap, f\$Freq)
rownames(hap_frq) <- f\$Hap

write.table(f, 
			file=freq_out, 
			col.names=T, 
			row.names=F, 
			sep="\t", 
			quote=F)
