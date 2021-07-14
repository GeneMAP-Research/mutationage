#!/usr/bin/env Rscript

tsv_out <- gsub(".csv",
				".transposed.tsv", 
				"${hapCsv}")

txt_out <- gsub(".csv",
				".transposed.txt", 
				"${hapCsv}")

f <- read.table("${hapCsv}",
				h=F, 
				sep=",")
ft <- t(f)

write.table(ft, 
			file=tsv_out, 
			col.names=F, 
			row.names=F, 
			quote=F, 
			sep="\t")

write.table(ft, 
			file=txt_out, 
			col.names=F, 
			row.names=F, 
			quote=F, 
			sep="")