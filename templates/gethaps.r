#!/usr/bin/env Rscript

f <- read.table("${haplotypeFiles}", h=F); 
f\$cm <- c(0)

#--- Compute genetic map, each position compared to the previous (makes sense biologically)
#for (i in 1:(length(f\$V1)-1)) { f\$cm[i+1] <- (f\$V1[i+1]-f\$V1[i])/1000000}


#--- Compute genetic map, each position compared to the first (gives increasing order)
for (i in 1:(length(f\$V1)-1)) { f\$cm[i+1] <- (f\$V1[i+1]-f\$V1[1])/1000000}

f\$V1 <- formatC(signif(f\$cm,digits=5), digits=5,format="fg", flag="#")

f\$cm <- NULL

write.table(t(f), file="${haplotypeFiles.baseName}.haps.transposed", col.names=F, row.names=F, quote=F, sep =" ")

