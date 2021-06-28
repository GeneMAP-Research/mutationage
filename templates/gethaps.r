#!/usr/bin/env Rscript


#					COMPUTE GENETIC MAP FOR DMLE INPUT
#				==========================================
#
#	This template is called from the `getTransposedHaplotypes` function.
#	It takes as input a .hap file emitted by the function above and computes
#	a 'DMLE' genetic map:
#
#		*	The map starts with the first variant position as zero (0). 
#		*	The distance between the first variant and each subsequent 
#			position is computed by subtracting the first position from
#			each subsequent position.
#		*	The distance is then converted to centi Morgans (cM) by using
#			the analogy 
#
#					1 cM = 1 Mb
#
#		*	Example:
#					if 		1 cM --> 1,000,000 bp
#					then	X cM --> 5,248,232 bp
#
#							1,000,000 bp * X cM = 5,248,232 bp * 1 cM
#
#							X cM = ( 5,248,232 bp * 1 cM ) / 1,000,000 bp
#
#
################################################################################



f <- read.table("${haplotypeFiles}", h=F)

f\$cm <- c(0)

for (i in 1:(length(f\$V1)-1)) { 

	f\$cm[i+1] <- (f\$V1[i+1]-f\$V1[1])/1000000

}

f\$V1 <- formatC(signif(f\$cm,digits=5), 
				 digits=5,
				 format="fg", 
				 flag="#")

f\$cm <- NULL

write.table(t(f), 
			file="${haplotypeFiles.baseName}.haps.transposed", 
			col.names=F, 
			row.names=F, 
			quote=F, 
			sep =" ")
