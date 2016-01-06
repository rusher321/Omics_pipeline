args <- commandArgs(T)
if (length(args) != 6){
	stop("Rscript permanova.R [input_table] [pheno_table] <case> <control> <phenotype> [output]
	input_table	: Distance/abundance Matrix (The head column are sample IDs and the first raw are variable names)
	pheno_table	: The table of phenotypes
	phenotype	: The specific phenotype for wilcox grouping
	output    	: output file
	case		: phe factor for case
	control		: phe factor for control
	* make sure sample of table keep equal with that of the list's
	* make suer your input file contain HEADER. Or your first line will be read as one")
}

profile.f <- args[1]
    phe.f <- args[2]
     case <- args[3]
  control <- args[4]
      phe <- sub("-",".",args[5])
    out.f <- args[6]

#reading
profile.tab <- as.matrix(read.table(profile.f, head=T, row.names=1))
    phe.tab <- read.table(phe.f, head=T, row.names=1)
    phe.lst <- phe.tab[,colnames(phe.tab)==phe,drop=F]

# organize the phe factor.
phe.lst[,1] <- as.factor(phe.lst[,1])
fac1 <- case #levels(phe.lst[,1])[1]
fac2 <- control #levels(phe.lst[,1])[2]

# calculating
out <- matrix(nrow=nrow(profile.tab),ncol=5)
colnames(out) <- list(paste(fac1,"median",sep="_"), paste(fac2,"median",sep="_"), "Wilcox_statistic", "Wilcox_pvalue", "Wilcox_FDR")
rownames(out) <- rownames(profile.tab)

library(coin)
for(i in 1:nrow(profile.tab)){
	value <- as.numeric(profile.tab[i,])
	frame <- data.frame(sample=rownames(phe.lst),phe=phe.lst[,1],value=value)
	frame.cc <- frame[frame$phe == fac1 | frame$phe == fac2,]
	w.res <- wilcox_test(value~phe,frame.cc)

	out[i,1] <- median(frame.cc$value[frame.cc$phe == fac1])
	out[i,2] <- median(frame.cc$value[frame.cc$phe == fac2])
	out[i,3] <- statistic(w.res)
	out[i,4] <- pvalue(w.res)
}

# FDR
out[,5] <- p.adjust(out[,4],method="BH")

# print
write.table(out, out.f, quote = F, sep = "\t", col.names = NA, row.names = T)
