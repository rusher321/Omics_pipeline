args <- commandArgs(T)
if (length(args) != 4){
	stop("Rscript permanova.R [input_table] [Group_list] <case> <control> [output]
	input_table	: Distance/abundance Matrix (The head column are sample IDs and the first raw are variable names)
	Group_list	: The first raw should define the sample ID and second the belonged group. Make sure the group contain at least 2 levels.
	output    	: output file
	case		: group factor for case
	control		: group factor for control
	* make sure sample of table keep equal with that of the list's
	* make suer your input file contain HEADER. Or your first line will be read as one")
}

profile.f <- args[1]
group.f <- args[2]
out.f <- args[3]
#case <- args[4]
#control <- args[5]
cc <- strsplit(args[4],"-vs-")
case <- cc[[1]][1]
control <- cc[[1]][2] 

library(coin)

#reading
profile.tab <- as.matrix(read.table(profile.f, head=T, row.names=1))
group.tab   <- read.table(group.f, head=T, row.names=1)

# organize the group factor.
group.tab[,1] <- as.factor(group.tab[,1])
fac1 <- case #levels(group.tab[,1])[1]
fac2 <- control #levels(group.tab[,1])[2]

# calculating
out <- matrix(nrow=nrow(profile.tab),ncol=5)
colnames(out) <- list("Wilcox_pvalue","Wilcox_FDR","statistic",paste(fac1,"median",sep="_"), paste(fac2,"median",sep="_"))
rownames(out) <- rownames(profile.tab)

for(i in 1:nrow(profile.tab)){
	value <- as.numeric(profile.tab[i,])
	frame <- data.frame(sample=rownames(group.tab),group=group.tab[,1],value=value)
	frame.cc <- frame[frame$group == fac1 | frame$group == fac2,]
	w.res <- wilcox_test(value~group,frame.cc)

	out[i,4] <- median(frame.cc$value[frame.cc$group == fac1])
	out[i,5] <- median(frame.cc$value[frame.cc$group == fac2])
	out[i,3] <- statistic(w.res)
	out[i,1] <- pvalue(w.res)
}

# FDR
out[,2] <- p.adjust(out[,4],method="BH")

# print
write.table(out, out.f, quote = F, sep = "\t", col.names = NA, append = F)
