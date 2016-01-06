args <- commandArgs(T)
if (length(args) != 8) {stop("Rscript wilcox.test.KW.R [input1] [input2]
input1: profile.tab_profile
input2: Phenotype")
}

profile.f <- args[1]
    phe.f <- args[2]
  control <- args[3]
    case1 <- args[4]
    case2 <- args[5]
      phe <- sub("-",".",args[6])
	  adj <- args[7]
    out.f <- args[8]

profile.tab <- as.matrix(read.table(profile.f,head=T))
profile.tab <- profile.tab[rowSums(profile.tab)!=0,]
profile.nrow <- nrow(profile.tab)

### phenotypes...
phe.tab <- read.table(phe.f, head=T, row.names=1)
phe.lst <- phe.tab[,colnames(phe.tab)==phe,drop=F]

# organize the phe factor
phe.lst[,1] <- as.factor(phe.lst[,1])
#bf <- as.factor(phe.tab$HT)
if(adj!="NA"){
	adj.lst <- phe.tab[,colnames(phe.tab)==adj,drop=F]
	adj.lst[,1] <- as.factor(adj.lst[,1])
}


# result...
out <- matrix(0,nrow=profile.nrow,ncol=11)
rownames(out) <- rownames(profile.tab)
colnames(out) <- c(paste0(control,"_median"),paste0(case1,"_median"),paste0(case2,"_median"),"KW.p",paste0(control,"_dunn.p"),paste0(case1,"_dunn.p"),paste0(case2,"_dunn.p"),paste0(control,"_Z.Stat"),paste0(case1,"_Z.Stat"),paste0(case2,"_Z.Stat"),"BH.FDR")
# test...

library(coin)
library(dunn.test)

for(i in 1:profile.nrow){
	value <- as.numeric(profile.tab[i,])
	if (adj=="NA"){
		frame <- data.frame(sample=rownames(phe.lst),phe=phe.lst[,1],value=value)
		frame.cc <- frame[frame$phe == control | frame$phe == case1 | frame$phe == case2,]
		kw.res <- kruskal_test(value ~ phe, frame.cc)
	}else {
		frame <- data.frame(sample=rownames(phe.lst),phe=phe.lst[,1],value=value,adjust=adj.lst[,1])
		frame.cc <- frame[frame$phe == control | frame$phe == case1 | frame$phe == case2,]
		kw.res <- kruskal_test(value ~ phe|adjust, frame.cc)
	}
	out[i,1] <- median(frame.cc$value[frame.cc$phe == control])
	out[i,2] <- median(frame.cc$value[frame.cc$phe == case1])
	out[i,3] <- median(frame.cc$value[frame.cc$phe == case2])
	out[i,4] <- pvalue(kw.res)

	dunn.res <- dunn.test(frame.cc$value, g=frame.cc$phe, method="bh", kw=FALSE,rmc=FALSE)
	out[i,5:10]<-c(dunn.res$P.adjusted[1],dunn.res$P.adjusted[2],dunn.res$P.adjusted[3],-dunn.res$Z[1],-dunn.res$Z[2],-dunn.res$Z[3])
}
# FDR
out[,11] <- p.adjust(out[,4],method="BH")
# write table...
write.table(out, out.f, sep="\t", quote=F, col.names=NA, row.names = T)
