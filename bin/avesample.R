args <- commandArgs(T) 
if (length(args) != 8) {stop("Rscript randomForest14.R [input1] [input2] [input3] [input4] [case1] [case2] [control] [profix]
	input1: groups info
	input2: picked list/significant compounds
	input3: profile
	case: the case group to be distinguish
	control: the control group to be distinguish
	prefix: prefix
	make sure sample of phe.pr and dist.pr have the same rownames")
}

group.file <- args[1]
todo <- sub("-",".",args[2])
pick.file <- args[3]
profile.file <- args[4]
case1 <- args[5]
case2 <- args[6]
control <- args[7]
st <- args[8]

library(sampling)


dat1 <- read.table(profile.file, head=T,sep="\t",row.names=1)
conf0 <- read.table(group.file,head=T,row.names=1,sep="\t")
conf1 <- conf0[,colnames(conf0)==todo,drop=F]
sign <- read.table(pick.file, head=T,sep="\t",row.names=1)
#sign <- dat1

cN.prof <- colnames(dat1)
rN.prof <- rownames(dat1)

rN.conf1 <- rownames(conf1)
rN.sign <- rownames(sign)
#if(pick.file != profile.file){
#	rN.sign <- rownames(sign)[which(sign[,1]<0.05)]
#}
rN.conf1 <- rownames(conf1)
sid <- intersect(cN.prof ,rN.conf1)
cid <- intersect(rN.prof ,rN.sign)

dat1 <- dat1[pmatch(cid,rN.prof),pmatch(sid,cN.prof)]
#conf1 <- conf1[pmatch(gid,rN.conf1),]
conf2 <-conf1
dat2 <- dat1

#dat2 <- dat1[,conf1[,1]== control | conf1[,1]== case ]
#conf2 <- conf1[conf1[,1]== control | conf1[,1]== case ,,drop=F]
#conf2[,1]=as.factor(as.character(conf2[,1]))
phe = conf1[,1]
#phe<-sub(control ,"0" ,phe)
#phe<-sub(case ,"1" ,phe)
phe<-as.factor(phe)

dat2t <- as.data.frame(t(dat2))
dat2t$phe <- phe
avesampling <- strata(dat2t,"phe",size=rep(70,length(table(dat2t$phe))),method = "srswor")
train<- getdata(dat2t,avesampling)[,1:ncol(dat2t)]

tr.fun <- function (case,control){
	train.x <- train[train$phe == control | train$phe == case,-ncol(train)]
	train.y <- train[train$phe == control | train$phe == case,ncol(train),drop=F]
	write.table(t(train.x),paste0(st,".",case,"-vs-",control,".train_x.xls"),quote=F,col.names=NA,sep="\t")
	write.table(train.y,paste0(st,".",case,"-vs-",control,".train_y.xls"),quote=F,col.names=NA,sep="\t")
	return(train.x)
}

train.x1 <- tr.fun(case1,control)
if (case2 != "F"){
	train.x2 <- tr.fun(case2,control)
	train.x3 <- tr.fun(case2,case1)
}

#########iest.set###########
dat3 <- as.data.frame(t(dat1))
dat3$phe <- conf1[,1]
#test.order <- setdiff(rownames(dat3),rownames(train.x))
#test <- dat3[pmatch(test.order, rownames(dat3)),]

#test.x <- test[,-ncol(test)]
#test.y <- test[,ncol(test),drop=F]
w.fun <- function (x,case,control){
	test.order <- setdiff(rownames(dat3),rownames(x))
	test <- dat3[pmatch(test.order, rownames(dat3)),]

	test.x <- test[,-ncol(test)]
	test.y <- test[,ncol(test),drop=F]
	write.table(t(test.x),paste0(st,".",case,"-vs-",control,".test_x.xls"),quote=F,col.names=NA,sep="\t")
	write.table(test.y,paste0(st,".",case,"-vs-",control,".test_y.xls"),quote=F,col.names=NA,sep="\t")
	return(test.y)
}

y1 <- w.fun(train.x1,case1,control)
if (case2 != "F"){
	y2 <- w.fun(train.x2,case2,control)
	y3 <- w.fun(train.x3,case2,case1)
}
