args <- commandArgs(T)

train.xF <- args[1]
train.yF <- args[2]
test.xF <- args[3]
test.yF <- args[4]
markerF <- args[5]
cv.fold <- as.numeric(args[6])
cv.step <- as.numeric(args[7])
cv.time <- as.numeric(args[8])
marker.num <- as.numeric(args[9])
prefix <- args[10]

# package
library(randomForest)
library(pROC)

# function
source("/home/fangchao/bin/RANDOMFOREST/rfcv1.R")
source("/home/xiehailiang/bin/ROC/ROC.R")

# crossvalidation
pdf.dir <- paste0(prefix, "_randomForest.pdf")
pdf(pdf.dir, width = 35, height = 7)
par(mfrow = c(1, 5))

save.file <- paste0(prefix,".training.Rdata")
#save.(save.file)
load(save.file)
path <- strsplit(prefix,"/")
watermarker <- paste0(path[[1]][11]," | ",path[[1]][12])
#=======================================================
error.cv <- sapply(train.cv, "[[", "error.cv")
error.cv.rm <- rowMeans(error.cv)
#id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
id <- error.cv.rm == min(error.cv.rm)
error.cv[id, ]
if (marker.num == 0) {
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
}
matplot(train.cv[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cv.time), 
        main = paste("select", marker.num, "Vars"), xlab = "Number of vars", ylab = "CV Error")
lines(train.cv[[1]]$n.var, error.cv.rm, lwd = 2)
abline(v = marker.num, col = "pink", lwd = 2)

mtext(watermarker,side=3,line=3)

# pick marker by corossvalidation
marker.t <- table(unlist(lapply(train.cv, function(x) {
  lapply(x$res, "[", 1:marker.num)
})))
marker.t <- sort(marker.t, d = T)
names(marker.t) <- colnames(train.x)[as.numeric(names(marker.t))]
marker.dir <- paste0(prefix, "_marker.txt")
write.table(marker.t, marker.dir, col.names = F, sep = "\t", quote = F)
marker.p <- names(marker.t)[1:marker.num]

marker.s.dir <- paste0(prefix, "_select.marker.txt")
write.table(marker.t[1:marker.num],marker.s.dir,col.names = F, sep = "\t", quote = F)

# train model
set.seed(0)
train.rf <- randomForest(train.x[, marker.p], train.y, importance = T)
train.p <- predict(train.rf, type = "prob")
boxplot(train.p[, 2] ~ train.y, col = 2:3, main = "Probability", names = train.l)
pr.dir <- paste0(prefix, "_train_probability.txt")
write.table(train.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)

# train ROC
plot_roc(train.y, train.p[, 2])

# test predict
test.p <- predict(train.rf, test.x[, marker.p], type = "prob")
pr.dir <- paste0(prefix, "_test_probability.txt")
write.table(test.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)

# predict plot
p.col <- ifelse(is.na(test.y), 4, as.numeric(test.y) + 1)
plot(rank(test.p[, 2]), test.p[, 2], col = p.col, pch = 16, xlab = "", ylab = "Probability", main = "Testset")
txt <- train.l
if (length(test.l) > 2) {
  txt <- c(txt, "the rest")
}
legend("bottomright", txt, col = 2:4, pch = 16)
abline(h = 0.5)

# test ROC
plot_roc(test.y, test.p[, 2])
dev.off() 
