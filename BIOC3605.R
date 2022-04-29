#The assignment4 of BIOC3605--R_Programming#
load("mat.3035563738.RData");
result<-list();
head(mat)
result[["Q1"]]<-dim(mat);
result[["Q2"]]<-apply(mat, 2, sum);
cpm<- t(t(mat)/apply(mat, 2, sum))*10^6;
hist(cpm);
lcpm <- log2(cpm+1);
hist(lcpm,breaks = 50);
boxplot(lcpm);
boxplot(lcpm,outline=F);
lcpm[700:703,]
ll<-apply(lcpm,1,function(x) length(which(x<1)));
ll[1:10]
table(ll)
result[["Q3"]]<-table(ll)
good<-which(ll<=37)
lcpm.good<-lcpm[good,]
dim(lcpm.good)
hist(lcpm.good)
boxplot(lcpm.good,outline=F);
med<-apply(lcpm.good, 2, median);
result[["Q4"]]<-med;
lcpm.norm<-t(t(lcpm.good)-med)+median(med);
b5<-boxplot(lcpm.norm, outline=F, las=2);
result[["Q5"]]<-b5;
hist(lcpm.norm,xlab = "median adjusted log2(cpm+1)", main = "Histogram")

#Do differential expression analysis#
tt<-lapply(rownames(lcpm.norm),function(g){
  t.test(lcpm.norm[g,21:40], lcpm.norm[g,1:20])
})
result[["Q6"]]<-tt[[1]];
p<-sapply(tt, function(x)x$p.value);
h1<-hist(p);
result[["Q7"]]<-h1;
fdr<-p.adjust(p,method="fdr");
h2<-hist(fdr);
result[['Q8']]<-h2;
fc<- sapply(tt,function(x){
  x$estimate[1]-x$estimate[2]
});

de<-which(fdr<0.05)
de.gene<-rownames(lcpm.norm)[de]
writeLines(de.gene)

result[["Q9"]]<-de.gene;

myCol<-array("grey",length(p))
myCol[de]<-"red";
plot(fc,-log10(fdr), pch=16, cex=1.2, col=myCol)
abline(v=c(-1,1), h=-log10(0.05), lty=2)
heatmap(cor(lcpm.norm, method="pearson"),scale="none")
hm<-heatmap(cor(lcpm.norm[de,], method="pearson"),scale="none")
result[["Q10"]]<-hm
heatmap(lcpm.norm[de,],scale='none')
summary(result)
save(result,file="result.3035563738.RData");