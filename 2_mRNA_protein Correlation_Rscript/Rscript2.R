##Example
library(data.table)

## protein ID convert to gene symbol (using UniProt Retrieve/ID mapping webtool)
LDLIST<-as.matrix(fread("ID_LIST.txt",sep="\t",header=F))
##Protein
Protein<-as.matrix(fread("ProteinExpression.txt",sep="\t",header=F))
##mRNA
Transcript<-as.matrix(fread("mRNAExpression.txt",sep="\t",header=F))


k<-1
temp<-NULL
while(k<=length(LDLIST[,1]))
{

      A_pro<-match(LDLIST[k,1],Protein[,1],nomatch=0)
      B_tran<-match(LDLIST[k,2],Transcript[,1],nomatch=0)

    if(A_pro!=0 && B_tran!=0)
    {
        temp<-rbind(temp,cbind(t(Protein[A_pro,]),t(Transcript[B_tran,])))
    }

k<-k+1
}


#write.table(temp,"Protein_RNA_data.txt",sep="\t",col.names=F,row.names=F,append=F,quote=F)

#####################Sample_Wise_correlation

Protein_RNA<-as.matrix(fread("Protein_RNA_data.txt",sep="\t",header=F))

##Group1
samA<-2:24
samB<-91:113

k<-1
Group1_sample<-NULL
while(k<=length(samA))
{

Group1_sample<-c(Group1_sample,cor.test(as.numeric(Protein_RNA[,samA[k]]), as.numeric(Protein_RNA[,samB[k]]), method = "spearm",exact=F)$ estimate)

k<-k+1
}


##Group2 
samA<-25:82
samB<-114:171

k<-1
Group2_sample<-NULL
while(k<=length(samA))
{

Group2_sample<-c(Group2_sample,cor.test(as.numeric(Protein_RNA[,samA[k]]), as.numeric(Protein_RNA[,samB[k]]), method = "spearm",exact=F)$ estimate)

k<-k+1
}

##Group3 
samA<-83:89
samB<-172:178

k<-1
Group3_sample<-NULL
while(k<=length(samA))
{

Group3_sample<-c(Group3_sample,cor.test(as.numeric(Protein_RNA[,samA[k]]), as.numeric(Protein_RNA[,samB[k]]), method = "spearm",exact=F)$ estimate)

k<-k+1
}


##Group4  all samples
samA<-2:89
samB<-91:178

k<-1
Group4_sample<-NULL
while(k<=length(samA))
{

Group4_sample<-c(Group4_sample,cor.test(as.numeric(Protein_RNA[,samA[k]]), as.numeric(Protein_RNA[,samB[k]]), method = "spearm",exact=F)$ estimate)

k<-k+1
}


library(ggpubr)
library(ggplot2)
library(gridExtra)
library(grid)




TTdata<-data.frame(class=c(rep("G1",length(Group1_sample)),rep("G2",length(Group2_sample)),rep("G3",length(Group3_sample)),rep("ALL",length(Group4_sample))),
value=c(Group1_sample,Group2_sample,Group3_sample,Group4_sample) )

Namess<-colnames(TTdata)[k]

my_comparisons <- list( c("G1", "G2"), c("G1", "G3"), c("G2", "G3"))

p1<-ggboxplot(TTdata, x = "class", y = "value"  ,size=0.7,  font.label = list(size = 4, color = "black"),repel=T,ylim=c(-0.2,0.7),
      palette="Set3",fill="class",ylab="Log2(T/N)",xlab="",outlier.shape = NA)+
         stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test",p.adjust.method="BH")

##double check p.adjust
aa<-wilcox.test(Group1_sample,Group2_sample)$p.value
bb<-wilcox.test(Group2_sample,Group3_sample)$p.value
cc<-wilcox.test(Group1_sample,Group3_sample)$p.value

FFq<-p.adjust(c(wilcox.test(Group1_sample,Group2_sample)$p.value ,
wilcox.test(Group2_sample,Group3_sample)$p.value,
wilcox.test(Group1_sample,Group3_sample)$p.value),method="BH")

FFq

pdf("Patient_Wise_Correc.pdf",height=5,width=5)
   p1
dev.off()

############ ANOVA and kruskal.test
ZZZ<-data.frame(rho=c(Group1_sample,Group2_sample,Group3_sample),
class=c(rep("A",23),rep("B",58),rep("C",7)))

kruskal.test(rho~class,data=ZZZ)
result=aov(rho~class,data=ZZZ)
summary(result)



#######################################################Gene_Wise_correlation

k<-1
Group1_sampleTT<-NULL
while(k<=length(Protein_RNA[,1]))
{

Group1_sampleTT<-c(Group1_sampleTT,cor.test(as.numeric(Protein_RNA[k,2:24]), as.numeric(Protein_RNA[k,91:113]), method = "spearm",exact=F)$ estimate)

k<-k+1
}


k<-1
Group2_sampleTT<-NULL
while(k<=length(Protein_RNA[,1]))
{

Group2_sampleTT<-c(Group2_sampleTT,cor.test(as.numeric(Protein_RNA[k,25:82]), as.numeric(Protein_RNA[k,114:171]), method = "spearm",exact=F)$ estimate)

k<-k+1
}



k<-1
Group3_sampleTT<-NULL
while(k<=length(Protein_RNA[,1]))
{

Group3_sampleTT<-c(Group3_sampleTT,cor.test(as.numeric(Protein_RNA[k,83:89]), as.numeric(Protein_RNA[k,172:178]), method = "spearm",exact=F)$ estimate)

k<-k+1
}


k<-1
Group4_sampleTT<-NULL
while(k<=length(Protein_RNA[,1]))
{

Group4_sampleTT<-c(Group4_sampleTT,cor.test(as.numeric(Protein_RNA[k,2:89]), as.numeric(Protein_RNA[k,91:178]), method = "spearm",exact=F)$ estimate)

k<-k+1
}




library(ggridges)
library(ggplot2)
theme_set(theme_minimal())


Tdata<-data.frame(rho=c(Group1_sampleTT,Group2_sampleTT,Group3_sampleTT,Group4_sampleTT), group=c(rep("A",length(Group1_sampleTT)),rep("B",length(Group2_sampleTT)),rep("C",length(Group3_sampleTT)),rep("D",length(Group4_sampleTT))) )
OutPutTT_ALLdata<-cbind(Protein_RNA[,1],cbind(Group1_sampleTT,rep("G1",length(Group1_sampleTT))), cbind(Group2_sampleTT,rep("G2",length(Group2_sampleTT)))   , cbind(Group3_sampleTT,rep("G3",length(Group3_sampleTT)))   , cbind(Group4_sampleTT,rep("ALL",length(Group4_sampleTT))))
write.table(OutPutTT_ALLdata,"Gene_wise_SpearsmanCorALL.txt",col.names=F,row.names=F,sep="\t",append=F,quote=F)

OutPutTT<-cbind(Protein_RNA[,1],Group1_sampleTT,Group2_sampleTT,Group3_sampleTT)

write.table(OutPutTT,"Gene_wise_SpearsmanCorG1G2G3.txt",col.names=F,row.names=F,sep="\t",append=F,quote=F)




Tdata<-data.frame(rho=c(Group1_sampleTT,Group2_sampleTT,Group3_sampleTT,Group4_sampleTT), group=c(rep("G1",length(Group1_sampleTT)),rep("G2",length(Group2_sampleTT)),rep("G3",length(Group3_sampleTT)),rep("ALL",length(Group4_sampleTT))) )

my_comparisons <- list( c("G1", "G2"), c("G1", "G3"), c("G2", "G3") )
p2<-ggboxplot(Tdata, x = "group", y = "rho"  ,size=0.6,  font.label = list(size = 4, color = "black"),repel=T,notch=F,ylim=c(-1,2),
        fill="group" ,ylab="Log2(T/N)",xlab="",outlier.shape = NA)+
         stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test",p.adjust.method="BH")

############ ANOVA and kruskal.test
ZZZ<-data.frame(rho=c(Group1_sampleTT,Group2_sampleTT,Group3_sampleTT),
class=c(rep("A",9009),rep("B",9009),rep("C",9009)))

kruskal.test(rho~class,data=ZZZ)
result=aov(rho~class,data=ZZZ)
summary(result)
#############


pdf("Gene_Wise_Correc.pdf",height=5,width=5)
   p2
dev.off()



#################################







