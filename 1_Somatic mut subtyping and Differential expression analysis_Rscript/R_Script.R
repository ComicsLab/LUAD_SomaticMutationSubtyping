#Example

library(data.table)
library(cluster)
library(ConsensusClusterPlus)

A_data<-fread("All_data_TopALLNew_data.txt",sep="\t",header=T)

B_data<-t(A_data)

B_data<-data.frame(B_data[2:length(B_data[,1]),1:length(B_data[1,])])
gower.dissimilarity.mtrx <- daisy(B_data, metric = c("gower"))


###################a ¡§decrease progressively¡¨ manner to discard unnecessary somatic mutation genes
#ALL_data
seed=12345
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="ALL_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top5000
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:5000])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top5000_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)



#top2500
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:2500])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top2500_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)



#top1200
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:1200])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top1200_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top600
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:600])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top600_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top300
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:300])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top300_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top150
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:150])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top150_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)

#top100
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:100])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top100_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top50
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:50])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top50_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top25
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:25])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top25_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top20
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:20])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top20_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)

#top15
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:15])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top15_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)


#top10
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:10])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top10_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)

#top5
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:5])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top5_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)

#top3
B_data1<-data.frame(B_data[1:length(B_data[,1]),1:3])
gower.dissimilarity.mtrx <- daisy(B_data1, metric = c("gower"))
dist<- as.matrix(gower.dissimilarity.mtrx)
res <- ConsensusClusterPlus(as.matrix(dist), maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1, title="Top3_NEwRank_SOMATIC_consensus_cluster",clusterAlg = "pam", plot="png", writeTable=T,seed=seed)



############## select TOP20 cluster result (provides stable and optimal efficiency of clustering)
A_data<-as.matrix(fread("All_data_TopALLNew_data.txt",sep="\t",header=T))

Top20Gene<-as.matrix(A_data[1:20,1])

#set file path 
classNeed<-as.matrix(fread("<--setPath-->/Top20_NEwRank_SOMATIC_consensus_cluster/Top20_NEwRank_SOMATIC_consensus_cluster.k=3.consensusClass.csv",sep=",",header=F))
#classNeed<-as.matrix(fread("Top20_NEwRank_SOMATIC_consensus_cluster.k=3.consensusClass.csv",sep=",",header=F))


k<-1
tempDaT<-NULL
tempDaVT<-NULL
while(k<=max(as.numeric(classNeed[,2])))
{
DaP<-grep(as.character(k),classNeed[,2],fixed=T)

clasNN<-rep(as.character(k),length(DaP))

tempDaT<-c(tempDaT,clasNN)
tempDaVT<-c(tempDaVT,DaP)

k<-k+1
}

  ERWE<-as.matrix(rbind(tempDaT,A_data[1:20,(tempDaVT)+1]))
  ERWE1<-cbind(as.matrix(c("Class",Top20Gene)),ERWE)
AQWout<-rbind(t(as.matrix(c("",colnames(A_data)[(tempDaVT)+1]) )), ERWE1)

#write.table(AQWout,"C3_NEWtop20_result.txt",sep="\t",col.names=F,row.names=F,append=F,quote=F)

dataLilst<-colnames(  ERWE)

clinicalDARA<-t(as.matrix(fread("Clinical_DARA.txt",sep="\t",header=F)))

k<-1
temp<-NULL
while(k<=length(dataLilst))
{

temp<-rbind(temp,clinicalDARA[ ,grep(dataLilst[k] , clinicalDARA[1,],fixed=T)] )


k<-k+1
}

Clinal_RR<-rbind(AQWout,cbind(clinicalDARA[,1],(t(temp))))
#output_class_clinicallabel_protein expression
write.table( Clinal_RR,"C3_newtop20_Clicalresult.txt",sep="\t",col.names=F,row.names=F,append=F,quote=F)





##################GetProteinexpression and perform Quantile Normalization
library(preprocessCore)
ProteinR_data<-(fread("ProteinRatio.txt",sep="\t",header=T))


nordata<-as.matrix(ProteinR_data[,2:length(ProteinR_data[1,])])
nordataNN<-NULL


X2<-normalize.quantiles(as.matrix(nordata))
ProteinR_data<-as.matrix(ProteinR_data)
ProteinR_data[,2:length(ProteinR_data[1,])]<-as.matrix(X2)

k<-1
temp<-NULL
while(k<=length(dataLilst))
{
temp<-cbind(temp,(as.matrix(ProteinR_data[ ,match(dataLilst[k] , colnames(ProteinR_data))])) )
k<-k+1
}

write.table(rbind(Clinal_RR,cbind(ProteinR_data[,1] ,temp)),"C3_newtop20_RatioProteinresult.txt",sep="\t",col.names=F,row.names=F,append=F,quote=F)


###############################################################start get differential protein of 3 class
 
library(data.table)

AllDATA<-as.matrix(fread("C3_newtop20_RatioProteinresult.txt",sep="\t",header=F))



tempPP<-NULL
tempPPdunnTest<-NULL

#class
 GroupQ=c(AllDATA[2,2:24],AllDATA[2,25:82],AllDATA[2,83:89])

k<-35
while(k<=length(AllDATA[,1]))
{
 ValueQ=as.numeric((c(AllDATA[k,2:24],AllDATA[k,25:82],AllDATA[k,83:89])))
SweetTest=data.frame(ValueQ,GroupQ)

#class region
if( length(which(is.na(AllDATA[k,2:24])))/length(AllDATA[k,2:24])<0.5 &&
 length(which(is.na( AllDATA[k,25:82])))/length( AllDATA[k,25:82])<0.5 &&
 length(which(is.na(AllDATA[k,83:89])))/length(AllDATA[k,83:89])<0.5)
{


result=aov(ValueQ~GroupQ,data=SweetTest)
SS<-summary(result)
KP_R<-unlist(SS)[9]
T_KP_R<-KP_R

kruskal_KP_R<-kruskal.test(ValueQ~GroupQ,data=SweetTest)

	tempPP<-rbind(tempPP,cbind(k,T_KP_R,kruskal_KP_R$ p.value))


}

k<-k+1
}

WWEWET<-cbind(as.matrix(AllDATA[tempPP[,1],]),tempPP)


#write.table(cbind(as.matrix(AllDATA[tempPP[,1],1]),tempPP),"ST_newtop20_RResult.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)
write.table(WWEWET,"ST_newtop20_ProRResult.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)

#output: for PermutationFDR raw data
write.table(rbind( as.matrix(AllDATA[1:2,])  , WWEWET[,1:(length(WWEWET[1,])-3)]),"For_PermutationFDR_ST_newtop20_ProRResult.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)



##########################################PermutationFDR
#function
PermutationTest<-function(AllDATA,ClassLabelNow)
{
tempPP<-NULL
tempPPdunnTest<-NULL

u<-3
end<-length(AllDATA[,1])

for(k in u:end)
{
 ValueQ=as.numeric(AllDATA[k,])
SweetTest=data.frame(ValueQ,ClassLabelNow)
result=aov(ValueQ~ClassLabelNow,data=SweetTest)
SS<-summary(result)
KP_R<-unlist(SS)[9]
T_KP_R<-KP_R
kruskal_KP_R<-kruskal.test(ValueQ~ClassLabelNow,data=SweetTest)
tempPP<-rbind(tempPP,cbind(T_KP_R,kruskal_KP_R$ p.value))

}

return(tempPP)
}


######################### run random pair data for Permutation FDR
library(data.table)
set.seed(12345)
Data<-as.matrix(fread("For_PermutationFDR_ST_newtop20_ProRResult.txt",sep="\t",header=F))

ClassLabelNow<-Data[2,2:length(Data[1,])]
ExpressionNow<-Data[ 3:length(Data[,1]),2:length(Data[1,])  ]
XqClass<-table(ClassLabelNow)

p_Index<-1
P_store<-NULL
while(p_Index<=300) # Permutation count 300times
{
poolnow<-1:length(ClassLabelNow)
poolremain<-poolnow

k<-1
randomClass<-NULL
while(k<=length(XqClass))
{
 getN<-XqClass[k]
 getSam<-sample(poolremain,getN)
randomClass<-c(randomClass,getSam)
poolremain<- setdiff(poolremain,getSam)
k<-k+1
}


RandomDataA<-ExpressionNow[,sample(sample(poolnow))]
ZZZ<-PermutationTest(RandomDataA,ClassLabelNow)
P_store<-rbind(P_store,ZZZ)
p_Index<-p_Index+1
cat(p_Index)
}

##output ramdom class test p.value results
write.table(P_store,"NEW_top20_data_random_Pair.txt",col.names=F,row.names=F,quote=F,append=F)



RanmdomDataP<-as.matrix(fread("NEW_top20_data_random_Pair.txt",sep=" ",header=F))

###ANOVA.p
   TD_all<-c(as.numeric(tempPP[,2]),as.numeric(RanmdomDataP[,1]))
   TD_tag<- c(rep(0,length(as.numeric(tempPP[,2]))), rep(1,length(as.numeric(RanmdomDataP[,1]))))
     TD_all_order<-TD_all[order(TD_all)]
     TD_tag_order<-TD_tag[order(TD_all)]

k<-0
permutationN<-300
Temp_ANO<-NULL
while(k<=0.1)
{
    indexN<-which(TD_all_order<=k)
 
    Temp_ANO<-rbind(Temp_ANO,cbind(k,sum(TD_tag_order[indexN])/ ((length(indexN)-sum(TD_tag_order[indexN]))*permutationN),(length(indexN)-sum(TD_tag_order[indexN]))))
  
k<-k+0.0001

}


###KruskalWallis.p
   TD_all<-c(as.numeric(tempPP[,3]),as.numeric(RanmdomDataP[,2]))
   TD_tag<- c(rep(0,length(as.numeric(tempPP[,3]))), rep(1,length(as.numeric(RanmdomDataP[,2]))))
     TD_all_order<-TD_all[order(TD_all)]
     TD_tag_order<-TD_tag[order(TD_all)]

k<-0
permutationN<-300
Temp_KRU<-NULL
while(k<=0.1)
{
    indexN<-which(TD_all_order<=k)
 
    Temp_KRU<-rbind(Temp_KRU,cbind(k,sum(TD_tag_order[indexN])/ ((length(indexN)-sum(TD_tag_order[indexN]))*permutationN),(length(indexN)-sum(TD_tag_order[indexN]))))
  
k<-k+0.0001

}

write.table(Temp_ANO,"FDR_newtop20_Table_ANO.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)
#write.table( Temp_KRU,"FDR_newtop20_Table_Temp_KRU.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)


#passFDR 0.1

ProRResult<-as.matrix(fread("ST_newtop20_ProRResult.txt",sep="\t",header=F))

#ANOVA
CutOffYouNeed<-0.01135 #based on Temp_ANO data
PassFDR_indexa<-which(as.numeric(ProRResult[,91])< CutOffYouNeed)
#Outpu Final result Anova, FDR <0.1
write.table(rbind(cbind(as.matrix(AllDATA[1:2,]),"k","Anova.p","KruskalWallis.p"),ProRResult[PassFDR_indexa,]),"Top20_PASS_FDRcriteria_data_ANOVA.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)

#KRUSTA
CutOffYouNeed<-0.00151 #based on  Temp_KRU data
PassFDR_indexb<-which(as.numeric(ProRResult[,92])<= CutOffYouNeed)
#write.table(rbind(cbind(as.matrix(AllDATA[1:2,]),"k","Anova.p","KruskalWallis.p"),ProRResult[PassFDR_indexb,]),"Top20_PASS_FDRcriteria_data_KruskalWallis.txt",sep="\t",row.names=F,col.names=F,quote=F,append=F)


#############################perform ANOVA + TukeyHSD


library(data.table)

AllDATA<-as.matrix(fread("Top20_PASS_FDRcriteria_data_ANOVA.txt",sep="\t",header=F))

tempPP<-NULL
tempPPdunnTest<-NULL


  GroupQ=(AllDATA[2,2:(length(AllDATA[2,])-3)])

k<-3
while(k<=length(AllDATA[,1]))
{

ValueQ=as.numeric(AllDATA[k,2:(length(AllDATA[2,])-3)])
SweetTest=data.frame(ValueQ,GroupQ)

result=aov(ValueQ~GroupQ,data=SweetTest)
SS<-summary(result)
KP_R<-unlist(SS)[9]
T_KP_R<-KP_R
TukeyHSD<-TukeyHSD(result,conf.level=0.95)
	tempPP<-rbind(tempPP,cbind(k,as.matrix(T_KP_R),t(as.matrix(TukeyHSD$ GroupQ[,4]))))

k<-k+1
}

#write.table(cbind(as.matrix(AllDATA[3:length(AllDATA[,1]),1]),tempPP),"New_top20_Anova_TurkeryHSD.txt",col.names=T,row.names=F,append=F,quote=F,sep="\t")



Gdata<-cbind(as.matrix(AllDATA[3:length(AllDATA[,1]),1]),tempPP)

Gdata22<-Gdata

 Gdata22[which(as.numeric(Gdata[,4])<0.05),4]<-"A"
 Gdata22[which(as.numeric(Gdata[,5])<0.05),5]<-"B"
 Gdata22[which(as.numeric(Gdata[,6])<0.05),6]<-"C"
 Gdata22[which(as.numeric(Gdata[,4])>=0.05),4]<-"0"
 Gdata22[which(as.numeric(Gdata[,5])>=0.05),5]<-"0"
 Gdata22[which(as.numeric(Gdata[,6])>=0.05),6]<-"0"

QQQQc<-rbind(c("ID","k","ANOVA.p","G1_vs_G2","G1_vs_G3","G2_vs_G3"),rep("x",length(Gdata22[1,])),Gdata22)

write.table(cbind(AllDATA,QQQQc),"Final_New_top20_Gdata22_Anova_TurkeryHSD.txt",col.names=T,row.names=F,append=F,quote=F,sep="\t")



##########################################FDR versus DE gene/proteins graph


#FDR_newtop20_Table_ANO.txt
SomaticNumBinary<-read.table("FDR_newtop20_Table_ANO.txt",sep="\t",header=F)


dpi=600   #pixels per square inch
pdf("New_top_20_outputsGeneBased_0909_A.pdf", width=6, height=6)


plot(SomaticNumBinary[,2],SomaticNumBinary[,3],type="l",pch=".",xlim=c(0,1),ylab="# of DE proteins",xlab="Permutated FDR",main="Somatic mutation subtypes")
abline(v =0.1,col="red",lty=14)

text(0.15, 3000,"0.1",col="red")
text(0.70, 1600,"g1: EGFR and TP53 mutations",col="black" ,cex=1.15)
text(0.70, 1200,"g2: EGFR mutation only",col="black",cex=1.15)
text(0.70, 800,"g3: Multiple gene mutations",col="black",cex=1.15)

dev.off()



