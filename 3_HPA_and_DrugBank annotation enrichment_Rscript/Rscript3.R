#Example

library(data.table)

myproExp<-as.matrix(fread("Final_New_top20_Gdata22_Anova_TurkeryHSD.txt",sep="\t",header=T))

   inFile<-list.files(pattern="Anno_")
   outFile<-paste("DE_",list.files(pattern="Anno_"),sep="")

uu<-1
while(uu<=7)
{
   
  HPAprognostic<-as.matrix(fread(inFile[uu],sep="\t",header=T))

k<-3
temp<-NULL
while(k<=length(myproExp[,1]))
{

  tttIN<-grep(myproExp[k,1],HPAprognostic[,5],fixed=T)
   if(length(tttIN)>=1)
   {
         zz<-1
        while(zz<=length(tttIN))
        {
       aaas3<-t(as.matrix(HPAprognostic[tttIN[zz],])) 
       sRRD<-t(as.matrix(myproExp[k,]))
       temp<-rbind(temp,cbind(sRRD,aaas3))
       zz<-zz+1
       }
   }

k<-k+1
}

colnames(temp)<-c(myproExp[2,],colnames(HPAprognostic))
write.table(temp,outFile[uu],sep="\t",col.names=T,row.names=F,append=F,quote=F)

uu<-uu+1

}


#######################DrugBank


library(data.table)
library(plotrix)

DrugBank<-as.matrix(fread("Anno_DrugbankData.txt",sep="\t",header=T))

DEdata<-as.matrix(fread("Final_New_top20_Gdata22_Anova_TurkeryHSD.txt",sep="\t",header=T))

k<-3
temp<-NULL
while(k<=length(DEdata[,1]))
{
  tttIN<-grep(DEdata[k,1],DrugBank[,15],fixed=T)
   tttIN<-tttIN[which(DrugBank[tttIN,5]=="yes")]
   if(length(tttIN)>=1)
   {
       Actiona<-pasteCols(as.matrix(unique(DrugBank[tttIN,4])),sep=";")
       DRGUSS<-pasteCols(as.matrix(unique(DrugBank[tttIN,6])),sep=";")
        DRGUSSPPPS<-pasteCols(as.matrix(unique(DrugBank[tttIN,5])),sep=";")
       SSSq1<-pasteCols(as.matrix(DrugBank[tttIN,8]))
       DrugCategories<-pasteCols(as.matrix(c( grepl("approved",SSSq1),grepl("investigational",SSSq1),grepl("experimental",SSSq1),grepl("nutraceutical",SSSq1),grepl("illicit",SSSq1),grepl("withdrawn",SSSq1) )),sep="\t" )

       Catergysss<-pasteCols(as.matrix(DrugBank[tttIN,9]),"__")
       temp<-rbind(temp,cbind(DEdata[k,1],Actiona,DRGUSS,DRGUSSPPPS,DrugCategories,Catergysss,t(DEdata[k,])))
   }

k<-k+1
}

colnames(temp)<-c("ID","Action","Drugs","knownaction","approved\tinvestigational\texperimental\tnutraceutical\tillicit\twithdrawn","Categories",(DEdata[2,]))

write.table(temp,"DE_Anno_DrugbankData.txt",sep="\t",col.names=T,row.names=F,append=F,quote=F)


##get Drug Categories
  ALLcLAsaa<-unique(unlist(  strsplit(temp[,6],"__")))

k<-1
class_summ<-NULL
while(k<=length(ALLcLAsaa))
{
      class_summ<-rbind(class_summ,cbind(ALLcLAsaa[k], length(grep(ALLcLAsaa[k],temp[,6],fixed=T))))
     
k<-k+1

}
   SSSS1233<- class_summ[order(as.numeric(class_summ[,2]),decreasing = T),]
   colnames(SSSS1233)<-c("Categories","Count")
#write.table(SSSS1233,"DE_DrugBankclass_counts.txt",sep="\t",col.names=T,row.names=F,append=F,quote=F)



######
CheckLOO<-c("Antineoplastic Agents",
"Protective Agents",
"Anti-Infective Agents",
"Vasodilating Agents",
"Antineoplastic and Immunomodulating Agents",
"Central Nervous System Agents",
"Immunosuppressive Agents",
"Cardiovascular Agents",
"Myelosuppressive Agents",
"Peripheral Nervous System Agents",
"Hypotensive Agents",
"QTc Prolonging Agents",
"Sensory System Agents",
"Hematologic Agents",
"Anti-Inflammatory Agents",
"Antirheumatic Agents",
"Dental Agents",
"Gastrointestinal Agents",
"Cardiotoxic antineoplastic agents",
"Antidepressive Agents",
"Anti-Bacterial Agents",
"Neurotransmitter Agents",
"Potential QTc-Prolonging Agents",
"Antiparasitic Agents",
"Miscellaneous Therapeutic Agents",
"Serotonin Agents",
"Antiplatelet agents",
"Respiratory System Agents",
"Antihypertensive Agents",
"Gastrointestinal Acidifying Agents",
"Hepatotoxic Agents",
"Highest Risk QTc-Prolonging Agents",
"Bradycardia-Causing Agents",
"Anticarcinogenic Agents",
"Growth Substances",
"BCRP/ABCG2 Inhibitors",
"Antidotes",
"Analgesics",
"Supplements")
##########


INforCCC<-temp[,6]

k<-1
temp_Catee_eff<-NULL
while(k<=length(INforCCC))
{
      i<-1
      SStor_temp<-NULL
	while(i<=length(CheckLOO))
	{
           SStor_temp<-cbind(SStor_temp,grepl(CheckLOO[i],INforCCC[k],fixed=T))
      i<-i+1
	}
   temp_Catee_eff<-rbind(temp_Catee_eff,pasteCols(as.matrix(t(SStor_temp)),sep="\t"))

k<-k+1
}

Out_ppoi<-cbind(temp,temp_Catee_eff)
colnames(Out_ppoi)[105]<- pasteCols(as.matrix(CheckLOO),sep="\t")
write.table(Out_ppoi,"DE_Anno_DrugBank_with_category.txt",sep="\t",col.names=T,row.names=F,append=F,quote=F)


#Based on above data of each categories, it can generate 2x2  contingency tables. then perfom statistical test.
##############################################  statistical test
##Cancer-related prognostic gene set enrichment analysis 
#example   Fisher's exact test
fisher_A1<-fisher.test(rbind(c(25,1171),c(142,17639)), alternative ="greater")$p.value

##Disease-related and drug target gene set enrichment analysis
#example Pearson's chi-squared 
chisq_A1<-chisq.test(rbind(c(15,186),c(735,18041)),correct=F)$p.value

##Adjust P-values for Multiple Comparisons
#example
adjustedP<-as.matrix(p.adjust(c(4.810056e-05,0.04240692,0.5892773,0.4333686,0.0001763425,0.4864509),method="fdr"))

############################



