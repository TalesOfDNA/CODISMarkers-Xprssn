#FUNCTIONS
changeDiscrete<- function(x){ # function that gives specific genotypes a numeric value
  if (x=="0|0"){
    x=0
  }
  else if (x== "1|1"){
    x = 2
  }
  else{
    x= 1
  }
  return(x)
}

syntheticSnpFunction <- function(x){ # function that gives specific str dosages a numeric value
  if (x >= 16 & x < 22){ #gives a 0 value for str dosages between 16-22
    x = 0
  }
  else if (x >= 22 & x < 28){ #gives a 1 value for str dosages between 22-28
    x = 1
  }
  else if (x >= 28 & x <= 34){ #gives a 2 value for str dosages between 28-34
    x = 2
  }
  else{ #returns NA should value be anything else 
    x= NA
  }
  return(x)
}

#SNP DATA
snpData <- read.table("STR_D12S391_1000Genomes.vcf", header = TRUE, skip = 252, comment.char = "") #reads vcf file
#header is to make the first line the column names, skip = 252 to not read in metadata, comment.char turns off the read of comments
rawSnpData <- snpData #save unedited SNP data  

genotypes <- snpData[-c(1:9)] #getting only genotypes

snpNames <- snpData[c(3)] #isolate snp names 
syntheticSnpDf <- data.frame("Synthetic") #create synthetic df
colnames(syntheticSnpDf)[1] <- "ID" #names column
snpAndSynthetic <- rbind(syntheticSnpDf,snpNames) #binds snp names with synthetic snp

position <- snpData[c(2)] #isolate snp positions
syntheticPosition <- data.frame(12449954) #creates synthetic position, at the lowest to make high cor obvious
colnames(syntheticPosition)[1] <- "Synthetic" #names column
positionAndSynthetic <- rbind(as.integer(syntheticPosition), position) #binds snp positions with synthetic snp

snpDataTranspose<- as.data.frame(t(genotypes)) #turn genotype data into a dataframe and transpose
rowNames <- row.names(snpDataTranspose) #create "rowNames" variable
row.names(snpDataTranspose) <- NULL #make row names NULL
snpDataFinal <- cbind(rowNames, snpDataTranspose) #combine the names and genotypes into one dataframe
colnames(snpDataFinal)[1] <- "Sample ID" #rename "names" as "sample ID" leaving a data frame of
#the names of all individuals and their genotypes

#STR DOSAGES DATA
library(readr)
strDosageRaw <- read_table2("D12S391_dose_exp.txt") #import updated average STR dosage data

strDosage <- strDosageRaw[c(1:2)] #remove gene expression

SyntheticDosageNumbers <- strDosage[c(2)]
colnames(SyntheticDosageNumbers)[1] <- "Synthetic"

colnames(strDosage)[1] <- "Sample ID" #rename first column 
colnames(strDosage)[2] <- "STR Dosage" ##rename second column

bind <- cbind(strDosage, SyntheticDosageNumbers)

plot(SyntheticDosageNumbers$Synthetic, bind$`STR Dosage`) #plot against each other to see if they are EXACTLY correlated
sum(is.na(SyntheticDosageNumbers)) #count the total number of NAs; total = 0

for (i in 1:nrow(bind)){
  for (j in 3:ncol(bind)){
    bind[i,j]<-syntheticSnpFunction(bind[i,j]) #applies function to df
  }
}

strSnpCombined <- merge(bind, snpDataFinal, by = "Sample ID") #combine genotypes and STR dosages by matching up Sample IDs 
#and disregarding nonmatches

oldData <- read.table("D12S391_100kb.vcf", header = FALSE) #get STR data from old data 
strPosition <- subset(oldData, V3 == "D12S391") #subset STR data to find the position of the STR; == 12449954

strSnpCombinedMRaw <- as.matrix(strSnpCombined) #transforming combo into matrix 
strSnpCombinedM <- strSnpCombinedMRaw #assigning new matrix so that old data can be checked if needed

#applying the function to the new matrix to convert gentotypes into values
for (i in 1:nrow(strSnpCombinedM)){
  for (j in 4:ncol(strSnpCombinedM)){
    strSnpCombinedM[i,j]<-changeDiscrete(strSnpCombinedM[i,j])
  }
}

strSnpCombinedDf<-as.data.frame(strSnpCombinedM) #putting  transformed matrix into new dataframe
str(strSnpCombinedDf) #check structure of data
DosageAndGenotypes<- strSnpCombinedDf[-c(1)] #remove the sample IDs
str(DosageAndGenotypes) #check structure of data


#make elements in STR dosage and genotype columns into numeric
DosageAndGenotypes[] <- lapply(DosageAndGenotypes, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
str(DosageAndGenotypes) #check structure of new data

cors <- rep(NA, 20156) #setting an empty list to put calculated correlation later on

#use for loop to calculate correlation
for (i in 2:20157){
  cors[i-1]<- cor(DosageAndGenotypes[,1],DosageAndGenotypes[,i])
}
correlation <- cors #turn cors into a variable

cor(DosageAndGenotypes[,1],DosageAndGenotypes[,2])
cor(DosageAndGenotypes[,1],DosageAndGenotypes[,20157])
cor(DosageAndGenotypes[,1],DosageAndGenotypes[,10522])

#creating the new dataframe with first row as SNP name
finalCorRaw <- cbind(positionAndSynthetic, snpAndSynthetic, correlation) #combine the position of the SNP, name of SNP, 
#and it's correlation with the STR dosage into a single dataframe

sum(is.na(finalCorRaw$correlation)) #count the total number of NAs; total = 7509 
finalCor <- na.omit(finalCorRaw)

#plot the data
pdf(file="STR_D12S391_corplotSynthetic.pdf", width=50, height=20) #create a pdf of the plot
plot(finalCor$POS, finalCor$correlation, xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()

max(finalCor$correlation)

#plot within 100kb window
hundredKb <- subset(finalCor, POS < 12499954 & POS > 12399954, select = c(POS, ID, correlation))
#subset positions within 100kb of the STR and their correlation

pdf(file="STR_D12S391_corplot100kbSythetic.pdf", width=20, height=10) #create a pdf of the plot
plot(hundredKb$POS, hundredKb$correlation, main = "Correlation with Synthetic SNP", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()

#Plot without Synthetic SNP
combined <- merge(strDosage, snpDataFinal, by = "Sample ID")

combinedMRaw <- as.matrix(combined) #transforming combo into matrix 
combinedM <- combinedMRaw  #assigning new matrix so that old data can be checked if needed

#applying the function to the new matrix to convert gentotypes into values
for (i in 1:nrow(combinedM)){
  for (j in 3:ncol(combinedM)){
    combinedM[i,j]<-changeDiscrete(combinedM[i,j])
  }
}

combinedMDf<-as.data.frame(combinedM) #putting  transformed matrix into new dataframe
str(combinedMDf) #check structure of data
combinedNoID<- combinedMDf[-c(1)] #remove the sample IDs
str(combinedNoID) #check structure of data

combinedNoID[] <- lapply(combinedNoID, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
str(combinedNoID)

corsFinal <- rep(NA, 20155)

for (i in 2:20156){
  corsFinal[i-1]<- cor(combinedNoID[,1],combinedNoID[,i])
}
correlationFinal <- corsFinal

cor(combinedNoID[,1],combinedNoID[,2])
cor(combinedNoID[,1],combinedNoID[,20156])
cor(combinedNoID[,1],combinedNoID[,20061])

CorRaw <- cbind(position, snpNames, correlationFinal) #combine the position of the SNP, name of SNP, 

sum(is.na(CorRaw$correlationFinal)) #count the total number of NAs; total = 7509 
Cor <- na.omit(CorRaw)

#plot the data
pdf(file="STR_D12S391_corplot.pdf", width=50, height=20) #create a pdf of the plot
plot(Cor$POS, Cor$correlationFinal, xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()

max(Cor$correlationFinal)

#plot within 100kb window
hundredKbFinal <- subset(Cor, POS < 12499954 & POS > 12399954, select = c(POS, ID, correlationFinal))
#subset positions within 100kb of the STR and their correlation

pdf(file="STR_D12S391_corplot100kb.pdf", width=20, height=10) #create a pdf of the plot
plot(hundredKbFinal$POS, hundredKbFinal$correlationFinal, main = "Correlation without Synthetic SNP", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()


