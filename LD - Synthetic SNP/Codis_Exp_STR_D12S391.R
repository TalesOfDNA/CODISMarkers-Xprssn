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
raw1000genomesData <- read.table("STR_D12S391_1000GenomesNew.vcf", header = TRUE, skip = 252, comment.char = "") #reads vcf file
#header is to make the first line the column names, skip = 252 to not read in metadata, comment.char turns off the read of comments
#NOTE THAT IT'S NOT JUST SNP DATA HERE, INDELS ARE INCLUDED AS WELL.  A MORE CLEAR VARIABLE NAME MIGHT BE raw1000genomesData.  THEN YOU CAN JUST NOT MODIFY THAT AND YOU DON'T NEED THE ADDITIONAL VARIABLE BELOW

genotypes <- raw1000genomesData[-c(1:9)] #getting only genotypes

variantNames <- raw1000genomesData[c(3)] #isolate snp names 
#AGAIN, NOT ONLY SNPS
syntheticSnpDf <- data.frame("Synthetic") #create synthetic df
colnames(syntheticSnpDf)[1] <- "ID" #names column
variantAndSyntheticNames <- rbind(syntheticSnpDf, variantNames) #binds snp names with synthetic snp 
#VARIABLE NAME SHOULD CONVEY THAT THESE ARE ONLY NAMES, NOT GENOTYPES (LIKE SNPandSyntheticNames) [ALSO, AGAIN, NOT ONLY SNPS]

position <- raw1000genomesData[c(2)] #isolate snp positions
syntheticPosition <- data.frame(12449954) #creates synthetic position, at the lowest to make high cor obvious
colnames(syntheticPosition)[1] <- "Synthetic" #names column
positionAndSynthetic <- rbind(as.integer(syntheticPosition), position) #binds snp positions with synthetic snp

genotypesTranspose<- as.data.frame(t(genotypes)) #turn genotype data into a dataframe and transpose
sampleNames <- row.names(genotypesTranspose) #create "rowNames" variable
#A MORE DESCRIPTIVE VARIABLE NAME WILL MAKE YOUR CODE MORE CLEAR.  WHAT ROW NAMES?  YOU COULD SAY sampleNames OR indivNames INSTEAD
row.names(genotypesTranspose) <- NULL #make row names NULL
variantNamesAndGenotypes <- cbind(sampleNames, genotypesTranspose) #combine the names and genotypes into one dataframe
#I KNOW I'M HARPING ON VARIABLE NAMES, BUT I THINK BY NAMING THEM MORE CLEARLY, YOU'LL HAVE AN EVEN STRONGER HANDLE ON YOUR CODE.  FIRST THING, THEY'RE NOT ONLY SNPS.  SECOND, CALLING SOMETHING 'FINAL' IS DANGEROUS BECAUSE YOU OFTEN END UP MODIFYING IT AGAIN IN WAYS YOU DIDN'T ANTICIPATE.  TRY TO USE A MORE DESCRIPTIVE TERM FOR EXACTLY HOW YOU HAVE THE DATA SET UP IN THIS VARIABLE
colnames(variantNamesAndGenotypes)[1] <- "Sample ID" #rename "names" as "sample ID" leaving a data frame of
#the names of all individuals and their genotypes

#STR DOSAGES DATA
library(readr)
strDosageRaw <- read_table2("~/Documents/GitHub/CODISMarkers-Xprssn/Data/D12S391_dose_exp.txt") #import updated average STR dosage data
#strDosageRaw <- read.table("~/Documents/GitHub/CODISMarkers-Xprssn/Data/D12S391_dose_exp.txt", header=TRUE, row.names=1) #import updated average STR dosage data
#USING THE VERY SIMILAR FUNCTION READ.TABLE, YOU CAN ACCOUNT FOR BOTH ROW AND COLUMN NAMES, SO YOU DON'T GET THAT WARNING AND PROPERLY READ IN THE WHOLE FILE.  I'LL COMMENT OUT MY VERSION FOR NOW SINCE YOUR CODE DEPENDS ON GETTING IT THE WAY YOU DO

strDosage <- strDosageRaw[c(1:2)] #remove gene expression

SyntheticDosageNumbers <- strDosage[c(2)]
colnames(SyntheticDosageNumbers)[1] <- "Synthetic"

colnames(strDosage)[1] <- "Sample ID" #rename first column 
colnames(strDosage)[2] <- "STR Dosage" ##rename second column

strDosageAndSynthetic <- cbind(strDosage, SyntheticDosageNumbers)
#MORE DESCRIPTIVE/SPECIFIC VARIABLE NAME

plot(strDosageAndSynthetic$Synthetic, strDosageAndSynthetic$`STR Dosage`) #plot against each other to see if they are EXACTLY correlated
sum(is.na(SyntheticDosageNumbers)) #count the total number of NAs; total = 0

for (i in 1:nrow(strDosageAndSynthetic)){
  for (j in 3:ncol(strDosageAndSynthetic)){
    strDosageAndSynthetic[i,j]<-syntheticSnpFunction(strDosageAndSynthetic[i,j]) #applies function to df
  }
}

idDosageAndGenoWithSynthetic <- merge(strDosageAndSynthetic, variantNamesAndGenotypes, by = "Sample ID") #combine genotypes and STR dosages by matching up Sample IDs 
#and disregarding nonmatches

bdspData <- read.table("D12S391_100kb.vcf", header = FALSE) #get STR data from old data 
#WHAT OLD DATA?  MORE SPECIFIC NAMING
strPosition <- subset(bdspData, V3 == "D12S391") #subset STR data to find the position of the STR; == 12449954

idDosageAndGenoWithSyntheticM <- as.matrix(idDosageAndGenoWithSynthetic) #transforming combo into matrix 
#TRANSFORMING YOUR DATA BETWEEN DATAFRAMES AND MATRICES GIVES YOU A LOT OF SIMILAR, BUT DISTINCT VARIABLES.  MIGHT BE USEFUL TO REDUCE THE NUMBER OF DATA TYPE CONVERSIONS IF POSSIBLE
#yes, it is necessary to transform into a matrix bc otherwise an NA value is generated for all genotypes instead. It is also much faster. 

#applying the function to the new matrix to convert gentotypes into values
for (i in 1:nrow(idDosageAndGenoWithSyntheticM)){
  for (j in 4:ncol(idDosageAndGenoWithSyntheticM)){
    idDosageAndGenoWithSyntheticM[i,j]<-changeDiscrete(idDosageAndGenoWithSyntheticM[i,j])
  }
}

idDosageAndGenoWithSyntheticDf<-as.data.frame(idDosageAndGenoWithSyntheticM) #putting  transformed matrix into new dataframe
str(idDosageAndGenoWithSyntheticDf) #check structure of data

#make elements in STR dosage and genotype columns into numeric
idDosageAndGenoWithSyntheticDf[-c(1)] <- lapply(idDosageAndGenoWithSyntheticDf[-c(1)], function(x) {  #remove the sample IDs
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
str(idDosageAndGenoWithSyntheticDf[-c(1)]) #check structure of new data
#3312
cors <- rep(NA, 3311) #setting an empty list to put calculated correlation later on

#use for loop to calculate correlation
for (i in 3:3313){
  cors[i-2]<- cor(idDosageAndGenoWithSyntheticDf[,2], idDosageAndGenoWithSyntheticDf[,i])
}

#cor(idDosageAndGenoWithSyntheticDf[,2],idDosageAndGenoWithSyntheticDf[,3]) to double check correlation
#cor(idDosageAndGenoWithSyntheticDf[,2],idDosageAndGenoWithSyntheticDf[,7])
#cor(idDosageAndGenoWithSyntheticDf[,2],idDosageAndGenoWithSyntheticDf[,3311])

#creating the new dataframe with first row as SNP name
variantAndSyntheticPosNamesCor<- cbind(positionAndSynthetic, variantAndSyntheticNames, cors) #combine the position of the SNP, name of SNP, 
#and it's correlation with the STR dosage into a single dataframe
#POSITIONANDSYNTHETIC IS ACTUALLY THE POSITIONS INCLUDING SYNTHETIC, RIGHT?  SNPANDSYNTHETIC IS  THE BI-ALLELIC VARIANT _NAMES_ INLUDING SYNTHETIC, RIGHT?  MORE ACCURATE VARIABLE NAMES WILL HELP MAKE IT EASIER TO UNDERSTAND WHAT'S BEING ACCOMPLISHED IN THIS LINE.  I'D ALSO RECONSIDER THE NAMING OF THIS VARIABLE
#IS IT AT ALL POSSIBLE THAT THE VARIANT ORDER ISN'T THE SAME BETWEEN POSITIONANDSYNTHETIC, SNPANDSYNTHETIC, AND CORRELATION?  I DON'T THINK SO, BUT THE BOOK-KEEPING MAKES IT A BIG OPAQUE

sum(is.na(variantAndSyntheticPosNamesCor$cors)) #count the total number of NAs; total = 1903
naOmitVariantAndSyntheticPosNamesCor <- na.omit(variantAndSyntheticPosNamesCor)

#plot the data
pdf(file="STR_D12S391_corplotSyntheticUpdate.pdf", width=50, height=20) #create a pdf of the plot
plot(naOmitVariantAndSyntheticPosNamesCor$POS, naOmitVariantAndSyntheticPosNamesCor$cors, main = "Correlation with Synthetic SNP", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()


#plot within 100kb window
hundredKbSynth <- subset(naOmitVariantAndSyntheticPosNamesCor, POS < 12499954 & POS > 12399954, select = c(POS, ID, cors))
#subset positions within 100kb of the STR and their correlation

pdf(file="STR_D12S391_corplot100kbSytheticUpdate.pdf", width=20, height=10) #create a pdf of the plot
plot(hundredKbSynth$POS, hundredKbSynth$cors, main = "Correlation with Synthetic SNP Within 100 Kb", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()

#Plot without Synthetic SNP
idDosageAndGenoNoSynthetic <- merge(strDosage, variantNamesAndGenotypes, by = "Sample ID")

idDosageAndGenoNoSyntheticM <- as.matrix(idDosageAndGenoNoSynthetic) #transforming combo into matrix 

#applying the function to the new matrix to convert gentotypes into values
for (i in 1:nrow(idDosageAndGenoNoSyntheticM)){
  for (j in 3:ncol(idDosageAndGenoNoSyntheticM)){
    idDosageAndGenoNoSyntheticM[i,j]<-changeDiscrete(idDosageAndGenoNoSyntheticM[i,j])
  }
}

idDosageAndGenoNoSyntheticDf<-as.data.frame(idDosageAndGenoNoSyntheticM) #putting  transformed matrix into new dataframe

idDosageAndGenoNoSyntheticDf[-c(1)] <- lapply(idDosageAndGenoNoSyntheticDf[-c(1)], function(x) { #remove the sample IDs
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

cors <- rep(NA, 3310) #setting an empty list to put calculated correlation later on

#use for loop to calculate correlation
for (i in 3:3312){
  cors[i-2]<- cor(idDosageAndGenoNoSyntheticDf[,2], idDosageAndGenoNoSyntheticDf[,i])
}

#cor(idDosageAndGenoNoSyntheticDf[,2],idDosageAndGenoNoSyntheticDf[,3]) to double check correlation
#cor(idDosageAndGenoNoSyntheticDf[,2],idDosageAndGenoNoSyntheticDf[,3310])
#cor(idDosageAndGenoNoSyntheticDf[,2],idDosageAndGenoNoSyntheticDf[,3312])

variantsPosNamesCor <- cbind(position, variantNames, cors) #combine the position of the SNP, name of SNP, 

sum(is.na(variantsPosNamesCor$cors)) #count the total number of NAs; total = 1903 
naOmitvariantsPosNamesCor<- na.omit(variantsPosNamesCor)

#plot the data
pdf(file="STR_D12S391_corplotUpdate.pdf", width=50, height=20) #create a pdf of the plot
plot(naOmitvariantsPosNamesCor$POS, naOmitvariantsPosNamesCor$cors, main = "Correlation without Synthetic SNP", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()


#plot within 100kb window
hundredKbNoSynth <- subset(naOmitvariantsPosNamesCor, POS < 12499954 & POS > 12399954, select = c(POS, ID, cors))
#subset positions within 100kb of the STR and their correlation

pdf(file="STR_D12S391_corplot100kbUpdate.pdf", width=20, height=10) #create a pdf of the plot
plot(hundredKbNoSynth$POS, hundredKbNoSynth$cors, main = "Correlation without Synthetic SNP within 100 Kb", xlab= "Position", ylab= "Correlation")#plot by correlation vs position of the SNP 
abline(v = 12449954, col = "red" )
dev.off()
 

