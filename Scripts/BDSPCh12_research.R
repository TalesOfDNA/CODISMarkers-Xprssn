#functions
change_discrete<- function(x){
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


num_split<- function(x){
  y<-strsplit(x, "[|]")
  return(y)
}


#Import

ch12<- readLines("D12S391_100kb.vcf")
ch12_data<- read.table("D12S391_100kb.vcf", stringsAsFactors = FALSE)


#Subsetting 2 chromosome dataset into snip and STR
ch12_snip<- subset(ch12_data, V4== "A" | V4== "C"|V4== "T"| V4 =="G")
location<-  ch12_snip[c(2)]

#subtracting snip data from original data
library(dplyr)
ch12_str<- ch12_data[ !(ch12_data$V4 %in%  ch12_snip$V4), ]

#only taking snp name and genotype
df12<- ch12_snip[-c(1:2, 4:9)]

#transforming df2 into matrix and use substring to take only first three characters 
#for the genotype
m12<- as.matrix(df12)
for (i in 2:ncol(m12)){
  m12[,i]<- substr(m12[,i], start=1, stop=3)
}
#assigning new matrix so that old data can be checked
m212<- m12

#applying the function to the new matrix to count the ones
for (i in 1:nrow(m212)){
  for (j in 2:ncol(m212)){
    m212[i,j]<-change_discrete(m212[i,j])
  }
}
#putting  transformed matrix into new dataframe
s12<-as.data.frame(m212)
str(s12)
#using loop to make elements in genotype columns into numeric
for (i in 2:ncol(s12)){
  s12[,i]<-as.numeric(as.character(s12[,i]))
}
str(s12)
#creating an empty vector to store the frequency proportion of 1
freq_prop<-c()
#calculating the row sums and divide by 2 tims the number of individuals
freq_prop<- rowSums(s12[,2:ncol(s12)])/((ncol(s12)-1)*2)
summary(freq_prop)
#putting the freq_prop vector into the exisiting s2 dataframe
s12$freq_prop<-freq_prop
#visuals
library(ggplot2)
ggplot(aes(x=freq_prop), data= s12)+geom_histogram(color= "black", fill= "purple")+
  scale_x_continuous()+labs(title = "Chromosome 12", x= "frequency",y= "count")+
  theme(plot.title = element_text(hjust = 0.5))

#subsetting D2S441 from STR data
mstr12<- subset(ch12_str, V3 == "D12S391")
#Getting only STR name and genotype
mstr12<- mstr12[-c(1:2, 4:9)]
mmstr12<- as.matrix(mstr12)

#getting only numeric part of the STR genotype
for (i in 2:ncol(mmstr12)){
  mmstr12[,i]<-  gsub(":.*","",mmstr12[,i])
}
str<- mmstr12
dim(mmstr12)

#creating the loop to calculate average
for ( i in 2:ncol(mmstr12)){
  b<-num_split(mmstr12[1,i])
  c<-matrix(unlist(b))
  df<- as.data.frame(c)
  df$V1<- as.numeric(as.character(df$V1))
  avg<-( df[1,1]+df[2,1])/2
  mmstr12[1,i]<- avg
}

#changing the elements into numeric
dstr12<- as.data.frame(mmstr12)
for (i in 2:ncol(dstr12)){
  dstr12[,i]<-as.numeric(as.character(dstr12[,i]))
}

#removing the freq proportion from the latest SNP file
new_snp<- s12[-c(2506)]

#creating the new file by combining latest SNP and STR files
combo<- rbind( dstr12, new_snp)

#removing the names columns to calculate the correlation
new_combo<- combo[-c(1)]
cm<- as.matrix(new_combo)

#setting an empty list to put the row wise correlation
cors<- rep(NA, 1254)
#using for loop to calculate the row wise correlation
for (i in 2:1255){
  cors[i-1]<- cor(cm[1,], cm[i,])
}
#creating the new dataframe with first column as SNP and STR names
latest <- ch12_snip[c(2)]
which(is.na(cors))

#putting the correlation vector into new dataframe
latest$cors<- cors

#ignore the NA's
which(is.na(cors))

#picture
pdf(file="oldCorplotCh12.pdf", width=10, height=3)
plot(latest$V2, latest$cors, xlab= "location", ylab= "correlation")
dev.off()