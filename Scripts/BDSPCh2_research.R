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

ch2<- readLines("D2S441_100kb.vcf")
ch2_data<- read.table("D2S441_100kb.vcf", stringsAsFactors = FALSE)


#Subsetting 2 chromosome dataset into snip and STR
ch2_snip<- subset(ch2_data, V4== "A" | V4== "C"|V4== "T"| V4 =="G")
location<-  ch2_snip[c(2)]

#subtracting snip data from original data
library(dplyr)
ch2_str<- ch2_data[ !(ch2_data$V4 %in%  ch2_snip$V4), ]

#only taking snp name and genotype
df2<- ch2_snip[-c(1:2, 4:9)]

#transforming df2 into matrix and use substring to take only first three characters 
#for the genotype
m2<- as.matrix(df2)
for (i in 2:ncol(m2)){
  m2[,i]<- substr(m2[,i], start=1, stop=3)
}
#assigning new matrix so that old data can be checked
m22<- m2

#applying the function to the new matrix to count the ones
for (i in 1:nrow(m22)){
  for (j in 2:ncol(m22)){
    m22[i,j]<-change_discrete(m22[i,j])
  }
}
#putting  transformed matrix into new dataframe
s2<-as.data.frame(m22)
str(s2)
#using loop to make elements in genotype columns into numeric
for (i in 2:ncol(s2)){
  s2[,i]<-as.numeric(as.character(s2[,i]))
}
str(s2)
#creating an empty vector to store the frequency proportion of 1
freq_prop<-c()
#calculating the row sums and divide by 2 tims the number of individuals
freq_prop<- rowSums(s2[,2:ncol(s2)])/((ncol(s2)-1)*2)
summary(freq_prop)
#putting the freq_prop vector into the exisiting s2 dataframe
s2$freq_prop<-freq_prop
#visuals
library(ggplot2)
ggplot(aes(x=freq_prop), data= s2)+geom_histogram(color= "black", fill= "purple")+
  scale_x_continuous()+labs(title = "Chromosome 2", x= "frequency",y= "count")+
  theme(plot.title = element_text(hjust = 0.5))

#subsetting D2S441 from STR data
mstr2<- subset(ch2_str, V3 == "D2S441")
#Getting only STR name and genotype
mstr2<- mstr2[-c(1:2, 4:9)]
mmstr2<- as.matrix(mstr2)

#getting only numeric part of the STR genotype
for (i in 2:ncol(mmstr2)){
 mmstr2[,i]<-  gsub(":.*","",mmstr2[,i])
}
str<- mmstr2
dim(mmstr2)

#creating the loop to calculate average
for ( i in 2:ncol(mmstr2)){
  b<-num_split(mmstr2[1,i])
  c<-matrix(unlist(b))
  df<- as.data.frame(c)
  df$V1<- as.numeric(as.character(df$V1))
  avg<-( df[1,1]+df[2,1])/2
  mmstr2[1,i]<- avg
}

#changing the elements into numeric
dstr2<- as.data.frame(mmstr2)
for (i in 2:ncol(dstr2)){
  dstr2[,i]<-as.numeric(as.character(dstr2[,i]))
}

#removing the freq proportion from the latest SNP file
new_snp<- s2[-c(2506)]

#creating the new file by combining latest SNP and STR files
combo<- rbind( dstr2, new_snp)

#removing the names columns to calculate the correlation
new_combo<- combo[-c(1)]
cm<- as.matrix(new_combo)

#setting an empty list to put the row wise correlation
cors<- rep(NA, 958)
#using for loop to calculate the row wise correlation
for (i in 2:959){
  cors[i-1]<- cor(cm[1,], cm[i,])
}
#creating the new dataframe with first column as SNP and STR names
latest <- ch2_snip[c(2)]
which(is.na(cors))

#putting the correlation vector into new dataframe
latest$cors<- cors

#ignore the NA's
which(is.na(cors))

#picture
pdf(file="corplot.pdf", width=10, height=3)
plot(latest$V2, latest$cors, xlab= "location", ylab= "correlation")
dev.off()