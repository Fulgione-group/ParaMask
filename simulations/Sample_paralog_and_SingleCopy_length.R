# exponential distribution
# Specify x-values

# Plot dexp values 


seql<-1000000
ratio<-1/10

seqpl<- ratio*seql


# Sample length for paralogs
suml<-0
SV_length_sample<- c()
while(suml < seqpl){
  SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=1/1000)))
  suml<- sum(SV_length_sample)
}
seqscl<-seql-suml
sumscl<-0

# sample length for 

SC_length_sample<- c()

seqscl_rest<-seqscl
i<-1
for(i in 1:length(SV_length_sample)){
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+2)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[i]
}


SC_length_sample[length(SC_length_sample)+1]<- seqscl_rest

##convert to positions
SVtab <- c()
pos<- SC_length_sample[1]
SVtab<- c(1, pos, "sc")
i<-1
for(i in 1:length(SV_length_sample)){
  pos2<-pos + SV_length_sample[i]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
  pos <- pos2 
  pos2<- pos + SC_length_sample[i+1]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
  pos <- pos2 
}
SVtab <- data.frame(SVtab,stringsAsFactors = F)
colnames(SVtab)<- c("start", "end", "type")

write.table(x = SVtab, file = "~/Data/SeDuS/SVtab.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = SC_length_sample, file = "~/Data/SeDuS/SC_length_sample.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(x = SV_length_sample, file = "~/Data/SeDuS/SV_length_sample.txt", sep = "\t", row.names = F, quote = F, col.names = F)


##different ratios

seql<-1000000
ratio<-1/2

seqpl<- ratio*seql


# Sample length for paralogs
suml<-0
SV_length_sample<- c()
while(suml < seqpl){
  SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=1/1000)))
  suml<- sum(SV_length_sample)
}
seqscl<-seql-suml
sumscl<-0

# sample length for 

SC_length_sample<- c()

seqscl_rest<-seqscl
i<-1
for(i in 1:length(SV_length_sample)){
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+2)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[i]
}


SC_length_sample[length(SC_length_sample)+1]<- seqscl_rest

##convert to positions
SVtab <- c()
pos<- SC_length_sample[1]
SVtab<- c(1, pos, "sc")
i<-1
for(i in 1:length(SV_length_sample)){
  pos2<-pos + SV_length_sample[i]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
  pos <- pos2 
  pos2<- pos + SC_length_sample[i+1]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
  pos <- pos2 
}
SVtab <- data.frame(SVtab,stringsAsFactors = F)
colnames(SVtab)<- c("start", "end", "type")

write.table(x = SVtab, file = "~/Data/SeDuS/SVtab_0.5.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = SC_length_sample, file = "~/Data/SeDuS/SC_length_sample_0.5.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(x = SV_length_sample, file = "~/Data/SeDuS/SV_length_sample_0.5.txt", sep = "\t", row.names = F, quote = F, col.names = F)


##
seql<-1000000
ratio<-1/20

seqpl<- ratio*seql


# Sample length for paralogs
suml<-0
SV_length_sample<- c()
while(suml < seqpl){
  SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=1/1000)))
  suml<- sum(SV_length_sample)
}
seqscl<-seql-suml
sumscl<-0

# sample length for 

SC_length_sample<- c()

seqscl_rest<-seqscl
i<-1
for(i in 1:length(SV_length_sample)){
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+2)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[i]
}


SC_length_sample[length(SC_length_sample)+1]<- seqscl_rest

##convert to positions
SVtab <- c()
pos<- SC_length_sample[1]
SVtab<- c(1, pos, "sc")
i<-1
for(i in 1:length(SV_length_sample)){
  pos2<-pos + SV_length_sample[i]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
  pos <- pos2 
  pos2<- pos + SC_length_sample[i+1]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
  pos <- pos2 
}
SVtab <- data.frame(SVtab,stringsAsFactors = F)
colnames(SVtab)<- c("start", "end", "type")

write.table(x = SVtab, file = "~/Data/SeDuS/SVtab_0.05.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = SC_length_sample, file = "~/Data/SeDuS/SC_length_sample_0.05.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(x = SV_length_sample, file = "~/Data/SeDuS/SV_length_sample_0.05.txt", sep = "\t", row.names = F, quote = F, col.names = F)


##
seql<-1000000
ratio<-1/4

seqpl<- ratio*seql


# Sample length for paralogs
suml<-0
SV_length_sample<- c()
while(suml < seqpl){
  SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=1/1000)))
  suml<- sum(SV_length_sample)
}
seqscl<-seql-suml
sumscl<-0

# sample length for 

SC_length_sample<- c()

seqscl_rest<-seqscl
i<-1
for(i in 1:length(SV_length_sample)){
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+2)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[i]
}


SC_length_sample[length(SC_length_sample)+1]<- seqscl_rest

##convert to positions
SVtab <- c()
pos<- SC_length_sample[1]
SVtab<- c(1, pos, "sc")
i<-1
for(i in 1:length(SV_length_sample)){
  pos2<-pos + SV_length_sample[i]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
  pos <- pos2 
  pos2<- pos + SC_length_sample[i+1]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
  pos <- pos2 
}
SVtab <- data.frame(SVtab,stringsAsFactors = F)
colnames(SVtab)<- c("start", "end", "type")

write.table(x = SVtab, file = "~/Data/SeDuS/SVtab_0.25.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = SC_length_sample, file = "~/Data/SeDuS/SC_length_sample_0.25.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(x = SV_length_sample, file = "~/Data/SeDuS/SV_length_sample_0.25.txt", sep = "\t", row.names = F, quote = F, col.names = F)



##
seql<-1000000
ratio<-0.75

seqpl<- ratio*seql


# Sample length for paralogs
suml<-0
SV_length_sample<- c()
while(suml < seqpl){
  SV_length_sample<-c(SV_length_sample,round(rexp(1, rate=1/1000)))
  suml<- sum(SV_length_sample)
}
seqscl<-seql-suml
sumscl<-0

# sample length for 

SC_length_sample<- c()

seqscl_rest<-seqscl
i<-1
for(i in 1:length(SV_length_sample)){
  SC_length_sample<-c(SC_length_sample,rgeom(1, prob = (length(SV_length_sample)-i+2)/seqscl_rest))
  sumscl<- sum(SC_length_sample)
  seqscl_rest<- seqscl_rest-SC_length_sample[i]
}


SC_length_sample[length(SC_length_sample)+1]<- seqscl_rest

##convert to positions
SVtab <- c()
pos<- SC_length_sample[1]
SVtab<- c(1, pos, "sc")
i<-1
for(i in 1:length(SV_length_sample)){
  pos2<-pos + SV_length_sample[i]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "par"))
  pos <- pos2 
  pos2<- pos + SC_length_sample[i+1]
  SVtab <- rbind(SVtab, c((pos+1),pos2, "sc"))
  pos <- pos2 
}
SVtab <- data.frame(SVtab,stringsAsFactors = F)
colnames(SVtab)<- c("start", "end", "type")

write.table(x = SVtab, file = "~/Data/SeDuS/SVtab_0.75.txt", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(x = SC_length_sample, file = "~/Data/SeDuS/SC_length_sample_0.75.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(x = SV_length_sample, file = "~/Data/SeDuS/SV_length_sample_0.75.txt", sep = "\t", row.names = F, quote = F, col.names = F)




