####
##
#
### Takes an OTU table and gives read count of each phylum/class/order/genus/species/OTU for each sample group ("AlfalfaC1", "GrassC1", etc.)




library(dplyr)
library(datasets)
library(multcomp)
library(lme4)


otutable <- read.table('alltaxa_sums5990.csv',sep=",", header=T,row.names = 1)

for (m in c("Phylum", "Class", "Order", "Family", "Genus", "Species")){
## get function lets you call something in the function 	
#m <- "OTUID"
## group all columns by phylum
aaa <- aggregate(cbind(AlfalfaC1,AlfalfaC1.1,AlfalfaC1.2,AlfalfaC1.3,AlfalfaC1.4,AlfalfaCecum,AlfalfaCecum.1,AlfalfaCecum.2,AlfalfaCecum.3,AlfalfaCecum.4,AlfalfaIleum,AlfalfaIleum.1,AlfalfaIleum.2,AlfalfaIleum.3,AlfalfaIleum.4,AlfalfaJej,AlfalfaJej.1,AlfalfaJej.2,AlfalfaJej.3,AlfalfaLI,AlfalfaLI.1,AlfalfaLI.2,AlfalfaLI.3,AlfalfaLI.4,GrassC1,GrassC1.1,GrassC1.2,GrassC1.3,GrassC1.4,GrassCecum,GrassCecum.1,GrassCecum.2,GrassCecum.3,GrassCecum.4,GrassIleum,GrassIleum.1,GrassIleum.2,GrassIleum.3,GrassIleum.4,GrassJej,GrassJej.1,GrassJej.2,GrassLI,GrassLI.1,GrassLI.2,GrassLI.3,GrassLI.4) ~ get(m), otutable,sum)

ccc <- matrix(nrow=dim(aaa)[1], ncol=11)
ccc <- data.frame(ccc, stringsAsFactors = F)
ccc$X1 <- aaa$`get(m)`

ccc$X2 <- by(aaa[,c(2:6)], aaa$`get(m)`, FUN=rowSums)
ccc$X3 <- by(aaa[,c(17:20)], aaa$`get(m)`, FUN=rowSums)
ccc$X4 <- by(aaa[,c(12:16)], aaa$`get(m)`, FUN=rowSums)
ccc$X5 <- by(aaa[,c(7:11)], aaa$`get(m)`, FUN=rowSums)
ccc$X6 <- by(aaa[,c(21:25)], aaa$`get(m)`, FUN=rowSums)
ccc$X7 <- by(aaa[,c(26:30)], aaa$`get(m)`, FUN=rowSums)
ccc$X8 <- by(aaa[,c(41:43)], aaa$`get(m)`, FUN=rowSums)
ccc$X9 <- by(aaa[,c(36:40)], aaa$`get(m)`, FUN=rowSums)
ccc$X10 <- by(aaa[,c(31:35)], aaa$`get(m)`, FUN=rowSums)
ccc$X11 <- by(aaa[,c(44:48)], aaa$`get(m)`, FUN=rowSums)


write.csv(ccc, paste(m, "sums.csv", sep=""))
rm(ccc,aaa)
}


