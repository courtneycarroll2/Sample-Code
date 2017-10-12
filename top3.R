####
# Reads in an OTU table and outputs the 3 most abundant taxa (phylum through OTU) for each sample group.
####


library(dplyr)
library(datasets)
library(multcomp)
library(lme4)


otutable <- read.table('table_even5990.csv',sep=",", header=T,row.names = 1)

ccc <- matrix(c("AC1","# Reads", "# Otus", "AJej", "# Reads", "# Otus", "AIle", "# Reads", "# Otus", "ACec", "# Reads", "# Otus", "ALI", "# Reads", "# Otus", "GC1", "# Reads", "# Otus", "GJej", "# Reads", "# Otus", "GIle", "# Reads", "# Otus", "GCec", "# Reads", "# Otus", "GLI", "# Reads", "# Otus"),ncol=30)
ccc <- data.frame(ccc, stringsAsFactors = F)

for (m in c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")){
## get function lets you call something in the function 	
#m <- "OTUID"
## group all columns by phylum
aaa <- aggregate(cbind(AlfalfaC1,AlfalfaC1.1,AlfalfaC1.2,AlfalfaC1.3,AlfalfaC1.4,AlfalfaCecum,AlfalfaCecum.1,AlfalfaCecum.2,AlfalfaCecum.3,AlfalfaCecum.4,AlfalfaIleum,AlfalfaIleum.1,AlfalfaIleum.2,AlfalfaIleum.3,AlfalfaIleum.4,AlfalfaJej,AlfalfaJej.1,AlfalfaJej.2,AlfalfaJej.3,AlfalfaLI,AlfalfaLI.1,AlfalfaLI.2,AlfalfaLI.3,AlfalfaLI.4,GrassC1,GrassC1.1,GrassC1.2,GrassC1.3,GrassC1.4,GrassCecum,GrassCecum.1,GrassCecum.2,GrassCecum.3,GrassCecum.4,GrassIleum,GrassIleum.1,GrassIleum.2,GrassIleum.3,GrassIleum.4,GrassJej,GrassJej.1,GrassJej.2,GrassLI,GrassLI.1,GrassLI.2,GrassLI.3,GrassLI.4) ~ get(m), otutable,sum)
#colnames(aaa[2:6]) ## AC1
#colnames(aaa[17:20]) ## AJej
#colnames(aaa[12:16]) ## AIle
#colnames(aaa[7:11]) ## ACec
#colnames(aaa[21:25]) ## ALI
#colnames(aaa[26:30]) ## GC1
#colnames(aaa[41:43]) ## GJej
#colnames(aaa[36:40]) ## GIle
#colnames(aaa[31:35]) ## GCec
#colnames(aaa[44:48]) ## GLI


aaa$AC1sums <- rowSums(aaa[2:6])
aaa$AJejsums <- rowSums(aaa[17:20])
aaa$AIlesums <- rowSums(aaa[12:16])
aaa$ACecsums <- rowSums(aaa[7:11])
aaa$ALIsums <- rowSums(aaa[21:25])
aaa$GC1sums <- rowSums(aaa[26:30])
aaa$GJejsums <- rowSums(aaa[41:43])
aaa$GIlesums <- rowSums(aaa[36:40])
aaa$GCecsums <- rowSums(aaa[31:35])
aaa$GLIsums <- rowSums(aaa[44:48])

aaa <- droplevels(subset(aaa, !aaa$`get(m)`==""))

sortedAC1 <- aaa[order(-aaa$AC1sums),]
sortedAJej <- aaa[order(-aaa$AJejsums),]
sortedAIle <- aaa[order(-aaa$AIlesums),]
sortedACec <- aaa[order(-aaa$ACecsums),]
sortedALI <- aaa[order(-aaa$ALIsums),]
sortedGC1 <- aaa[order(-aaa$GC1sums),]
sortedGJej <- aaa[order(-aaa$GJejsums),]
sortedGIle <- aaa[order(-aaa$GIlesums),]
sortedGCec <- aaa[order(-aaa$GCecsums),]
sortedGLI <- aaa[order(-aaa$GLIsums),]


subsetotus <- data.frame(otutable) 
#p <- "GLI"


for (i in 1:3){
	vec <- vector()
	for (p in c("AC1", "AJej", "AIle", "ACec", "ALI", "GC1", "GJej", "GIle", "GCec", "GLI")){
		if (p == "AC1"){
			j <- subsetotus[,c(1:5, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 1
		}else if(p == "AJej"){
			j <- subsetotus[,c(16:19, 49:55)]
			j$sum <- rowSums(j[1:4])
			num <- 2
		}else if(p == "AIle"){
			j <- subsetotus[,c(11:15, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 3
		}else if(p == "ACec"){
			j <- subsetotus[,c(6:10, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 4
		}else if(p == "ALI"){
			j <- subsetotus[,c(20:24, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 5
		}else if (p == "GC1"){
			j <- subsetotus[,c(25:29, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 6
		}else if(p == "GJej"){
			j <- subsetotus[,c(40:42, 49:55)]
			j$sum <- rowSums(j[1:3])
			num <- 7
		}else if(p == "GIle"){
			j <- subsetotus[,c(35:39, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 8
		}else if(p == "GCec"){
			j <- subsetotus[,c(30:34, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 9
		}else if(p == "GLI"){
			j <- subsetotus[,c(43:47, 49:55)]
			j$sum <- rowSums(j[1:5])
			num <- 10
		}
	
		j <- droplevels(subset(j, !j$sum==0))

		x <- get(paste("sorted", p, sep=""))
		a <- droplevels(subset(j, j[m]==as.character(x[i,1])))
		a1 <- dim(a)[1]
		vec[num] <- a1
		rm(a, a1, j, x, num)
	}

		ccc <- rbind(ccc, c(as.character(sortedAC1$`get(m)`[i]), as.character(sortedAC1$AC1sums[i]), as.character(vec[1]), as.character(sortedAJej$`get(m)`[i]), as.character(sortedAJej$AJejsums[i]), as.character(vec[2]), as.character(sortedAIle$`get(m)`[i]), as.character(sortedAIle$AIlesums[i]), as.character(vec[3]), as.character(sortedACec$`get(m)`[i]), as.character(sortedACec$ACecsums[i]), as.character(vec[4]), as.character(sortedALI$`get(m)`[i]), as.character(sortedALI$ALIsums[i]), as.character(vec[5]), as.character(sortedGC1$`get(m)`[i]), as.character(sortedGC1$GC1sums[i]), as.character(vec[6]), as.character(sortedGJej$`get(m)`[i]), as.character(sortedGJej$GJejsums[i]), as.character(vec[7]), as.character(sortedGIle$`get(m)`[i]), as.character(sortedGIle$GIlesums[i]), as.character(vec[8]), as.character(sortedGCec$`get(m)`[i]), as.character(sortedGCec$GCecsums[i]), as.character(vec[9]), as.character(sortedGLI$`get(m)`[i]), as.character(sortedGLI$GLIsums[i]), as.character(vec[10])))
}

ccc <- rbind(ccc, c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""))
}

write.csv(ccc, "5990Stats.csv")
