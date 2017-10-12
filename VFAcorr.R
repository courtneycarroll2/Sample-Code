## 
#Tests for differential abundance between C1 samples from alfalfa hay vs grass hay-fed alpacas; lines 25-171 of this script were sent to me by my professor and I adjusted it for my needs.
#Taxa with a significant diet effect are tested for correlation with volatile fatty acid abundances and plots are produced- I wrote this myself.
##


notaxa <- otus[,1:18]
sumofrows <- rowSums(notaxa)
sumofrows <- sort(sumofrows, decreasing = T)
hist(sumofrows,breaks = 200)
sumofrows  	## Use 200?

library(datasets)
library(dplyr)
library(lme4)
library(multcomp)

###
##
# FUNCTIONS
##
###


cor_prep <- function (xxx, yyy) {
  
  ## prepares the data for running correlations (i.e. puts the OTU table into formats at all taxonomic levels): INPUT: OTU Table; OUTPUT: OTU table, with data clustered at specified taxonomic levels
  ## xxx OTU table
  ## yyy taxonomic level
  
  ## aggregate it all and move data around; then transpose and add 'animal' and 'body' columns to split the name into two colums
  a1 <- aggregate(cbind(A1.C1,A2.C1,A3.C1,A4.C1,A5.C1,G1.C1,G2.C1,G4.C1,G1.C1L,G2.C1L,G3.C1L,G4.C1L,G5.C1L,A1.C1L,A2.C1L,A3.C1L) ~ get(yyy),xxx,sum)
  
#G3.C1,G2.C1,G3.C1L,G4.Ce,G1.C1,A5.C1,A1.Ce,A1.C1,G1.Ce,A2.Ce,A5.Ce,A3.Ce,A4.Ce,G4.C1,A3.C1L,A2.C1,A4.C1,G2.C1L,A2.C1.1,A3.C1,G3.Ce,A1.C1L,G5.Ce,G3.C1.1,G1.C1L,G2.Ce,G5.C1L,G4.C1L
  rownames(a1) <- a1[,1]
  a1 <- a1[,2:dim(a1)[2]]
  a1 <- t(a1)
  
  a1 <- data.frame(a1)
  a1$animal <- "xx"
  a1$body <- "yy"
  
  ## split the names into two columns (animal and body)
  for (i in 1:dim(a1)[1]) {
    xyz <- unlist(strsplit(rownames(a1)[i],split = "\\."))
    unlist(xyz)
    a1$animal[i] <- xyz[1]
    a1$body[i] <- xyz[2]
  }
  
  a1
}
compare_AG_abundance_all <- function (xxx,yyy,zzz) {
  
  ## xxx 	taxonomic level
  ## yyy 	otus file
  ## zzz 	column in the deg file
#   ## i		body site
#   
#  xxx <- genus
#  yyy <- AG
#  zzz <- "C1NDF"
#   i <- "C1L"
#   #  'KPCOFGS',AG,"AP-PP-BP",'C1'
    
  ## get the right data out of the OTU table (F= AG; T = AGQH)
  king <- droplevels(subset(xxx, is.na(xxx$body)==F))
  
  ## get the column out of the phenotype file
  Cec <- cbind(as.character(yyy$ID),as.character(yyy$Diet),as.numeric(as.character(yyy[,which(colnames(yyy)==zzz)])))
  Cec
  ## merge the two for the correlation
  b1 <- merge(Cec, king, by.x="V1",by.y="animal",all=F)
 
  ## make an output file
  output <- matrix(rep(0,(dim(b1)[2]-3)*6),ncol=6)
  dim(output)
  
  ## get rid of all other body sites
  c1 <- b1
  
  ## set the column name
  c1$V3 <- as.numeric(as.character(c1$V3))
  ## for every taxonomic group in file C1 (columns 3 to two before the end)...
  for (j in 4:(dim(c1)[2]-1)) {
    if (sum(c1[,j]) > 100) {
#      if (sum(c1[,j]>0)>5) {
        ## do the correlation		
        valA <- droplevels(subset(c1, c1$V2=="A"))
        valB <- droplevels(subset(c1, c1$V2=="G"))
        vala <- valA[,j]
        valb <- valB[,j]
        val <- try(t.test(vala, valb)$p.value)
        val3 <- try((t.test(vala, valb)$p.value)*dim(output)[1])
        val2 <- try(t.test(vala, valb)$statistic)
        ## get the column number
        #o <- which(names(table(list(b1$body)))==i)
        ## subtract 2 from the end
        p <- (j-2)
        ## the 'o-...' arranges the results within the right-sized matrix
        output[p,1] = "all"
        output[p,2] = colnames(b1)[j]		
        output[p,3] = zzz
        output[p,4] = as.numeric(as.character(val))
        output[p,5] = as.numeric(as.character(val3))
        output[p,6] = as.numeric(as.character(val2))
    }
  }
  ## outputs the output
  #	foo[order(foo$V1),]
  output3 <- output
  output2 <- output[order(output[,4]),]
  output2 <- data.frame(output2)
  output2 <- droplevels(subset(output2, output2$X1!=0))
  output2$X5 <- as.numeric(as.character(output2$X4))*dim(output2)[1]
  output2$fdr <- as.numeric(as.character(p.adjust(as.numeric(as.character(output2$X4)),method = "fdr")))
  output2
}

#----

###
# for using the body sites file 
###
## started with the file Kyle sent me, including correcting ht enames on the vials

otus <- read.csv('GA6350otus_filtered.csv', header=T, row.names=1, na.strings=".")
otus[1:10,1:10]

## make OTU tables, clustered at each taxonomic level
for (i in c("Phylum","Class","Order","Family", "Genus", "Species")) {
	assign(i,value = cor_prep(otus,i))
}

Phylum$V1 = NULL
Class$V1 = NULL
Order$V1 = NULL
Family$V1 = NULL
Species$V1 = NULL


otus$otu <- rownames(otus)
i="otu"
assign(i,value=cor_prep(otus,i))

## read in phenotype data
AG <- read.csv('AG_plusPCA.csv')

###
##
#	TARGETED BODY SITES - which microbes are signficantly correlated with VFA abundance in the C1 
##
###


for (i in c("Phylum","Class","Order","Family", "Genus", "Species","otu")) {
  for (j in c("C1NDF")) {
    
    ## make the file
    assign(paste("AGdif",i,j,sep="_"), compare_AG_abundance_all(get(i),AG,j))  
    
    ## make the file with only corrected values - remove the starting values, get the file made above, remove everything that is greater than a corrected value of 0.05, and then write out what remains with "sig" in front of it
    rm(aaa,bbb)
    aaa <- get(paste("AGdif",i,j,sep="_"))
    bbb <- droplevels(subset(aaa, aaa$fdr<0.05))
    if (dim(bbb)[1] > 0) {
      assign(paste("sigAGdif",i,j,sep="_"), bbb)  
    	#write.csv(get(paste("sigAGdif",i,j,sep="_")), paste("sigAGdif",i,j,".csv",sep="_"))
    }
  }
  }

########################### 
# sig_family_C1BUT
#lev <- "Family"; vvv <- t(Family)
#zzz <- "C1AC"
#taxon <-"Bacteroidaceae"
#num <- 1
#sigAGdif <- sigAGdif_Family_C1NDF
#xxx <- get(lev)
#vecC1AC <- vector()
#vecC1ACRho <- vector()
#vecC1BUT <- vector()
#vecC1BUTRho <- vector()
#vecC1PROP <- vector()
#vecC1PROPRho <- vector()
yyy <- AG
lev <- "Genus"
zzz <- "C1AC"
for (lev in c("Family", "Genus", "otu")){
	xxx <- get(lev)
	sigAGdif <- get(paste("sigAGdif", lev, "C1NDF", sep="_"))
	vecC1AC <- vector()
	vecC1ACRho <- vector()
	vecC1BUT <- vector()
	vecC1BUTRho <- vector()
	vecC1PROP <- vector()
	vecC1PROPRho <- vector()
	vecA <- vector()
	vecG <- vector()
	for (num in 1:dim(sigAGdif)[1] ){
		taxon <- sigAGdif[num,2]
		for (zzz in c("C1AC", "C1BUT", "C1PROP")){
			if (zzz=="C1AC"){
				vfa <- "Acetate"
			} else if (zzz =="C1BUT"){
				vfa <- "Butyrate"
			} else if (zzz =="C1PROP"){
				vfa <- "Propionate"
			}
			king <- droplevels(subset(xxx, (is.na(xxx$body)==F)))
			Cec <- cbind(as.character(yyy$ID),as.numeric(as.character(yyy[,which(colnames(yyy)==zzz)])))
			b1 <- merge(Cec, king, by.x="V1",by.y="animal",all=F)
			something<-cor.test(x=as.numeric(as.character(b1$V2)), y=as.numeric(as.character(b1[[paste(taxon)]])), method="spearman", exact=F)
			something
			pval <- something$p.value
			rho <- something$estimate
			AGmeans <- t.test(b1[[paste(taxon)]][1:8],b1[[paste(taxon)]][9:16])
			#t = 2.4619, df = 7.6409, p-value = 0.04056
			# Alfalfa 3.000    	Grass  122.125 
			forplot <- b1 
			forplot$diet<-c(rep("alfalfa",8),rep("grass",8))
	
			xlabel <- paste(vfa, "Extraction (mM)", sep=" ")
			ylabel <- paste(taxon, "Abundance", sep = " ")
			jpeg(file = paste(taxon, "_", vfa, ".jpeg", sep = "")) 
			plot(x=as.numeric(as.character(forplot$V2)), y=as.numeric(as.character(forplot[[paste(taxon)]])), xlab=xlabel , ylab= ylabel, cex.lab = 1.1, cex=2.2, cex.axis=1.0, col=ifelse(forplot$diet=="alfalfa", "red", "blue"))
		dev.off()
			AMean <- AGmeans$estimate[1]
			GMean <- AGmeans$estimate[2]
			if (zzz=="C1AC"){
				vecC1AC[num] <- as.numeric(as.character(pval))
				vecC1ACRho[num] <- as.numeric(as.character(rho))
				vecG[num] <- as.numeric(as.character(GMean))
				vecA[num] <- as.numeric(as.character(AMean))
			} else if (zzz =="C1BUT"){
				vecC1BUT[num] <- as.numeric(as.character(pval))
				vecC1BUTRho[num] <- as.numeric(as.character(rho))
			} else if (zzz =="C1PROP"){
				vecC1PROP[num] <- as.numeric(as.character(pval))
				vecC1PROPRho[num] <- as.numeric(as.character(rho))
			}
			rm(king, Cec, b1, something, pval, rho, AGmeans, forplot, xlabel, ylabel)
			## from here, VFA changes
		}
		rm(vfa)
		## from here, taxon changes
	}
	sigAGdif$Acetate <- vecC1AC
	sigAGdif$AcRho <- vecC1ACRho
	sigAGdif$Butyrate <- vecC1BUT
	sigAGdif$ButRho <- vecC1BUTRho
	sigAGdif$Propionate <- vecC1PROP
	sigAGdif$PropRho <- vecC1PROPRho
	sigAGdif$Amean <- vecA
	sigAGdif$Gmean <- vecG
	sigAGdif$fdrAC <- as.numeric(as.character(p.adjust(as.numeric(as.character(sigAGdif$Acetate)),method = "fdr")))
	sigAGdif$fdrBUT <- as.numeric(as.character(p.adjust(as.numeric(as.character(sigAGdif$Butyrate)),method = "fdr"
	)))
	sigAGdif$fdrPROP <- as.numeric(as.character(p.adjust(as.numeric(as.character(sigAGdif$Propionate)),method = "fdr")))
	filename <- paste(lev, "Stats", sep="")
	write.csv(sigAGdif, filename)
	rm(sigAGdif, filename, taxon, num, vecC1AC, vecC1ACRho, vecC1BUT, vecC1BUTRho, vecC1PROP, vecC1PROPRho, lev, vecA, vecG)
	## from here, level changes
}

	