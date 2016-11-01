#Open minus closed

#libraries
library(stringr)

source("/data/SHAKER_POPC_TYLER_Closed3/analysis_waters_GofR_TerminalAtom_new/auto_read_wats_TerminalAtom_closed.r")
source("/data/SHAKER_POPC_TYLER_Open/analysis_waters_GofR_TerminalAtom_new/auto_read_wats_TerminalAtom_open.r")

#keep only means
rm(list=setdiff(ls(), grep("*mean*", ls(), value=T)))




#make sure that there is the same number of means in both lists
res_mean_O <- grep("*mean*_O", ls(), value=T)
res_mean_C <- grep("*mean*_C", ls(), value=T)

eq_size <- function(A,B){
  if (!isTRUE(length(A)==length(B))){
    stop("the two means do not have the same number of residues" )
  } else {
    message("Equal number of residues in the vectors!")
  }
}

eq_size(res_mean_O,res_mean_C)
SimNum <- as.numeric(str_extract(res_mean_O, "[0-9]+"))
AA <- (str_extract(res_mean_O, "[aA-zZ]+"))

ShakerNum <- vector()
for (var in SimNum) {
	if (var <= 133){
		ShakerNum <- c(ShakerNum, var+204)
	}
	if (var > 133) {
		ShakerNum <- c(ShakerNum, var+218)
	}
}



#put everything in data frame
Residue <- gsub("_wat_mean_O","",res_mean_O)
shaker_data_OC <- data.frame(Residue,AA,SimNum,ShakerNum)
x <-vector()
y <- vector()


for (obj in res_mean_O) x <- c(x, get(obj))

shaker_data_OC$Open_mean <- x

for (obj in res_mean_C) y <- c(y, get(obj))

shaker_data_OC$Closed_mean <- y


#Get difference and percent increase
shaker_data_OC$ABS_diff <- shaker_data_OC$Open_mean-shaker_data_OC$Closed_mean

shaker_data_OC$Percent_Incr <- (shaker_data_OC$Open_mean-shaker_data_OC$Closed_mean)/shaker_data_OC$Closed_mean

shaker_data_OC <- shaker_data_OC[order(ShakerNum),]

shaker_data_OC_S1toS4 <- subset(shaker_data_OC, ShakerNum >= 224 & ShakerNum <= 377)

shaker_data_OC_S1 <- subset(shaker_data_OC, ShakerNum >= 224 & ShakerNum <= 248)

shaker_data_OC_S2 <- subset(shaker_data_OC, ShakerNum >= 278 & ShakerNum <= 300)

shaker_data_OC_S3 <- subset(shaker_data_OC, ShakerNum >= 311 & ShakerNum <= 329)

shaker_data_OC_S4 <- subset(shaker_data_OC, ShakerNum >= 336 & ShakerNum <= 377)