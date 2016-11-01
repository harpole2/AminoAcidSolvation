library(plyr)

library(ggplot2)

#working directory
setwd("/data/SHAKER_POPC_TYLER_Open/analysis_waters_GofR_TerminalAtom/")

test <- "LeuProAspAsnGluLysGlnArgLysValTrpLeuLeuPheGluTyrProGluSerSer
GlnAlaAlaArgValValAlaIleIleSerValPheValIleLeuLeuSerIleValIle
PheCysLeuGluThrLeuProGluPheLysHisTyrLysValPheAsnThrThrThrAsn
GlyThrLysIleGluGluAspGluValProAspIleThrAspProPhePheLeuIleGlu
ThrLeuCysIleIleTrpPheThrPheGluLeuThrValArgPheLeuAlaCysProAsn
LysLeuAsnPheCysArgAspValMetAsnValIleAspIleIleAlaIleIleProTyr
PheIleThrLeuAlaThrValValAlaGluGluGluAspSerSerAsnGlnAlaMetSer
LeuAlaIleLeuArgValIleArgLeuValArgValPheArgIlePheLysLeuSerArg
HisSerLysGlyLeuGlnIleLeuGlyArgThrLeuLysAlaSerMetArgGluLeuGly
LeuLeuIlePhePheLeuPheIleGlyValValLeuPheSerSerAlaValTyrPheAla
GluAlaGlySerGluAsnSerPhePheLysSerIleProAspAlaPheTrpTrpAlaVal
ValThrMetThrThrValGlyTyrGlyAspMetThrProValGlyValTrpGlyLysIle
ValGlySerLeuCysAlaIleAlaGlyValLeuThrIleAlaLeuProValProValIle
ValSerAsnPheAsnTyrPheTyrHisArgGluThr"

test <- gsub('([[:upper:]])', ' \\1', test)

test <- toupper(test)

test <- gsub("[\r\n]", "", test)

#trim leading trailing space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
test <- trim(test)
test <- strsplit(test, split=" ")

seq <- vector()
for (does in test[[1]]) seq <- c(seq, does)

i <- 0
j <- 272
k <- 544
l <- 816

for (AA in seq) {
	tmp1 <-paste("Water-within-GofR-of-residue",i,"-",AA,"-notbackbone.dat", sep="")
	tmp2 <-paste("Water-within-GofR-of-residue",j,"-",AA,"-notbackbone.dat", sep="")
	tmp3 <-paste("Water-within-GofR-of-residue",k,"-",AA,"-notbackbone.dat", sep="")
	tmp4 <-paste("Water-within-GofR-of-residue",l,"-",AA,"-notbackbone.dat", sep="")

	tmp_tab1 <- read.table(tmp1, header=T)
	tmp_tab2 <- read.table(tmp2, header=T)
	tmp_tab3 <- read.table(tmp3, header=T)
	tmp_tab4 <- read.table(tmp4, header=T)

	tab1 <- paste(AA,i,"_A_O", sep="")
	tab2 <- paste(AA,i,"_B_O", sep="")
	tab3 <- paste(AA,i,"_C_O", sep="")
	tab4 <- paste(AA,i,"_D_O", sep="")

	assign(tab1, tmp_tab1)
	assign(tab2, tmp_tab2)
	assign(tab3, tmp_tab3)
	assign(tab4, tmp_tab4)

	data_all <- data.frame(tmp_tab1["frame"], tmp_tab1["Corrected_water"], 
			tmp_tab2["Corrected_water"], tmp_tab3["Corrected_water"], tmp_tab4["Corrected_water"])

	data_all <- rename(data_all, c("Corrected_water"="water_1", 
                                 "Corrected_water.1"="water_2", 
                                 "Corrected_water.2"="water_3", 
                                 "Corrected_water.3"="water_4"))

	data_all$mean <- rowMeans(data_all[,-1])
	
	data_mean <- paste(AA,i,"_wat_mean_O", sep="")
	data_sd <- paste(AA,i,"_wat_sd_O", sep="")

	assign(data_mean, mean(data_all$mean))
	assign(data_sd, sd(data_all$mean))
	i <- i+1
	j <- j+1
	k <- k+1
	l <- l+1
}
