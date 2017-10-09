#As coded, this script is designed to implement the scanone() function of Rqtl with the EM algorithm. 
#See A Guide to Mapping with R/QTL by Karl Broman (ISBN: 978-0-387-92124-2) for a comprehensive tutorial and guidelines.
#Also a helpful forum moderated by Dr. Broman (as of the date of this publication) 
#can be found at http://groups.google.com/group/rqtl-disc/

#This script is built on v6 of Spindel et al., 2013. Theoretical and Applied Genetics. 2013, 126, 2699–2716
#A version of this script, alongside a sample file, is available from SolGenomics.net: https://github.com/solgenomics/snpbinner
#In this version (Gonda et al., 2017) we have added genotype exploration steps to plot and extract recombination frequencies data.
#The LOD scores from the model part are generated into files.
#This version is fully parallelized (excluding step 21) for a use with multiple CPUs server.
##It is highly recommended that you check your phenotype and genotype data before running the automated model selection.

####LOAD THE CROSS####

#(1) Load R/qtl. Make sure you have the following R packages installed

##########################################################################################
library(qtl)
library(rlecuyer)
library(snow)
library(foreach)  
library(doSNOW)  
cl<-makeCluster(2) #change the 2 to your number of CPU cores  
registerDoSNOW(cl)
##########################################################################################

#(2) Read in the genotype and phenotype data. Be sure the file structure is of class "bc"
#for a RIL or a backcross population. As in AA=A, and BB=H, AB=NA. Otherwise structure the 
#data as an intercross (AA=A, BB=B, and AB=H). Hets are treated as missing data in a selfed 
#RIL population. See ?read.cross() for more information. 

#READING THE DATA IN WITH THIS STATEMENT SIMULTANEOUSLY ESTIMATES THE GENETIC MAP POSITIONS 
#OF EACH OF THE MARKERS IF NOT ALREADY PRESENT IN YOUR DATA FILE. DEPENDING ON THE NUMBER OF 
#INDIVIDUALS AND THE NUMBER OF MARKERS, THIS FUNCTION COULD TAKE A WHILE. IF POSISTION INFORMATION 
#IN YOUR DATA FILE REPRESENTS PHYSICAL DISTANCES YOU MAY CONSIDER USING rescalemap() TO CONVERT 
#PHYSICAL POSITIONS TO CENTIMORGANS. BE SURE TO USE A SCALE APPROPRIATE FOR YOUR ORGANISM. 
#SCALE ARGUMENT SHOULD BE SET TO 1/(#bp	in a cM). SEE ?rescalemap() FOR MORE INFO.

#Note the "csvr" format. If your data is not formatted as a rotated comma separated file,
#then the format argument should be changed to reflect it's format. See ?read.cross() for more info.
#Sample files can be find at this link: https://github.com/solgenomics/snpbinner

##########################################################################################
cross<-read.cross(format ="csvr", file="your_cross_file.csv", estimate.map=F)
##########################################################################################

#(3) Change population class from bc to RIL if you are working with a RIL population. 
#If you are working with a bc or an intercross and coded your data appropriately, you may skip this step. 

#R/QTL has strange support for RILs, it assumes you have an F2 if you code the genotypes correctly. 
#Hence why your file structure has to simulate a BC. In a BC there are only two genotype classes (AA and AB). 
#Likewise you only have two genotype classes in a RIL population (AA and BB). 
#So if in the data file you marked all of your BB individuals with an H (ie. made your RILs look like a BC) 
#this next step will tell R/QTL that those *H* genotypes are actually BB and not hets, and account 
#for the map expansion that occurs in a RIL population. Use convert2risib(cross) for RILs developed by sib-mating. 

##########################################################################################
cross<-convert2riself(cross)
newmap=est.map(cross, map.function="kosambi", n.cluster=12, tol=0.01,
               maxit=1000)
cross= replace.map(cross, newmap)
##########################################################################################

#(4) Slightly adjust marker positions to to avoid identical positions for markers on different chromosomes.

##########################################################################################
cross<-jittermap(cross)
##########################################################################################

#(5) Calculate the underlying genotype probabilities using a kosambi recombination model.

##########################################################################################
cross<-calc.genoprob(cross, map.function="kosambi")
##########################################################################################


############################
####CHECK PHENOTYPE DATA####
############################


##########################################################################################

#(6) Generate Histogram plots to visualize each phenotypic distribution.

##########################################################################################
pdf(file="1_Phenotype Histograms_transformed.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)) {
	plot.pheno(cross, pheno.col=i)
	}
dev.off()	
##########################################################################################

#(7) Generate .csv file for Shapiro-wilk normality test. Low p-values indicate departures from normality.

##########################################################################################
norm<-{}
norm2<-{}
for(i in 1:nphe(cross)){
	x<-cross$pheno[i]
	x<-x[,1]
	y<-shapiro.test(x)
	ShapiroWilk.pvalue<-y$p.value
	Phenotype<-colnames(cross$pheno[i])
	norm<-cbind(Phenotype, ShapiroWilk.pvalue)
	norm2<-rbind(norm2, norm)
}
write.csv(file="2_Normality_test_p-values_post_transformation.csv", norm2)
##########################################################################################

#(8) Generate pdf of plot to check for batch effects.

##########################################################################################
pdf(file="3_Batch effects.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)){
	par(mfrow=c(1,2), las=1, cex=0.8)
	means<-apply(cross$pheno[i], 1, mean)
	plot(means)
	plot(sample(means), xlab="Random Index", ylab="means", main=colnames(cross$pheno[i]))
}
dev.off()
##########################################################################################


###########################
####CHECK GENOTYPE DATA####
###########################


#(9) Generate a graphical genotype.

##########################################################################################
pdf(file="4_Graphical Genotype.pdf", width=11, height=8.5)
geno.image(cross)
dev.off()
##########################################################################################

#(10) Generate a genetic map and save as a pdf. 
#Look for map expansion to identify markers whose map position may be in error.

##########################################################################################
pdf(file="5_Estimated Genetic Map.pdf", width=11, height=8.5)
par(cex=0.6)
plot.map(cross, show.marker.names=FALSE)
dev.off()
##########################################################################################

#(11.1) plot the recombination frequencies.
#Estimated recombination fractions (upper left) and LOD scores (lower right) for all pairs of markers.

##########################################################################################
png(file="6_recombination frequncies.png", width=2400, height=1900, pointsize = 48)
cross <- est.rf(cross)
plot.rf(cross, col.scheme=c("redblue"))
dev.off()
##########################################################################################

#(11.2) Plot the recombination frequencies for the single chromosomes
#If you prepared your cross file in excell, make sure chromosome column type is well defined (best is to pase format from sample file)

##########################################################################################
for(i in 1:nchr(cross)){
file<-paste("6_recombination frequncies_ch_", i, ".png")
png(filename=file, width=2400, height=1900, pointsize = 48)
plot.rf(cross, chr=c(i), col.scheme=c("redblue"))
dev.off()
}
##########################################################################################

#(11.3) create csv files for the recombination frequencies
#A very large file size

##########################################################################################
{
rf <- pull.rf (cross)
write.csv(file="6_rf_all_chr.csv", rf)
}
##########################################################################################

#(11.4) create csv files for the recombination frequencies for each chromosome

##########################################################################################
for(i in 1:nchr(cross)){
file<-paste("6_rf_ch_", i, ".csv")
rf <- pull.rf(subset(cross, chr=i))
write.csv(file=file, rf)
}
##########################################################################################

#(11.5) create csv files for the recombination frequencies
#A very large file size

##########################################################################################
{
lod <- pull.rf(cross, what="lod")
write.csv(file="6_rf_lod_all_chr.csv", lod)
}
##########################################################################################

#(11.6) create csv files for the recombination frequencies LOD scores

##########################################################################################
for(i in 1:nchr(cross)){
file<-paste("6_rf_lod_ch_", i, ".csv")
lod <- pull.rf(subset(cross, chr=i), what="lod")
write.csv(file=file, lod)
}
##########################################################################################

#(12) Generate a .csv file with p-values for a X^2 test for segregation distortion 
#i.e. departures from mendelian expectation.

##########################################################################################
sd<-geno.table(cross)
sd<-sd[ sd$P.value < 1e-5, ]
write.csv(file="7_Chi-square for segregation distortion.csv", sd)
##########################################################################################

#(13) Generate a pdf of a histogram to compare genotypes for each pair of individuals 
#and identify pairs that have unusually similiar genotypes. i.e. sample mix up.

##########################################################################################
pdf(file="8_Histogram of genotype comparisons.pdf")
genotype.comparisons<-comparegeno(cross)
hist(genotype.comparisons, breaks=200, xlab="Proportion of identical genotypes")
rug(genotype.comparisons)
dev.off()
##########################################################################################

#(14) Generate csv file to identify which individuals are unusually similar and their proportion of similiarity

##########################################################################################
cg.high<-which(genotype.comparisons>0.95, arr.ind=TRUE)
proportion<-{}
for(i in 1:nrow(cg.high)){
	x<-genotype.comparisons[cg.high[i,1], cg.high[i,2]]
	proportion<-c(proportion, x)
	}
Genotype1<-cg.high[,1]
Genotype2<-cg.high[,2]
cg.high<-cbind(Genotype1, Genotype2, proportion)
write.csv(file="9_Unusually similiar genotypes.csv", cg.high)
##########################################################################################

#(15) Identify markers where >15% of the population is NA and export as csv. 
#To change this threshold change the value of x<-x*0.15 to the desired tolerance.

##########################################################################################
x<-(nind(cross))
x<-x*0.15
missing.ind<-nmissing(cross, what="mar")
missing.ind<-missing.ind[missing.ind > x]
nmi<-"Number of Individuals without marker data"
missing.ind<-c(nmi, missing.ind)
write.csv(file="10_Poorly typed markers.csv", missing.ind)
##########################################################################################

#(16) Identify Individuals where >15% of the markers is NA and export as csv. 
#To change this threshold change the value of x<-x*0.15 to the desired tolerance.

##########################################################################################
x<-(totmar(cross))
x<-x*0.15
number.missing.mar<-nmissing(cross, what="ind")
Line.ID<-1:nind(cross)
missing.mar<-cbind(Line.ID, number.missing.mar)
missing.mar<-missing.mar[missing.mar[,2]>x,]
write.csv(file="11_Poorly typed individuals.csv", missing.mar)
##########################################################################################

############################################
####INITIAL GENOME SCAN FOR ADDITIVE QTL####
############################################

##########################################################################################
#YOU MAY CONSIDER SPLITTING NORMAL FROM NON-NORMAL AND/OR BINARY PHENOTYPES AND RUNNING SEPARATELY.

#FOR BINARY TRAITS CHANGE THE model="normal" ARGUMENTS TO model="binary" FOR STEPS #(17), #(18), AND #(26) ONLY. 

#FOR CONSIDERATION OF NON-PARAMETRIC TRAITS CHANGE model="normal" TO model="np" ONLY FOR STEPS #(17), AND #(18).

##########################################################################################

#(17) Perform single QTL Scan for all phenotypes using Haley-Knott Regression. 
#This may produce a warning if there are individuals with missing phenotype data. 
#If this occurs it does not change the output of the code. 

#*NOTE*# THIS STEP MAY GIVE A NON-CONVERGENCE ERROR FOR BINARY TRAITS.

##########################################################################################
cross.sc1<-scanone(cross, pheno.col=1:nphe(cross), method="em", use="all.obs", model="normal")
##########################################################################################

#(18) Calculate LOD significance thresholds based on 1000 permutations<--THIS COULD TAKE A WHILE. 
#If you have a multi-core processor, you may consider invoking the SNOW package ' library(snow) ' 
#and specifying the argument n.cluster with the number of parallel computations you'd like to do 
#(usually the number of chromosomes in the genome). 

#*NOTE*# THIS STEP MAY GIVE A NON-CONVERGENCE ERROR FOR BINARY TRAITS.

##########################################################################################
cross.sc1.perms<-scanone(cross, pheno.col=1:nphe(cross), method="em", n.perm=1000, verbose=TRUE, model="normal", n.cluster=12, tol=0.001)
##########################################################################################

#(19) Produce a text file indicating the most signficant marker on each chromosme is extracted (they may or may not cross the significance threshold determined by the permutation test.), 
#the drop 1.5 LOD intervals, LODs, and pvalues for QTL detected.

##########################################################################################
sum<-summary(cross.sc1, threshold=0, perms=cross.sc1.perms, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE)
space<-" "
foreach(i=1:nphe(cross)) %dopar% {
	x<-capture.output(sum[[i]])
	cat(colnames(cross$pheno[i]), file="12_Initial QTL hits by phenotype.txt", sep="\n", append=TRUE)
	cat(x, file="12_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
	cat(space, file="12_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
}
##########################################################################################

#(20) Export csv file with LOD scores at every marker. 

##########################################################################################
data <- data.frame(cross.sc1)

foreach(i=1:nphe(cross)) %dopar% {
	library(qtl)
	ilods <- data[,c(1,2,i+2)]
	phenotype<-i
	pheno<-colnames(cross$pheno[phenotype])
	file<-paste("13_LOD scores for every marker_", pheno,".txt")
 	write.table(ilods, file=file, col.names=NA, sep="\t")
}
##########################################################################################

#(21) Produce PDF file indicating genome scans with alpha=.05 threshold marked

##########################################################################################
z<-summary(cross.sc1.perms, alpha=.05)
pdf(file="14_QTL Plots.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)){
plot(cross.sc1, lodcolumn=i, lwd=1.5, gap=0, bandcol="gray70", incl.markers=TRUE, main=colnames(cross$pheno[i]), xlab=c("Threshold for alpha=.05 using 1000 permutations", z[i]))
add.threshold(cross.sc1, perms=cross.sc1.perms, alpha=0.05, lodcolumn=i, gap=0)
}
dev.off()
##########################################################################################

#(22) Produce a csv file with the map positions for each marker.

#*NOTE*# This may produce a warning consistent with the number of chromosomes. 
#This is expected as the program is making the user aware that it is appending the column names to the file. 

##########################################################################################
newmap<-pull.map(cross)
for(i in 1:length(names(newmap))){
	snps<-names(newmap[[i]])
	gm<-c(snps, newmap[[i]])
	gm2<-matrix(gm, ncol=2)
	write.table(file="15_Genetic Map Positions.csv", sep=",", append=TRUE, gm2)
}
##########################################################################################

#(23)run sim.geno (necessary for some of the following functions)

##########################################################################################
cross2<-sim.geno(cross)
##########################################################################################

######################################################################################
########																	##########
######## Steps #(24) - #(34) represent a loop that will scan for additional ##########
######## qtl, run model selection, and fit the qtl model for each phenotype ##########
######## in turn automatically. Run the entire code (Steps #(24) - #(34) )  ##########
######## simultaneously.                                                    ##########
########																	##########
######################################################################################


######################################
##Model Selection and QTL refinement##
######################################

#(24) Specifies the phenotype you wish to work with.
#If you are troubleshooting or otherwise want to run each step individually for a given phenotype rather than the whole things a loop,
#just replace k in ' phenotype<-k ' with the number of the phenotype you want to use, ignore the line ' for(k in 2:nphe(cross)){ '  and run each step individually. 
#If you run multiple phenotypes you do not need to change anything.
#After last loop (phenotype) finished, an Error on failed task will probably apear "no loop for break/next, jumping to top level"

##########################################################################################
foreach(k=1:nphe(cross)) %dopar% {

library(qtl)

phenotype<-k

pheno<-colnames(cross$pheno[phenotype])
##########################################################################################

#(25) Make a qtl object containing the ALL the qtl from file="12_Initial QTL hits by phenotype" 
#for the phenotype you specified.  

##########################################################################################
chromo<-sum[[pheno]]
chr<-{}
pos<-{}
for(i in 1:nrow(chromo)){
chr1<-chromo[i,1]
chr2<-as.numeric(as.character(chr1))
chr<-c(chr, chr2)
}
for(i in 1:nrow(chromo)){
pos1<-chromo[i,2]
pos<-c(pos, pos1)
}
if(is.na(chr[1])==FALSE) {
qtl<-makeqtl(cross, chr=chr, pos=pos, what="prob")
createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
##########################################################################################

#(26) Scan for additional linked QTL conditioning on the QTL already detected. .

#*NOTE*# You may receive warning messages about dropping individuals with missing phenotype data 
#and/or that the column names in scanone input do not match those in perms input. These are both 
#expected under some circumstances and do not effect the output of the code. Also you may recieve 
#a warning that there is no chromosome number NA if no additional QTL are to be found.

##########################################################################################
cross.aq<-addqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", model="normal")
sub.perms<-subset(cross.sc1.perms, lodcolumn=phenotype)
xx<-capture.output(summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE))
xxy<-c(pheno,xx)
cat(xxy, file="16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
cat(space, file="16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
##########################################################################################

#(27) Add additional QTL (if any) in the file "14_Additional QTL hits by phenotype.txt" 
#to the qtl object and update the formula object.

##########################################################################################
sum.aq<-summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=0.95, expandtomarkers=TRUE)
chr.aq<-{}
pos.aq<-{}
for(i in 1:nrow(sum.aq$lod)){
chr.aq1<-sum.aq$lod[i,1]
chr.aq2<-as.numeric(as.character(chr.aq1))
chr.aq<-c(chr.aq, chr.aq1)
}
if(is.na(chr.aq)==TRUE) { chr.aq<-chr } else { chr.aq<-c(chr, chr.aq) }

for(i in 1:nrow(sum.aq$lod)){
pos1.aq<-sum.aq$lod[i,2]
pos.aq<-c(pos.aq, pos1.aq)
}
if(is.na(pos.aq)==TRUE) { pos.aq<-pos } else { pos.aq<-c(pos,pos.aq) }

qtl<-makeqtl(cross, chr=chr.aq, pos=pos.aq, what="prob")
createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
##########################################################################################

#(28) Use forward selection and backward elimination model selection to probe the model space 
#for the best fit QTL model explaining your data.
#NOTE: If outliers were not removed, stepwiseqtl migth give artificial, adjacent QTLs with opposite effects

##########################################################################################
pen<-summary(sub.perms)
cross.sw<-stepwiseqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", penalties=pen, model="normal", additive.only=TRUE)
swQTL<-capture.output(print(cross.sw))
swQTL2<-c(pheno, swQTL)
cat(swQTL2, file="17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
cat(space, file="17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
##########################################################################################

#(29) Add additional QTL (if any) to the qtl object found by stepwise model selection.

##########################################################################################
	sum.sw<-summary(cross.sw)
if (length(sum.sw)==0){
	next
}
else{
	chr.sw<-{}
	pos.sw<-{}
	for(i in 1:nrow(sum.sw)){
	chr.sw1<-sum.sw[i,2]
	chr.sw2<-as.numeric(as.character(chr.sw1))
	chr.sw<-c(chr.sw, chr.sw2)
	}
	for(i in 1:nrow(sum.sw)){
	pos1.sw<-sum.sw[i,3]
	pos.sw<-c(pos.sw, pos1.sw)
	}

qtl2<-makeqtl(cross, chr=chr.sw, pos=pos.sw, what="prob")
createqtl<- paste("Q", 1:qtl2$n.qtl, sep="")
formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
rqtl<-refineqtl(cross, pheno.col=phenotype, qtl=qtl2, method="hk", model="normal")
}
##########################################################################################
#(30) write a file with the LOD scores of every marker of the model part. ONLY for chromosome with QTLs

##########################################################################################
model_lods <- list()
  lodprof <- attr(rqtl,"lodprofile")
  
  for (i in 1:length(lodprof)){
    model_lods[[i]] <- lodprof[[i]]
  }
model_lods.df <- as.data.frame(do.call(rbind, model_lods))
file<-paste("18_LOD scores at every marker of model part ", pheno, ".txt")
write.table(model_lods.df, file=file, append=TRUE, col.names=NA, sep="\t")
##########################################################################################
#(31) Write a file containing the peak marker and the flanking markers representing the 1.5 LOD interval of each QTL.

##########################################################################################
Q<-"Q"
space<-" "
for (i in 1:rqtl$n.qtl){
	interval<-capture.output(lodint(rqtl, qtl.index=i, drop=1.5, expandtomarkers=TRUE))
	q<-paste(Q, i, sep="")
	interval.new<-c(pheno, q, interval, space)
	cat(interval.new, file="19_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
	}
##########################################################################################

#(32) write a csv file containing the results of ANOVA for the full and reduced models, 
#the % variance explained by each QTL and the estimated effect size 
#(half the distance between the means for each genotype class in the case of RILs).

#*NOTE*# You may receive a warning here about dropping individuals with missing phenotypes. 
#This is expected if such a case exists and does not effect the output of the code.

##########################################################################################
cross.ests<-fitqtl(cross, pheno.col=phenotype, qtl=rqtl, formula=formula, method="hk", dropone=TRUE, get.ests=TRUE, model="normal")
ests<-capture.output(summary(cross.ests))
ests<-c(pheno, ests)
write(ests, file="20_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
##########################################################################################

#(32) Generate PDF of the effect plots and .txt with means and standard error for each genotype class. 
#Red circles represent imputed data by sim.geno (step 23)

#*NOTE*# This will likely produce a warning indicating that column names are being appended to file.

##########################################################################################
file<-paste("21_Marker effect plots", pheno, ".pdf")
pdf(file=file, width=11, height=8.5)
for(i in 1:length(chr.sw)) {
	b<-paste("Q",i, sep="")
	
	mar<-find.marker(cross2, chr=chr.sw[i], pos=pos.sw[i])
	
	plot.pxg(cross, marker=mar, pheno.col=phenotype)
	
	phenoqtl<-paste(pheno, b)
	pheno.eff<-c(phenoqtl, space)
	means<-effectplot(cross2, pheno.col=phenotype, mname1=mar, draw=FALSE)
	cat(pheno.eff, file="22_means and SE.txt", sep="\n", append=TRUE)
	write.table(means$Means, file="22_means and SE.txt", sep=",", col.names="Means", row.names=TRUE, append=TRUE)
	cat(space, file="22_means and SE.txt", sep="\n", append=TRUE)
	write.table(means$SEs, file="22_means and SE.txt", sep=",", col.names="Standard Error", row.names=TRUE, append=TRUE)
	cat(space, file="22_means and SE.txt", sep="\n", append=TRUE)
}
dev.off()
##########################################################################################

#(34) Contingency statement if no QTL are detected in the initial genome scan for a given phenotype. 

##########################################################################################
} else { null<-"	There were no LOD peaks above the threshold"
	cat(pheno, file="16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
	cat(null, file="16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
	cat(space, file="16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
	cat(pheno, file="17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
	cat(null, file="17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
	cat(space, file="17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
	cat(pheno, file="19_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
	cat(null, file="19_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
	cat(space, file="19_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
	cat(pheno, file="20_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
	cat(null, file="20_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
	cat(space, file="20_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
	cat(pheno, file="22_means and SE.txt", sep="\n", append=TRUE)
	cat(null, file="22_means and SE.txt", sep="\n", append=TRUE)
	cat(space, file="22_means and SE.txt", sep="\n", append=TRUE) }
}
##########################################################################################
