###################################
#
# Workflow to generate genetic map
#
###################################

# load required pacakges
require(ggplot2)
require(data.table)
require(qtl)
require(ASMap)
require(vioplot)
setwd("./GBS") # make sure you are already in ./GBS!
hedo <- function(x, len=10){
  #created function to see head in both directions
  print(x[1:len,1:len])
}
# specify name of the data to read in and process
# choose from relaxed, restricted01 and restricted001
dataname <- "restricted01"
resdirectory <- paste("./results/", dataname, "/", sep="")
borders <- c(0.1, 0.9)
auto_input <- T


# read in .geno file (angsd output)
geno <- fread(paste("results/", dataname, "/", dataname, ".geno", sep=""), colClasses="character")
# create vector with genome ids
ids <- c("site", paste(rep("ex", 3), seq(1:3), sep=""), paste(rep("ax", 3), seq(1:3), sep=""), paste(rep("RIL", 195), seq(1:195), sep=""))
# create vector to store infos
infos <- numeric()
# merge columns 1 and 2 to have the full scaffold name and location
geno[,V1 := .(paste(V1, V2, sep='_'))]
geno[,V2 := NULL]

# ######### --- for doGeno 11
# # check if last column is empty and delete if yes
# if(geno[1,ncol(geno)] == ""){
#   geno[,ncol(geno) := NULL]
# }
# # depending on doGeno, the genofile contains more than 1 table:
# # doGeno 11 : 3 "files" in one, four times the whole thing
# # called genotypes are in the ? one.
# # every forth column contains values in .geno format.
# unique(subset(geno,select=4))
# unique(subset(geno,select=205))
# unique(subset(geno,select=406))
# unique(subset(geno,select=607))
# unique(subset(geno,select=8))
# unique(subset(geno,select=12))
# # select columns in .geno format and save to geno
# genorows <- 4* seq(1:201)
# # save allele identity from geno
# alleles <- subset(geno, select=c(1, 2,3))
# geno <- subset(geno, select = c(1, genorows))
# # colnames must be the same as in other formats

########## --- for doGeno 3
# check if last column is empty and delete if yes
if(ncol(geno) == 205){
  geno[, 205 := NULL]
  print("set last column to NULL as it is empty")
}
ids <- c("site","major", "minor", paste(rep("ex", 3), seq(1:3), sep=""), paste(rep("ax", 3), seq(1:3), sep=""), paste(rep("RIL", 195), seq(1:195), sep=""))
########## --
names(geno) <- ids

# store some informations in the info vector
infos["fromgeno"] <- nrow(geno)
geno <- geno[!geno[,ex1 == "-1" & ex2 == "-1" & ex3 == "-1"], ] # take out rows that have -1 in all exserta columns
infos["allEXmissing"] <- infos["fromgeno"] - nrow(geno)
geno <- geno[!geno[,ax1 == "-1" & ax2 == "-1" & ax3 == "-1"], ] # take out rows that have -1 in all axillaris columns
infos["allAXmissing"] <- (infos["fromgeno"] - infos["allEXmissing"]) - nrow(geno)


# take out missing data
# only accept markers if they contain data from > 50% of RILs and
# only accept individuals that contain data from > 30% of the markers
rils <- paste(rep("RIL", 195), seq(1:195), sep="")
testo <- copy(geno[, .SD, .SDcols = rils])
#
# Find the columns wich contain less than 30% of marker information
marker_treshold <- 0.3
infos["individuals_lost_due_to_missing_data"] <- sum(apply(geno[, .SD, .SDcols = rils], 2, function(x) sum(x != -1)) < nrow(geno)*marker_treshold)
# delete the given individuals
remove <- names(which(apply(geno[, .SD, .SDcols = rils], 2, function(x) sum(x != -1)) < nrow(geno)*marker_treshold))
geno[, (remove):= NULL ]
# update rils vector (delete the deleted columns from there)
rils <- rils[! rils %in% remove]
#
# Find the rows which contain less than 50% of individuals
individual_treshold <- 0.5
infos["markers_lost_due_to_missing_data"] <- sum(apply(geno[, .SD, .SDcols = rils], 1, function(x) sum(x != -1)) < length(rils)*individual_treshold)
remove <- which(apply(geno[, .SD, .SDcols = rils], 1, function(x) sum(x != -1)) < length(rils)*individual_treshold) # store row indices to remove
testo <- copy(geno)
keep <- seq(1,nrow(geno))[!seq(1,nrow(geno)) %in% remove]
geno <- geno[keep, ]
#
rm(keep, remove, marker_treshold, individual_treshold); gc()


# find major genotype in AXILLARIS and its frequency
acols <- c("ax1", "ax2", "ax3")
amajor <- rep(NA, nrow(geno))
afreq_major <- rep(NA, nrow(geno))
aminor <- rep(NA, nrow(geno))
afreq_minor <- rep(NA, nrow(geno))
for(i in 1:nrow(geno)){
  amajor[i] <- names(sort(table(t(geno[i,.SD, .SDcols=acols])),decreasing=T)[1])
  afreq_major[i] <- as.integer(sort(table(t(geno[i,.SD, .SDcols=acols])),decreasing=T)[1])
  aminor[i] <- names(sort(table(t(geno[i,.SD, .SDcols=acols])),decreasing=T)[2])
  afreq_minor[i] <- as.integer(sort(table(t(geno[i,.SD, .SDcols=acols])),decreasing=T)[2])
}
geno <- cbind(amajor,afreq_major, aminor, afreq_minor, geno)
rm(amajor, afreq_major, aminor, afreq_minor) ; gc()

# find major genotype in EXSERTA and its frequency
ecols <- c("ex1", "ex2", "ex3")
# geno[,unique(.SD), .SDcols=ecols] # see unique column combinations
emajor <- rep(NA, nrow(geno))
efreq_major <- rep(NA, nrow(geno))
eminor <- rep(NA, nrow(geno))
efreq_minor <- rep(NA, nrow(geno))
for(i in 1:nrow(geno)){
  emajor[i] <- names(sort(table(t(geno[i,.SD, .SDcols=ecols])),decreasing=T)[1])
  efreq_major[i] <- as.integer(sort(table(t(geno[i,.SD, .SDcols=ecols])),decreasing=T)[1])
  eminor[i] <- names(sort(table(t(geno[i,.SD, .SDcols=ecols])),decreasing=T)[2])
  efreq_minor[i] <- as.integer(sort(table(t(geno[i,.SD, .SDcols=ecols])),decreasing=T)[2])
}
geno <- cbind(emajor,efreq_major, eminor, efreq_minor, geno)
rm(emajor, efreq_major, eminor, efreq_minor) ; gc()

# more filtering
# take out cases where amajor != emajor to have a smaller dataset
geno <- geno[geno[,amajor] != geno[,emajor],]
infos["majorNOTsameasminor"] <- nrow(geno)
# take out cases where ex or ax major is heterozygote "1"
infos["AXmajorisHET"] <- nrow(geno[amajor == "1"])
infos["EXmajorisHET"] <- nrow(geno[emajor == "1"])
infos["EXmajorisNA"] <- nrow(geno[emajor == "-1"])
infos["AXmajorisNA"] <- nrow(geno[amajor == "-1"])
geno <- geno[amajor != "1"]
geno <- geno[emajor != "1"]
geno <- geno[emajor != "-1"]
geno <- geno[amajor != "-1"]

# major and minor genotype frequency among all RIL population
rilfreq_ex <- rep(NA, nrow(geno))
rilfreq_ax <- rep(NA, nrow(geno))
rilfreq_het <- rep(NA, nrow(geno))
rilfreq_NA <- rep(NA, nrow(geno))
control <- rep(NA, nrow(geno))
for(i in 1:nrow(geno)){
  rilfreq_NA[i] <- length(which(geno[i, .SD, .SDcols = rils] == "-1"))
  rilfreq_ex[i] <- length(which(geno[i, .SD, .SDcols = rils] == geno[,emajor][i])) / (length(rils) - rilfreq_NA[i])
  rilfreq_ax[i] <- length(which(geno[i, .SD, .SDcols = rils] == geno[,amajor][i])) / (length(rils) - rilfreq_NA[i])
  rilfreq_het[i] <- length(which(geno[i, .SD, .SDcols = rils] == "1")) / (length(rils) - rilfreq_NA[i])
  control[i] <- round(sum(rilfreq_ex[i], rilfreq_ax[i], rilfreq_het[i]), digits=10)
  
  if(control[i] != 1 & rilfreq_ax[i] != rilfreq_ex[i] & !(geno[,emajor][i] %in% c("1", "-1")) & !(geno[,amajor][i] %in% c("1", "-1"))){
      stop("there must be an error, allele frequencies do not sum to 1")
  }
}
geno <- cbind(rilfreq_NA, rilfreq_ex, rilfreq_het, rilfreq_ax, geno)
rm(control, rilfreq_NA, rilfreq_ex, rilfreq_het, rilfreq_ax) ; gc()

# calc some stats:
infos["ex3_ax3"] <- nrow(geno[afreq_major == 3 & efreq_major == 3,]) 
infos["ex3_ax2NA"] <- nrow(geno[efreq_major == 3 & afreq_major == 2 & aminor == -1])   
infos["ex2NA_ax3"] <- nrow(geno[efreq_major == 2 & afreq_major == 2 & eminor == -1])  
infos["ex2NA_ax2NA"] <- nrow(geno[efreq_major == 2 & afreq_major == 2 & eminor == -1 & aminor == -1])  

# FILTER : take the cases above (in "infos" vector)
a <- geno[afreq_major == 3 & efreq_major == 3,]
b <- geno[efreq_major == 3 & afreq_major == 2 & aminor == -1]
c <- geno[efreq_major == 2 & afreq_major == 2 & eminor == -1]
d <- geno[efreq_major == 2 & afreq_major == 2 & eminor == -1 & aminor == -1]

geno <- rbind(a,b,c,d)
rm(a); rm(b); rm(c); rm(d); gc()


# POLARISATION 

# The exserta alleles and the axillairs alleles must be named accordingly.
# Exserta alleles are named 'B'
# axillaris alleles are named 'A'
for(col in rils) set(geno, i=which(geno[[col]] == geno[,emajor]), j=col, value = 'B')
for(col in rils) set(geno, i=which(geno[[col]] == geno[,amajor]), j=col, value = 'A')
for(col in rils) set(geno, i=which(geno[[col]]=='1'), j=col, value='H')
for(col in rils) set(geno, i=which(geno[[col]]=='-1'), j=col, value='U')
# unique(geno[,ax2]) # control the polarisation process


# calc allele frequencies among RILS
# count number of homozygous_ex * 2 + heterozygous to get number of exserta alleles
allele_frequencies <- data.table(ex=rep(100, nrow(geno)), ax=rep(100, nrow(geno)), nas=rep(100, nrow(geno)))

for(i in 1:nrow(allele_frequencies)){
  allele_frequencies[i,"ex"] <- (length(which(geno[, .SD, .SDcols = rils][i,] == "B")) * 2 + length(which(geno[, .SD, .SDcols = rils][i,] == "H"))) / (length(rils) * 2)
  allele_frequencies[i,"ax"] <- (length(which(geno[, .SD, .SDcols = rils][i,] == "A")) * 2 + length(which(geno[, .SD, .SDcols = rils][i,] == "H"))) / (length(rils) * 2)
  allele_frequencies[i,"nas"] <- (length(which(geno[, .SD, .SDcols = rils][i,] == "U")) * 2) / (length(rils) * 2)
}


# TODO : stuck here with copying to masterfile!
# -------------METADATA TO SAVE ------------------------------------------------------------------------------------------

genoinfo <- copy(geno)
genoinfo[,(rils):= NULL]

par(ps = 12, cex = 0.9)

# ALLELE FREQUENCIES
pdf(paste("./results/",dataname, "/allele_frequencies.pdf", sep=""))
plot(density(allele_frequencies[,ex]), col="red",main=paste("Allele Frequency (" , dataname, ")"), xlab="allele frequency",
     sub = paste("missing (median) : ", median(allele_frequencies[,nas])),
     xlim=c(0,1), ylim = c(0, 3.5),
     lwd=4, font=2)
lines(density(allele_frequencies[,ax]), lwd=4)
# lines(density((1-allele_frequencies[,nas])-allele_frequencies[,ax])) # control : ex + ax + nas = 1
legend(0.8, 3.5, legend= c("ex", "ax"), col=c("red", "black"), pch=16)
dev.off()
# lines(density(allele_frequencies[,nas]), col = "lightgray", lwd=4)
# legend(0.77, 18, legend= c("ex", "ax", "NA"), col=c("red", "black", "lightgray"), pch=16)


# GENOTYPE FREQUENCIES
df = data.frame(genoinfo[,.(rilfreq_ex, rilfreq_het, rilfreq_ax)])
df.m <- reshape2::melt(df, id.vars = NULL)
pdf(paste("./results/",dataname, "/genotype_frequencies_violin.pdf", sep=""))
ggplot(df.m, aes(x = variable, y = value))  + geom_jitter(height = 0, width = 0.4, colour = "darkgray") + geom_violin(draw_quantiles = 0.5) +
  labs(x = "Genotype", subtitle = paste("homozygous ex = ", median(genoinfo[,rilfreq_ex]), "; heterozygous = " , 
                                                         median(genoinfo[,rilfreq_het])," ; \n homozygous ax = ", median(genoinfo[,rilfreq_ax]), 
                                                         "; missing (taken out) = ", median(genoinfo[,rilfreq_NA]) / length(rils)),
       y = "Frequency", title = paste("Genotype frequencies (", dataname, ")")) +
  theme(text = element_text(size = 10))
dev.off()



# RILS : how are the alleles distributed per RIL? ---------
strange <- as.data.frame(geno[which(allele_frequencies[,ex] > 0.8)])
hedo(strange)
strange2 <- as.data.frame(geno[which(allele_frequencies[,ax] < 0.2)])
# plot(strange2[,"rilfreq_ex"], strange2[,"rilfreq_ax"], pch=16,
#      col="gray", main = "Genotype frequencies of special allele frequency cases" , sub = "allele frequencies ex > 0.8 and ax < 0.2",
#      xlab = "exserta genotype frequency", ylab = "axillairs gentoype frequency", ylim = c(0,1), xlim = c(0,1))
# points(strange[,"rilfreq_ex"], strange[,"rilfreq_ax"], pch=16, col="red")
# x <- seq(0,1,0.1)
# points(x, 1-x, col="gray") # expected behaviour
# points(strange[,"rilfreq_ex"], strange[,"rilfreq_NA"]/ 193 + strange[,"rilfreq_ax"], col="green", pch = 16)
# 
# # The high allele freqeuncies are not due to many heterozygous (e.g. only hetero and exex)
# vioplot(strange[,"rilfreq_het"], strange2[,"rilfreq_het"])
# # there are many NAs, but still enough rils left to have a meaningful value.
# vioplot(strange[,"rilfreq_NA"], strange2[,"rilfreq_NA"])

# how many axillaris absolutes are left?
i <- seq(1,nrow(strange))
abs_ex <- apply(strange[i, rils] == "B", 1, sum)
abs_ax <- apply(strange[i, rils] == "A", 1, sum)
plot(abs_ex + strange[,"rilfreq_NA"] , abs_ax, pch = 16, col="darkgray",
     main = "Absolute genotype occurrences", ylab = "# axillaris genotype", xlab="# exserta genotype + # missing genotype",
     sub="subset : allele frequency ax < 0.2 and allele frequency ex > 0.8", xlim = c(0,188), ylim = c(0,188))
points(abs_ax + strange[,"rilfreq_NA"], abs_ex, pch = 16, col = "black")
# TODO find the logic flow of that again and filter accordingly.

# --------------------------------------------------------------------------------------------------IMPORTANT

if(auto_input){
  stop("take out some rows according to last plot")
}

# take out the cases where sequencing errors are probable (# ax < 10)
# geno <- geno[!(which(allele_frequencies[,ax] < 0.01)),]# [inds]

# # --------------------------------------------------------------------------------------------------------------------------------
# # ------------------------------------check after group meeting (if keep or not)
# # plot some stats
# # plot the marker frequency among all markers
# genoinfo <- copy(geno)
# genoinfo[,(rils):= NULL]
# hist(genoinfo[,rilfreq_ex], breaks = 1000)
# plot(density(genoinfo[,rilfreq_ex]), col="red",main=paste("allele frequency (" , dataname, ")"), xlab="allele/genotype? frequency")
# lines(density(genoinfo[,rilfreq_ax]))
# 
# barplot(c(median(genoinfo[,rilfreq_ex]), median(genoinfo[,rilfreq_het]), median(genoinfo[,rilfreq_ax])), names.arg=c("ex", "het", "ax"),
#         main="Genotype frequencies of RIL7 population", ylim=c(0,1))
# 
# stripchart(list("rilfreq_ex" = genoinfo[,rilfreq_ex], "rilfreq_het"= genoinfo[,rilfreq_het], "rilfreq_ax" = genoinfo[,rilfreq_ax]),
#            vertical=T, method="overplot", pch=19, offset = 0.8)
# 
# 
# Take out markers with genotype frequencies which are very high or very low (see vector borders).
# borders are defined at the beginning of the document
print("taking out very high and very low frequency markers. Borders : ") ; borders
genoinfo <- genoinfo[rilfreq_ex >= borders[1] & rilfreq_ex <= borders[2] & rilfreq_ax >= borders[1] & rilfreq_ax <= borders[2],]
geno <- geno[rilfreq_ex >= borders[1] & rilfreq_ex <= borders[2] & rilfreq_ax >= borders[1] & rilfreq_ax <= borders[2],]

# after cutting out all above/below borders
pdf(paste("./results/",dataname, "/allele_frequencies_after_cutting.pdf", sep=""))
plot(density(allele_frequencies[,ex]), col="red",main=paste("Allele Frequency (" , dataname, ")"), xlab="allele frequency",
     sub = paste("missing (median) : ", median(allele_frequencies[,nas])),
     xlim=c(0,1), ylim = c(0, 3.5),
     lwd=4, font=2)
lines(density(allele_frequencies[,ax]), lwd=4)
# lines(density((1-allele_frequencies[,nas])-allele_frequencies[,ax])) # control : ex + ax + nas = 1
legend(0.8, 3.5, legend= c("ex", "ax"), col=c("red", "black"), pch=16)
dev.off()


# 
# # Plots after taking out very high and very low genotype frequencies.
# hist(genoinfo[,rilfreq_ex], breaks = 1000, xlim = c(0,1))
# plot(density(genoinfo[,rilfreq_ex]), col="red",main=paste("allele frequency (" , dataname, ")"), xlab="allele/genotype? frequency",
#      xlim=c(0,1))
# lines(density(genoinfo[,rilfreq_ax]))
# 
# barplot(c(median(genoinfo[,rilfreq_ex]), median(genoinfo[,rilfreq_het]), median(genoinfo[,rilfreq_ax])), names.arg=c("ex", "het", "ax"),
#         main="Genotype frequencies of RIL7 population", ylim=c(0,1))
# 
# stripchart(list("rilfreq_ex" = genoinfo[,rilfreq_ex], "rilfreq_het"= genoinfo[,rilfreq_het], "rilfreq_ax" = genoinfo[,rilfreq_ax]),
#            vertical=T, method="overplot", pch=19, offset = 0.8)
# 

# -------------------------------------------------------------------------------------------------------------------------------------------
# take out all info columns from geno for later read-in to ASMap
geno <- geno[,.SD, .SDcols = c("site", rils)]

# store amount of markers left in the end
infos["markers_left_final"] <- nrow(geno)
infos["individuals_left_final"] <- ncol(geno) -1
# save the informations vector
write.table(infos, paste("./results/",dataname, "/infos.txt", sep=""))

## bring into CSV format for read.cross
tgeno <- t(geno)
colnames(tgeno) <- tgeno[1,]
tgeno <- tgeno[-1,]
# create the first row of chromosome information (fake) and a column of phenotypes (fake)
chrom <- rep(1,ncol(tgeno)) # better to use rep than seq, as there is a for loop over all
    # chromosomes in the package code of read cross.
tgeno <- rbind(chrom, tgeno)

pheno <- c("", seq(1,length(tgeno[,2])-1))
tgeno <- cbind(pheno, tgeno)
tgeno <- rbind(colnames(tgeno), tgeno)
rownames(tgeno) <- NULL
colnames(tgeno) <- NULL
write.csv(tgeno, file = paste(resdirectory, "readcross_format_", dataname, ".csv", sep=""), na = "", row.names = FALSE)

# delete all the things I don't need any more (only keep geno)
rm(chrom); rm(col); rm(ids); rm(pheno); rm(tgeno); rm(ecols); rm(strange2)
rm(i); rm(rils); rm(abs_ax); rm(abs_ex) ; rm(acols); rm(infos); rm(tres)
gc()



# ------------------------------------------------------------------------
# Genetic map formation with ASMap
# read in the data
# read.cross is far too slow! edit it!
source(paste(resdirectory, "read.cross.fast.R", sep=""))
# error.prob from genotype calling could be given
gtypes <- read.cross(format="csv", file=paste(resdirectory, "readcross_format_", dataname, ".csv", sep=""), na.strings="U", F.gen = 7, BC.gen=0)
# gtypes$geno$`1`$data[1:10,1:10] # show the genotype data
# gtps <- convertbcsft(gtypes, F.gen=7, BC.gen=0) # use this to specify F.gen and BC.gen, this is required for RILn populations

# try to plot RIL genotypes
geno.image(gtypes, sub = "genotypes AX (red), HET (blue), EX (green)")
# pull out uninformative markers
sum(nmar(gtypes))
gtypes2 <- pullCross(gtypes, type = "co.located")
sum(nmar(gtypes2))


# generate genetic maps
map1 <- mstmap.cross(gtypes2, pop.type = "RIL7", as.cross=TRUE, bychr = TRUE, trace = FALSE, p.value = 1e-10, id="pheno")
plot(map1)
heatMap(map1, lmax=20)
nmar(map1)

# # bychr = F, trace = T
# # meaning : 
# map2 <- mstmap(gtypes2, pop.type = "RIL7", as.cross=T, p.value = 1e-15, id="pheno", bychr = F, dist.fun = "kosambi", trace = TRUE)
# plot(map2)d
# heatMap(map2, lmax=60)
# nmar(map2)

saveRDS(map1, paste(resdirectory, "geneticmap_", dataname, ".rds", sep=""))
# map1 <- readRDS("geneticmap_restricted.rds")

# manually split linkage groups
# ...



############### august test #######################
# Exserta alleles are named 'B'
# axillaris alleles are named 'A'
aug <- read.cross(format = "csv", file = "results/restricted01/readcross_format_restricted01.csv", F.gen = 7, genotypes = c("A", "H", "B"), na.strings = "U")#, crosstype = "riself")
aug <- convert2riself(aug)
f7_K <- pullCross(aug, type = "co.located")
map1 <- mstmap.cross(f7_K, bychr = F, dist.fun = "kosambi", trace = FALSE, detectBadData = F, p.value = 1e-09, mvest.bc = T, return.imputed = T)

# order markers within linkage groups
map1 <- mstmap.cross(map1, bychr = T, dist.fun = "kosambi", trace = FALSE, detectBadData = F, p.value = 1e-09, mvest.bc = F, return.imputed = T)
summary(map1)
heatMap(map1, lmax=15)
profileMark(map1, stat.type = c("seg.dist", "dxo", "erf", "lod"), id = "Genotype", layout = c(1, 4), type = "l")
mean(countXO(map1))
map1 <- pushCross(map1, type="co.located")
# order again
map1 <- mstmap(map1, bychr = T, dist.fun = "kosambi", trace = TRUE, detectBadData = F, p.value = 1e-09, mvest.bc = F, return.imputed = T)
heatMap(map1, lmax=15)
profileMark(map1, stat.type = c("seg.dist", "dxo", "erf", "lod"), id = "Genotype", layout = c(1, 4), type = "l")
map1 <- cross2int(map1, id="Genotype", rem.mark=F) # rem.mark = F to not take out colocated markers from the map
map1 <- mergeCross(map1, merge = list("L.8" = c("L.8", "L.9", "L.10", "L.11" ,"L.12")))
write.cross(map1, format="qtab", filestem="geneticmap")
# 'qtab' format produces a lot of files I will not further use. Delete them. The only file saved will be 'geneticmap_location.qtab'
file.remove("geneticmap_phenotypes.qtab"); file.remove("geneticmap_founder.qtab"); file.remove("geneticmap_genotypes.qtab"); file.remove("geneticmap_symbols.qtab")

# not translated to optical map genome, because case seemed to be lost.