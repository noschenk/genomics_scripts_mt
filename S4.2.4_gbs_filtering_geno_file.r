###################################
#
# Filtering data from geno file
#
###################################

# ---------- READ IN & PREPARE DATASET ----------------------------------------------------
# read in .geno file (angsd output)
geno <- fread(paste("results/", dataname, "/", dataname, ".geno", sep=""), colClasses="character")
# create vector with genome ids
ids <- c("site", paste(rep("ex", 3), seq(1:3), sep=""), paste(rep("ax", 3), seq(1:3), sep=""), paste(rep("RIL", 195), seq(1:195), sep=""))
# create vector to store infos
infos <- numeric()
# merge columns 1 and 2 to have the full scaffold name and location
geno[,V1 := .(paste(V1, V2, sep='_'))]
geno[,V2 := NULL]

# check if last column is empty and delete if yes
if(ncol(geno) == 205){
  geno[, 205 := NULL]
  print("set last column to NULL as it is empty")
}
ids <- c("site","major", "minor", paste(rep("ex", 3), seq(1:3), sep=""), paste(rep("ax", 3), seq(1:3), sep=""), paste(rep("RIL", 195), seq(1:195), sep=""))

names(geno) <- ids

# store some informations in the info vector
infos["fromgeno"] <- nrow(geno)
geno <- geno[!geno[,ex1 == "-1" & ex2 == "-1" & ex3 == "-1"], ] # take out rows that have -1 in all exserta columns
infos["allEXmissing"] <- infos["fromgeno"] - nrow(geno)
geno <- geno[!geno[,ax1 == "-1" & ax2 == "-1" & ax3 == "-1"], ] # take out rows that have -1 in all axillaris columns
infos["allAXmissing"] <- (infos["fromgeno"] - infos["allEXmissing"]) - nrow(geno)
# -----------------------------------------------------------------------------------------


# ----------- TAKE OUT MISSING DATA -------------------------------------------------------
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
# -----------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------


# -----FILTER FOR MAJOR GENOTYPES ---------------------------------------------------------
# take out cases where amajor != emajor
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
# -----------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------


# ---------POLARISATION-----------------------------------------------------------------------
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
# -----------------------------------------------------------------------------------------