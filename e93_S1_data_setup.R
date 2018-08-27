#------------------------------------#
#------~~~~~ SECTION 1 ~~~~~---------#
#-----~~~~ E93 DATA SET-UP ~~~~------#
#------------------------------------#

# Authors: Annabel Smith & Jane Catford

# load function library:
source("/Users/annabelsmith/Documents/01_Current/PROJECTS/03_CATFORD_invasion_traits/DATA/AS_data_and_analysis/ANALYSIS_files/e93_script/e93_S0_function_library.R")

load("e93_wksp_data")

######################################
#       	PROCESS RAW DATA		 #
######################################
{
# Species names must be in quotes. If excel replaces one quote with three, find and replace with a single quote in a text editor before importing. The text file names match the sheet names in the excel files. Check that the notes columns (orange in excel) have been deleted.

# Data dir:
dat.dir<-"../e93_data"
dir(dat.dir)

# Species and enviro data:
{
sp_dat<-read.table(paste(dat.dir,"E93_species_abund.txt",sep="/"),header=T)
sp_info<-read.table(paste(dat.dir,"species_info.txt",sep="/"),header=T)
no_pl<-read.table(paste(dat.dir,"E93_no_planted.txt",sep="/"),header=T)
seeded<-read.table(paste(dat.dir,"E93_seeded_species.txt",sep="/"),header=T)
lcd_names<-read.table(paste(dat.dir,"LCD_names.txt",sep="/"),header=T)
topog<-read.table(paste(dat.dir,"E93_elevation.txt",sep="/"),header=T)
soil<-read.table(paste(dat.dir,"E93_soil.txt",sep="/"),header=T)
burn_precip<-read.table(paste(dat.dir,"E93_burn_precip.txt",sep="/"),header=T)
soil_dist<-read.table(paste(dat.dir,"E93_soil_disturbance.txt",sep="/"),header=T)

head(sp_dat); head(sp_info); head(no_pl); head(seeded); head(lcd_names); head(topog); head(soil); head(burn_precip); head(soil_dist)
} # close species & enviro data

# Trait data:
{
hgt_dat<-read.table(paste(dat.dir,"E93_E275_height.txt",sep="/"),header=T)
lf_dat<-read.table(paste(dat.dir,"E93_E275_leaf.txt",sep="/"),header=T)
seed_dat<-read.table(paste(dat.dir,"E93_E275_seed_mass.txt",sep="/"),header=T)
e133_reps<-read.table(paste(dat.dir,"E133_leaf_reps.txt",sep="/"),header=T)
e133_avg<-read.table(paste(dat.dir,"E133_leaf_avg.txt",sep="/"),header=T)
gaps<-read.table(paste(dat.dir,"trait_gaps.txt",sep="/"),header=T)

head(hgt_dat); head(lf_dat); head(seed_dat); head(e133_reps); head(e133_avg); head(gaps)
} # close trait data

# E275 DATA:
# These aren't used for the E93 analysis but are included here for processing:
{
dsize<-read.table(paste(dat.dir,"E275_dispersule_size.txt",sep="/"),header=T)
dshape<-read.table(paste(dat.dir,"E275_dispersule_shape.txt",sep="/"),header=T)
root<-read.table(paste(dat.dir,"E275_root_SRL.txt",sep="/"),header=T)
rootMD<-read.table(paste(dat.dir,"E275_root_max_depth.txt",sep="/"),header=T)
plantDM<-read.table(paste(dat.dir,"E275_plant_dry_mass.txt",sep="/"),header=T)

head(dsize); head(dshape); head(root); head(rootMD); head(plantDM)

} # close E275 data

#### ------ LCD names 
{
# Update names to LCDs:

# After reading data, enter names of data sets which need to have names updated to LCDs. Also read in the LCD data file (lcd_names) but don't include it in this line, only the data sets you want to update:
data.sets<-c("sp_dat","sp_info","seeded","hgt_dat","lf_dat","seed_dat","e133_reps","e133_avg","dsize","dshape","root","rootMD","plantDM")

# Reduce LCD data to necessary columns (i.e. keep only the CDR cols):
lcd<-lcd_names[,which(names(lcd_names) %in% c("CDR","CDR_LCD","CDR_abbrev","CDR_LCD_abbrev"))]
lcd<-lcd[order(lcd$CDR),]
lcd<-data.frame(apply(lcd,2,as.character),stringsAsFactors=F)
lcd$CDR_LCD[which(is.na(lcd$CDR_LCD))]<-lcd$CDR[which(is.na(lcd$CDR_LCD))]
lcd$match<-lcd$CDR==lcd$CDR_LCD
lcd$match_abb<-lcd$CDR_abbrev==lcd$CDR_LCD_abbrev

# Subset LCD on species which need to be updated:
lcd<-lcd[-which(lcd$match==T & lcd$match_abb==T),]
lcd<-tidy.df(lcd)

# Run the LCD update loop:
for (i in 1:length(data.sets)){

data.thisrun<-get(data.sets[i])
data.thisrun$species<-as.character(data.thisrun$species)
data.thisrun$sp<-as.character(data.thisrun$sp)

# The name changes are going to eliminate the uniqueness of some samples. E.g. in the height data, Pan_lan will become Pan_pra and we will have two sets of Pan_pra A1:A5. Keep the old name so these instances can be identified if necessary:
data.thisrun$CDR_name<-data.thisrun$species

# Get index of species to be updated and compile table of name changes specific to the current data set (data.sets[i]):
ind.update<-which(data.thisrun$species %in% lcd$CDR==T)
sp.update<-data.thisrun$species[ind.update]
sp.upd<-data.frame(CDR=unique(sp.update),stringsAsFactors=F)
sp.upd<-merge(sp.upd,lcd,by="CDR",all.x=T,all.y=F)

# Update species names:
data.thisrun$species[ind.update]<-as.character(sapply(data.thisrun$species[ind.update],function(x){sp.upd$CDR_LCD[match(x,sp.upd$CDR)]}))

# Update species abbreviations:
data.thisrun$sp[ind.update]<-as.character(sapply(data.thisrun$sp[ind.update],function(x){sp.upd$CDR_LCD_abbrev[match(x,sp.upd$CDR_abbrev)]}))

# Re-factorise species, sp and CDR_name:
data.thisrun$species<-as.factor(data.thisrun$species)
data.thisrun$sp<-as.factor(data.thisrun$sp)
data.thisrun$CDR_name<-as.factor(data.thisrun$CDR_name)

# Assign back to original data name:
assign(data.sets[i],data.thisrun)

} # close i LCD

} # close LCD names

#### ------ SPECIES INFO:
{
# sp_info now has duplicated rows because some were grouped as LCDs. Remove duplicates: 
# sp_info<-sp_info[,1:(length(sp_info)-1)]
head(sp_info)

# Check to see if all life form data, etc. are the same for the LCDs:
dup.sp<-data.frame(species=unique(sp_info[which(duplicated(sp_info$species)),]$species))
dup.sp<-merge(dup.sp,sp_info,by="species",all.x=T,all.y=F)

# Remove duplicates
sp_info<-sp_info[-which(duplicated(sp_info$species)),]
sp_info<-tidy.df(sp_info)

# All life form data are the same for the LCDs, except for Lactuca sp. L. canadensis=native; L. sp.='unknown'. JC advised to use unknown. 
sp_info$origin[sp_info$species=="Lactuca sp."]<-"Unknown"
sp_info<-sp_info[,1:(length(sp_info)-1)]
sp_info<-tidy.df(sp_info)

# Update genus names:
sp_info$genus<-as.character(sp_info$genus)
for (i in 1:length(sp_info[,1])){
sp_info$genus[i]<-substr(sp_info$species[i],1,gregexpr(" ",sp_info$species[i])[[1]][1]-1)
}
sp_info$genus<-as.factor(sp_info$genus)
head(sp_info)

# Make grass and legume categories:
levels(sp_info$lifeform)
check.rows(sp_info,4)

sp_info$grass<-0
sp_info$grass[sp_info$lifeform=="Grass"]<-1

sp_info$legume<-0
sp_info$legume[sp_info$lifeform=="Legume"]<-1

} # close species info

#### ------ PROCESS SOIL DATA:
{
# Mean soil variables at the subplot level (i.e. the level used for analysis):
soil$plsub<-paste(soil$plot, soil$subplot, sep="")

# calculate means:
soil.mean<-data.frame(plsub=names(tapply(soil$moisture, soil$plsub, mean)), moisture=as.numeric(tapply(soil$moisture, soil$plsub, mean)), NO3soil=as.numeric(tapply(soil$NO3soil, soil$plsub, mean)), NH4soil=as.numeric(tapply(soil$NH4soil, soil$plsub, mean)))
rownames(soil.mean)<-1:length(soil.mean[,1])

# Add plot 1 data (mean of all other plots):
levels(soil.mean$plsub)<-c(levels(soil.mean$plsub),"1W","1X","1Y","1Z")
pl1<-data.frame(plsub=c("1W","1X","1Y","1Z"),moisture=mean(soil$moisture), NO3soil=mean(soil$NO3soil), NH4soil=mean(soil$NH4soil))
soil.mean<-rbind(soil.mean,pl1)

# add plot and subplot data and tidy up:
soil.mean$plot<-as.numeric(substr(soil.mean$plsub,1, unlist(gregexpr("[A-Z]", soil.mean$plsub))-1))
soil.mean$subplot<-substr(soil.mean$plsub, unlist(gregexpr("[A-Z]", soil.mean$plsub)),unlist(gregexpr("[A-Z]", soil.mean$plsub)))
soil.mean<-soil.mean[order(soil.mean$plot,soil.mean$subplot),]
soil.mean<-tidy.df(soil.mean)
soil.mean$plsub<-as.factor(paste(soil.mean$plot, soil.mean$subplot, sep=""))

head(soil.mean)
tail(soil.mean)
} # close soil data

#### ------ PROCESS HEIGHT 
{
hgt_dat<-hgt_dat[,1:(length(hgt_dat)-1)]

# Exclude juveniles and babies:
hgt_dat$age<-substr(hgt_dat$rep,1,1)
hgt_dat<-hgt_dat[hgt_dat$age=="A",]
hgt_dat<-tidy.df(hgt_dat)

# There are no species that have > 1 condition, so no need to calc height means by condition.
# Lol_mul and Cas_fas are all M, everything else is P
table(hgt_dat$sp,hgt_dat$cond)

# Calculate mean and CV removing NAs:
hgt<-data.frame(sp=levels(hgt_dat$sp),hgt_mean=round(tapply(hgt_dat$height,hgt_dat$sp,mean,na.rm=T),2),hgt_sd=round(tapply(hgt_dat$height,hgt_dat$sp,sd,na.rm=T),2))
hgt$hgt_cv<-round(hgt$hgt_sd/hgt$hgt_mean,2)
row.names(hgt)<-1:length(hgt[,1])
# Remove sd before merging:
hgt$hgt_sd<-NULL
head(hgt,20)

# Add height means to species info:
sp_info<-merge(sp_info,hgt,by="sp",all.x=T)

# There are four species that are included in the height data file but have no mean: Ama_ret, Set_ita, Set_vir (adults and babies), Sis_alt. Most of these (except Set_vir which has NAs for adults) were all babies and are indicated as 'not viable' in the original data. Their means will come out as NA, so shouldn't affect the processing if they're left in. 
# sp_info[as.numeric(which(is.na(sp_info[,15])==T)),]
} # close height

#### ------ PROCESS SEED MASS 
{
seed_dat<-seed_dat[,-which(colnames(seed_dat)=="CDR_name")]

# Mass per seed:
seed_dat$mps<-seed_dat$mass/seed_dat$no_seeds
head(seed_dat)

# FOR THE SPECIES WITH > 1 SOURCE, USE A SINGLE SOURCE FOR THE OVERALL MEANS:
seedsp<-data.frame(sp=levels(seed_dat$sp))
seedsp$no.sources<-NA
head(seedsp)
for(i in 1:length(levels(seed_dat$sp))){
seedsp$no.sources[i]<-length(unique(seed_dat[seed_dat$sp==seedsp$sp[i],3]))
}
sources<-seedsp[which(seedsp$no.sources>1),]
sources<-droplevels(sources)

# Seed data for all species with more than one source (9 species total):
sd_mg<-merge(sources,seed_dat)
sd_mg<-droplevels(sd_mg)
head(sd_mg)

# Add number of sources to seed_dat:
sdmg<-merge(seed_dat,seedsp)
seed_dat<-sdmg
head(seed_dat)

# Remove the > 1 source data from the main seed data frame:
seed_dat<-seed_dat[-which(seed_dat$no.sources==2),]

# From NOTES in E93_E275_height_seeds.xlsx and email correspondence from JC. Remove one source for species with > 1 source. Use Troy/Field/Adam seed over Prairie moon as the former were all collected at Cedar Creek and the latter in the region but not locally. For species with > 1 source, where removing the PMOON source would leave < 3 samples, use all seed samples. The latter conditions occur in three species (Artemisia ludoviciana, Schizacharium scoparium and Stipa spartea). 

sdtb<-data.frame(table(sd_mg$species, sd_mg$source))
sd.names<-unique(sdtb$Var1[which(sdtb$Freq <3 & sdtb$Freq>0)])
odd.species<-rbind(sd_mg[sd_mg$species==sd.names[1],],sd_mg[sd_mg$species==sd.names[2],],sd_mg[sd_mg$species==sd.names[3],])

head(seed_dat)
head(sd_mg)

# Rearrange sd_mg so that the cols line up with seed_dat:
sd_mg<-sd_mg[,colnames(seed_dat)]

# Remove three odd species from sd_mg:
sd_mg<-sd_mg[-c(which(sd_mg$species==sd.names[1]),which(sd_mg$species==sd.names[2]),which(sd_mg$species==sd.names[3])),]

# Remove any remaining PMOON from sd_mg:
sd_mg<-sd_mg[sd_mg$source!="PMOON",]

# Stitch them back together:
seed_dat<-rbind(seed_dat,sd_mg,odd.species)

# 48 records in the >1 source data, 25 by PMOON. seed_dat after removal should be 456 (481-25) rows, plus seven PMOON records in odd.species = 463. 
dim(seed_dat)

seed_dat<-seed_dat[order(seed_dat$species,seed_dat$rep),]
row.names(seed_dat)<-1:length(seed_dat[,1])
head(seed_dat)
tail(seed_dat)

# Add seed mean and CV to species info (na.rm):
sd<-data.frame(sp=levels(seed_dat$sp),seed_mean=round(tapply(seed_dat$mps,seed_dat$sp,mean,na.rm=T),8),seed_sd=round(tapply(seed_dat$mps,seed_dat$sp,sd,na.rm=T),8))
sd$seed_cv<-round(sd$seed_sd/sd$seed_mean,8)
row.names(sd)<-1:length(sd[,1])
# Remove sd before merging:
sd$seed_sd<-NULL
head(sd,20)
sp_info<-merge(sp_info,sd,by="sp",all.x=T)
head(sp_info)

} # close seed mass

#### ------ PROCESS LEAF DATA 
{
lf_dat<-lf_dat[,-which(colnames(lf_dat)=="CDR_name")]

# Change "J" to "A" for Ber_inc and Ver_tha:
head(lf_dat)
lf_dat$stage[lf_dat$sp %in% c("Ber_inc","Ver_tha")]<-"A"
lf_dat[lf_dat$sp %in% c("Ber_inc","Ver_tha"),]

# Use adult data only:
lf_dat<-lf_dat[lf_dat$stage=="A",]
lf_dat<-tidy.df(lf_dat)

# Use whole leaf only:
lf_dat<-lf_dat[lf_dat$type=="whole",]
lf_dat<-tidy.df(lf_dat)

# There are no species that have > 1 condition, so no need to calc leaf means by condition.
# Cas_fas and Lol_mul are all M, everything else is P:
table(lf_dat$sp,lf_dat$conditions)

# All species come from one experiment or the other, not both, so no need to calc means by exp:
table(lf_dat$sp,lf_dat$exp)
which(table(lf_dat$sp,lf_dat$exp)[,1]>0 & table(lf_dat$sp,lf_dat$exp)[,2]>0)

# These have only one sample in leaf data. JC advised to use them as the value for these sepcies:
# Gen_and, Lia_pyc
table(lf_dat$sp)

# Add LDMC mean and CV to species info (na.rm):
ldmc.mean<-data.frame(sp=levels(lf_dat$sp),ldmc_mean=round(tapply(lf_dat$LDMC,lf_dat$sp,mean,na.rm=T),8),ldmc_sd=round(tapply(lf_dat$LDMC,lf_dat$sp,sd,na.rm=T),8))
ldmc.mean$ldmc_cv<-round(ldmc.mean$ldmc_sd/ldmc.mean$ldmc_mean,8)
row.names(ldmc.mean)<-1:length(ldmc.mean[,1])
# Remove sd before merging:
ldmc.mean$ldmc_sd<-NULL
head(ldmc.mean,20)

sp_info<-merge(sp_info,ldmc.mean,by="sp",all.x=T)
head(sp_info)

# Add SLA mean and CV to species info (na.rm):
sla.mean<-data.frame(sp=levels(lf_dat$sp),sla_mean=round(tapply(lf_dat$SLA,lf_dat$sp,mean,na.rm=T),8),sla_sd=round(tapply(lf_dat$SLA,lf_dat$sp,sd,na.rm=T),8))
sla.mean$sla_cv<-round(sla.mean$sla_sd/sla.mean$sla_mean,8)
row.names(sla.mean)<-1:length(sla.mean[,1])
# Remove sd before merging:
sla.mean$sla_sd<-NULL
head(sla.mean,20)

sp_info<-merge(sp_info,sla.mean,by="sp",all.x=T)
head(sp_info)
sp_info$hgt_mean<-as.numeric(sp_info$hgt_mean)
sp_info$hgt_cv<-as.numeric(sp_info$hgt_cv)
sp_info$seed_mean<-as.numeric(sp_info$seed_mean)
sp_info$seed_cv<-as.numeric(sp_info$seed_cv)
sp_info$ldmc_mean<-as.numeric(sp_info$ldmc_mean)
sp_info$ldmc_cv<-as.numeric(sp_info$ldmc_cv)
sp_info$sla_mean<-as.numeric(sp_info$sla_mean)
sp_info$sla_cv<-as.numeric(sp_info$sla_cv)

# E133 leaf data:

# Preferentially use averages (a more complete data set, with 128 species) then fill gaps with reps (89 species). Update later if we get our hands on the reps that were used to derive the averages (this is looking unlikely).

# update the e133_avg data to deal with the replicates arising from LDC grouping. 
e133_avg<-e133_avg[,1:(length(e133_avg)-1)]

e133_dups<-data.frame(sp=e133_avg$sp[which(duplicated(e133_avg$sp))])
e133_dups<-merge(e133_dups,e133_avg,by="sp",all.x=T,all.y=F)
e133_dups<-tidy.df(e133_dups)

e133_dups<-data.frame(species=levels(e133_dups$species),sp=levels(e133_dups$sp),fresh_mass=tapply(e133_dups$fresh_mass,e133_dups$sp,mean,na.rm=T),dry_mass=tapply(e133_dups$dry_mass,e133_dups$sp,mean,na.rm=T),area=tapply(e133_dups$area,e133_dups$sp,mean,na.rm=T),LDMC=tapply(e133_dups$LDMC,e133_dups$sp,mean,na.rm=T),SLA=tapply(e133_dups$SLA,e133_dups$sp,mean,na.rm=T),SLA_cm2_g=tapply(e133_dups$SLA_cm2_g,e133_dups$sp,mean,na.rm=T))

head(e133_avg)
e133_avg<-e133_avg[-which(e133_avg$sp %in% e133_avg$sp[which(duplicated(e133_avg$sp))]==T),]
e133_avg<-rbind(e133_avg,e133_dups)
e133_avg<-e133_avg[order(e133_avg$species),]
e133_avg<-tidy.df(e133_avg)

# Put data in from E133 averages. There are different sets of species with no LDMC and SLA, so these will need to be done separately:
noldmc<-data.frame(sp=sp_info$sp[which(is.na(sp_info$ldmc_mean))])
noldmc<-tidy.df(noldmc)

nosla<-data.frame(sp=sp_info$sp[which(is.na(sp_info$sla_mean))])
nosla<-tidy.df(nosla)

ldmc.avg<-e133_avg[,c("sp","LDMC")]
sla.avg<-e133_avg[,c("sp","SLA")]

new.ldmc<-merge(noldmc,ldmc.avg,by="sp",all.x=T, all.y=F)
new.sla<-merge(nosla,sla.avg,by="sp",all.x=T, all.y=F)

# Reduce cols in sp_info:

sp_info<-sp_info[,-grep("cv",colnames(sp_info))]
sp_info<-merge(sp_info,new.ldmc,by="sp",all.x=T, all.y=F)

sp_info$ldmc_mean<-rowSums(sp_info[,c("ldmc_mean","LDMC")],na.rm=T)
sp_info$ldmc_mean[which(sp_info$ldmc_mean==0)]<-NA
sp_info$LDMC<-NULL

sp_info<-merge(sp_info,new.sla,by="sp",all.x=T, all.y=F)
sp_info$sla_mean<-rowSums(sp_info[,c("sla_mean","SLA")],na.rm=T)
sp_info$sla_mean[which(sp_info$sla_mean==0)]<-NA
sp_info$SLA<-NULL

# Fill SLA gaps with reps:

head(sp_info)
head(e133_reps)

nosla2<-data.frame(sp=sp_info$sp[which(is.na(sp_info$sla_mean))])
nosla2<-tidy.df(nosla2)

sla.rep<-e133_reps[,c("sp","SLA")]

# There are only 4 species with missing SLA data that exist in the reps data:
nosla2<-data.frame(sp=nosla2$sp[which(nosla2$sp %in% levels(sla.rep$sp)==T)])
nosla2<-tidy.df(nosla2)
nosla2
sla.mg<-merge(nosla2,sla.rep,by="sp",all.x=F,all.y=F)
slamn<-data.frame(sp=levels(sla.mg$sp),SLA=tapply(sla.mg$SLA,sla.mg$sp,mean,na.rm=T))
slamn<-tidy.df(slamn)

sp_info<-merge(sp_info,slamn,by="sp",all.x=T, all.y=F)

sp_info$sla_mean<-rowSums(sp_info[,c("sla_mean","SLA")],na.rm=T)
sp_info$sla_mean[which(sp_info$sla_mean==0)]<-NA
sp_info$SLA<-NULL
} # close leaf data

#### ------ FILL TRAIT GAPS
{
# IMPORTANT: The congener means for the "sp" species and the invader trait gaps were calculated in R and added to the 'trait_gaps' sheet in excel which is used here to fill gaps. These would need to be updated if we get new trait data to full further gaps

# Tabulate data for filling gaps:
gxt<-xtabs(gaps$value~gaps$sp+gaps$trait, na.action=na.omit)
gap<-as.matrix(gxt[1:nrow(gxt),1:4])
gap[which(gap==0)]<-NA
gap<-as.data.frame.matrix(gap)
gap$sp<-rownames(gap)
rownames(gap)<-1:length(gap[,1])
gap<-gap[,c("sp","height","seed_mass","LDMC","SLA")]
head(gap)

# make sure this is TRUE TRUE:
dim(gap) == dim(sp_info[sp_info$sp %in% gap$sp,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")])

# make sure these are all TRUE:
gap$sp == sp_info[sp_info$sp %in% gap$sp,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]$sp

# This is the subset of sp_info that we're looking to fill with the data from gap (call it something so we can check to see that the right gaps got filled):
sp_info_orig<-sp_info[sp_info$sp %in% gap$sp,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]

# fill data with this interesting subset:
sp_info[sp_info$sp %in% gap$sp,c("hgt_mean","seed_mean","ldmc_mean","sla_mean")][is.na(gap[,2:length(gap)])==F]<-gap[,2:length(gap)][is.na(gap[,2:length(gap)])==F]

# then check that the correct data were filled:
sp_info[sp_info$sp %in% gap$sp,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]
gap

# and that the original data were not altered:
sp_info[sp_info$sp %in% gap$sp,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]
sp_info_orig

} # close trait gaps

#### ------ TRAIT METADATA
{
# Add trait metadata to species info:
trcols<-c("spdata","trdata","hgtdata","seeddata","ldmcdata","sladata","seeded")

# Make sure this is ordered by species:
sp_info<-sp_info[order(sp_info$species),]
is.unsorted(sp_info$species)

# Identify which trait data we have for each species:
sp_info$spdata<-levels(sp_info$species) %in% levels(sp_dat$species)
sp_info$trdata<-NA
sp_info$hgtdata<-1
sp_info$hgtdata[which(is.na(sp_info$hgt_mean))]<-0
sp_info$seeddata<-1
sp_info$seeddata[which(is.na(sp_info$seed_mean))]<-0
sp_info$ldmcdata<-1
sp_info$ldmcdata[which(is.na(sp_info$ldmc_mean))]<-0
sp_info$sladata<-1
sp_info$sladata[which(is.na(sp_info$sla_mean))]<-0
sp_info$seeded<-levels(sp_info$species) %in% levels(seeded$species)

# Change T/F to 1/0:
sp_info[which(sp_info[,trcols[1]]==T),trcols[1]]<-replace(which(sp_info[,trcols[1]]==T),TRUE,1)
sp_info[which(sp_info[,trcols[3]]==T),trcols[3]]<-replace(which(sp_info[,trcols[3]]==T),TRUE,1)
sp_info[which(sp_info[,trcols[4]]==T),trcols[4]]<-replace(which(sp_info[,trcols[4]]==T),TRUE,1)
sp_info[which(sp_info[,trcols[7]]==T),trcols[7]]<-replace(which(sp_info[,trcols[7]]==T),TRUE,1)

# Update the trait data column:
sp_info$trdata<-rowSums(sp_info[,trcols[3:5]],na.rm=T)
sp_info$trdata[sp_info$trdata>0]<-1
head(sp_info)

sp_info[sp_info$seeded==1,(c(1,8,length(sp_info)))]

# Check that we have full trait data for all seeded species (invaders):
head(sp_info)

invaders<-sp_info[sp_info$seeded==1,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]

} # close trait metadata

#### ------ PROCESS E275 DATA:
{
head(dsize)
head(dshape)
head(root)
head(rootMD)
head(plantDM)

dsize$dmass_new<-dsize$mass/dsize$no_dispersules
dsize$dmass_new<-round(dsize$dmass_new,6)
dsize$disp_mass<-round(dsize$disp_mass,6)

# check calculations (should all be true except for the three NAs):
dsize$disp_mass == dsize$dmass_new

# Calculate mean SRL for each species, for fine and very fine roots separately:
root$type<-as.factor(substr(root$sample,start=nchar(as.character(root$sample)),stop=nchar(as.character(root$sample))))

SRL.means<-aggregate(SRL~sp,data=root, FUN=mean)
SRL.vf<-aggregate(SRL~sp+type,data=root, FUN=mean)
SRL.sp<-data.frame(sp=SRL.means$sp)
SRL.sp<-data.frame(sp=SRL.means$sp,SRL_all=SRL.means$SRL,SRL_fine=merge(SRL.sp,SRL.vf[SRL.vf$type=="f",],by="sp",all.x=T)$SRL,SRL_veryfine=merge(SRL.sp,SRL.vf[SRL.vf$type=="v",],by="sp",all.x=T)$SRL)

range(SRL.sp$SRL_all,na.rm=T)
range(SRL.sp$SRL_fine,na.rm=T)
range(SRL.sp$SRL_veryfine,na.rm=T)


} # close E275 data

#### ------ TIDY SUBPLOT & SPECIES DATA:
{
# Remove subplots I, J and K from species abundance data (not sure what these are for and there is only one record each):
sp_dat<-sp_dat[,1:(length(sp_dat)-1)]
sp_dat<-sp_dat[-which(sp_dat$subplot=="I"),]
sp_dat<-sp_dat[-which(sp_dat$subplot=="J"),]
sp_dat<-sp_dat[-which(sp_dat$subplot=="K"),]
sp_dat$plsub<-as.factor(paste(sp_dat$plot, sp_dat$subplot, sep=""))
sp_dat<-tidy.df(sp_dat)
head(sp_dat)

seeded<-seeded[,1:(length(seeded)-1)]
seeded$plsub<-as.factor(paste(seeded$plot, seeded$subplot, sep=""))

# Separate seeded and unseeded species:
seeded.species<-as.character(levels(seeded$sp))
sd.sp<-data.frame(sp=seeded.species)
unseeded.species<-as.character(sp_info$sp[which(sp_info$sp %in% sd.sp$sp ==F)])
unsd.sp<-data.frame(sp=unseeded.species)

# Check the lengths match:
if(length(unsd.sp[,1]) + length(sd.sp[,1])!=length(levels(sp_info$sp))) "something is wrong"

# Make data frame of subplots. I'm taking this from the species data, so if nothing was ever recorded at a plot in a given year, it won't be included in the analysis:
subpls<-sp_dat[,c("year","plot","subplot")]
subpls$plsub<-as.factor(paste(subpls$plot,subpls$subplot,sep=""))
subpls<-unique(subpls)

# Add columns for subsetting on seeded and satellite:
subpls$seeded<-subpls$plsub %in% seeded$plsub
subpls$seeded[which(subpls$seeded==T)]<-1
subpls$satellite<-as.character(subpls$subplot)<LETTERS[9]
subpls$satellite[which(subpls$satellite==T)]<-1
subpls<-tidy.df(subpls)
head(subpls)

###-- DEFINE SUBPLOTS:

# Seeded subplots, colonisation phase:
sd.sub<-unique(as.character(subpls$plsub[subpls$seeded==1]))

# Unseeded subplots, naturalisation phase:
unsd.sub<-unique(as.character(subpls$plsub[subpls$seeded==0 & subpls$satellite==0]))

# Remove the W and Y subplots as we will only use two unseeded subplots per plot:
unsd.sub<-unsd.sub[-c(grep("W",unsd.sub),grep("Y",unsd.sub))]

###-- REMOVE TWO SEEDED SPECIES NEVER RECORDED:

# Remove these from:
# sd.sp
# seeded

# seeded species
head(sd.sp)
head(seeded)

# species data
head(sp_dat)

for.removal<-levels(seeded$sp)[which(!levels(seeded$sp) %in% levels(sp_dat$sp))]

# REMOVE two species never recorded (Ast_nov and Mim_rin) from sd.sp:
sd.sp<-data.frame(sp=sd.sp[-which(sd.sp$sp %in% for.removal),])
sd.sp<-tidy.df(sd.sp)

# REMOVE two species never recorded (Ast_nov and Mim_rin) from seeded:
seeded<-seeded[-which(seeded$sp %in% for.removal),]
seeded<-tidy.df(seeded)

# remove species name:
seeded<-seeded[,-which(colnames(seeded) %in% c("species"))]

## ~~ Determine no. seeded (subplot level). 

# This should range from zero (some "seeded subplots" did not have any species seeded) to 52 (54 seeded species minus the two species never recorded). There should be values only for 54 seeded subplots because those "seeded subplots" that didn't receive any species are not part of the analysis. Thus, there is nothing more to do here than to summarise np1 as follows: 

np1<-data.frame(plsub=names(tapply(seeded$sp,seeded$plsub,function(x)length(unique(x)))),no_seeded=tapply(seeded$sp,seeded$plsub,function(x)length(unique(x))))
np1$no_seeded<-as.numeric(np1$no_seeded)
np1<-tidy.df(np1)
head(np1)
dim(np1)
range(np1$no_seeded)
range(no_pl$no_planted)

# Note that NUMBER PLANTED (E93_no_planted.txt) and SEEDED (E93_seeded_species.txt) do not line up. For example 12W has 53 species seeded in the original file, but the number planted file indicates 54 species. This is not a problem as we do not use the number planted file at all:
seeded[seeded$plsub=="12W",]

} # close subplots

#### ------ TRAIT EXTREMES:
{
# Get the upper and lower values for each trait for use in the gap calcs:
head(hgt_dat,20)
head(seed_dat)
head(lf_dat)

# make a vector of sp for separating woodies:
woodies<-as.character(sp_info$sp[sp_info$func_grp=="W"])

seednowood<-seed_dat[-which(seed_dat$sp %in% woodies),]
seednowood<-tidy.df(seednowood)

trait_extr<-data.frame(trait=c("hgt","seed","ldmc","sla"),rbind(range(hgt_dat$height,na.rm=T),range(seednowood$mps,na.rm=T),range(lf_dat$LDMC,na.rm=T),range(lf_dat$SLA,na.rm=T)))
colnames(trait_extr)[2:3]<-c("smallest","largest")
trait_extr

} # close trait extremes

} # close process raw

######################################
#    PROBABILITY OF COLONISATION	 #
######################################
{

#### ------ BASE COLONISE DATA:
{
# Calculate probability of colonisation for all seeded species, removing species that were anywhere in the plot in 1991:

# seeded subplots
sd.sub

# seeded species
sd.sp
head(seeded)

# species data
head(sp_dat)

# survey years
syr<-unique(sp_dat$year)

# use 'seeded' as the base data. It contains every seeded subplot and the species that were seeded. This can be replicated for each survey year and a presence/absence column added from the species abundance data.

# Remove satellite plots
abund<-sp_dat[sp_dat$subplot %in% c("W","X","Y","Z"),]
abund<-tidy.df(abund)

# remove species name:
abund<-abund[,-which(colnames(abund) %in% c("species"))]

####--- 
# Remove species x subplot combinations where that species was anywhere on the plot in 1991:
dat91<-abund[abund$year==1991,]
dat91<-dat91[,-which(colnames(dat91) %in% c("species","value"))]
dat91<-tidy.df(dat91)
s.rm<-seeded

head(s.rm)
# create species x subplot variable:
s.rm$sp_plot<-as.factor(paste(s.rm$sp,s.rm$plot,sep="_"))
dat91$sp_plot<-as.factor(paste(dat91$sp,dat91$plot,sep="_"))

# simplify data for merging:
s.rm<-s.rm[,c("sp","plsub","sp_plot")]
dat91<-dat91[,c("std_cover","sp_plot")]

# we only need one record of each species being on a plot in 1991, so remove dups:
dat91<-dat91[-which(duplicated(dat91$sp_plot)),]
dat91<-tidy.df(dat91)

# Merge to summarise presence in 1991:
s.add<-merge(s.rm,dat91,by=c("sp_plot"),all.x=T,all.y=F)
NAbefore<-length(which(is.na(s.add$std_cover)))
s.add$std_cover[which(is.na(s.add$std_cover))]<-0
# This should be T (just making sure there were no zeros in the std_cover variable):
length(which(s.add$std_cover==0)) == NAbefore
s.add$std_cover[which(s.add$std_cover!=0)]<-1
# This should be T:
(length(which(s.add$std_cover==1)) + NAbefore) == dim(s.add)[1]
s.add<-s.add[,c("sp","plsub","std_cover")]
colnames(s.add)[colnames(s.add) %in% "std_cover"]<-"was_it_there1991"

# CHECK it:
sp.now<-as.character(sample(s.add$sp,1))
plsub.now<-as.character(sample(s.add$plsub,1))
s.add[s.add$sp==sp.now & s.add$plsub==plsub.now,]
# If the result from the prev line == 0, this should come up with nothing, or with a list of records that doesn't include 1991. If the result is 1, then you should see 1991 in the list:
abund[abund$sp==sp.now & abund$plot==substr(plsub.now,1,unlist(gregexpr("[A-Z]",plsub.now))-1),]

# remove those that were there in 1991:
s<-s.add[-which(s.add$was_it_there1991==1),]
s$was_it_there1991<-NULL
s$plot<-as.numeric(substr(s$plsub,1,unlist(gregexpr("[A-Z]",s$plsub))-1))
s<-s[order(s$plot,s$plsub,s$sp),]
s<-tidy.df(s)
s<-s[,-which(colnames(s) %in% c("plot"))]

# remove plot and subplot:
abund<-abund[,-which(colnames(abund) %in% c("plot","subplot","value"))]
head(s)
head(abund)

# list for the output:
out.list<-list()

for (i in 2:length(syr)){

# do each year separately:
data.thisrun<-abund[abund$year==syr[i],]

# remove rows that are not seeded subplots:
data.thisrun<-data.thisrun[data.thisrun$plsub %in% sd.sub,]

head(data.thisrun)
head(s)

# remove rows that are not seeded species:
data.thisrun<-data.thisrun[data.thisrun$sp %in% s$sp,]
data.thisrun<-tidy.df(data.thisrun)

# add std_cover to base data:
mg.dat<-merge(s,data.thisrun,by=c("sp","plsub"),all.x=T,all.y=F)

# update year
mg.dat$year<-syr[i]

# turn std_cover into presence absence:
mg.dat$std_cover[which(!is.na(mg.dat$std_cover))]<-1
mg.dat$std_cover[which(is.na(mg.dat$std_cover))]<-0

# add to list:
out.list[[i]]<-mg.dat

} # close i survey year

colonise<-data.frame(do.call(rbind, out.list))

# random checks:
sp.now<-sample(sd.sp$sp,1)
sub.now<-sample(sd.sub,1)
yr.now<-sample(syr,1)
colonise[colonise$plsub==sub.now & colonise$year==yr.now,]
sp_dat[sp_dat$plsub==sub.now & sp_dat$year==yr.now,]

# update colname:
colnames(colonise)[which(colnames(colonise)=="std_cover")]<-"pr"
colonise$prev_present<-NULL

# add plot and subplot back in (removed in loop for merging):
colonise$plot<-substr(colonise$plsub,1,unlist(gregexpr("[A-Z]",colonise$plsub))-1)
colonise$subplot<-substr(colonise$plsub,unlist(gregexpr("[A-Z]",colonise$plsub)),unlist(gregexpr("[A-Z]",colonise$plsub)))

# Reduce species info to species with trait data and relevant cols: 
spinf<-sp_info[sp_info$trdata==1,]
spinf<-tidy.df(spinf)
spinf<-spinf[,c("sp","hgt_mean","seed_mean","ldmc_mean","sla_mean")]
head(spinf)
head(colonise)

} # close base colonise data

#### ------ FUNCTIONAL DISTANCE:
{
# Calculate continuous trait variables, for each seeded species (s), in each year, in each subplot:

# 1. RAW = actual trait value

# 2. GAP = gap size; nearest neighbour high - nearest neighbour low

# 3. relFD = relative functional difference; trait(s) - trait(community weighted mean):

# 4. relNN = relative neighbour; trait(s) - trait(nearest neighbour):

# Before we begin, work out 'raw' factorial variables. Note we need to do this for all species in the community so we can work out relative values. The spinf data has been subsetted on those species with continuous trait data, so don't use that either, go back to sp_info

head(sp_info)

ftraits<-sp_info[,c("sp","duration","func_grp","grass","legume")]
colnames(ftraits)[colnames(ftraits)  %in% "duration"]<-"lifespan"

# Re-group trait categories
ftraits$lifespan[-which(ftraits$lifespan=="A")]<-"P"

# From JC: I think it'll be OK to just look at four community chatracteristics for this (e.g. %C3 etc) but look at five groups as the invaders. it'll just mean that we'll see how woody success varies with %C3, % c4 etc; we won't know how it'll vary with %woody. [the plots were meant to just be herbaceous. Some woody spp have since grown, but there are few woody spp present overall, so not a big deal. I also suspect that we'll just ignore the woody invaders too, though probably bes to look at initially at least).
# So keep functional group as is.
ftraits<-tidy.df(ftraits)

# Add trait data to colonise data:
ft<-ftraits
head(ft)

colonise<-merge(colonise, ft, by="sp", all.x=T, all.y=F)
head(colonise)

# check:
sp.now<-sample(levels(colonise$sp),1)
ftraits[ftraits$sp==sp.now,]
head(colonise[colonise$sp==sp.now, ])

# WARNING: this is a SLOW (about 2 min) row-by-row loop at each row of the colonise data. 

# New variables:

# The prop_native and prop_intro variables both relate to the proportion native, the difference is where the unknown species were grouped (into native in prop_native and introduced in prop_intro):
colonise$prop_peren<-NA
colonise$prop_grass<-NA
colonise$prop_legume<-NA

colonise$hgt_cwm<-NA
colonise$seed_cwm<-NA
colonise$ldmc_cwm<-NA
colonise$sla_cwm<-NA

colonise$hgt_raw<-NA
colonise$seed_raw<-NA
colonise$ldmc_raw<-NA
colonise$sla_raw<-NA

colonise$hgt_gap<-NA
colonise$seed_gap<-NA
colonise$ldmc_gap<-NA
colonise$sla_gap<-NA

colonise$hgt_relFD<-NA
colonise$seed_relFD<-NA
colonise$ldmc_relFD<-NA
colonise$sla_relFD<-NA

colonise$hgt_relNN<-NA
colonise$seed_relNN<-NA
colonise$ldmc_relNN<-NA
colonise$sla_relNN<-NA

traits<-c("hgt","seed","ldmc","sla")

for (i in 1:length(colonise[,1])){

# define species, subplot and year for each row:
sp<-as.character(colonise$sp[i])
plsub<-as.character(colonise$plsub[i])
year<-as.character(colonise$year[i])

# get community abundance data for the current subplot and year:
comm.data<-abund[abund$plsub==plsub & abund$year==year,]

# remove the target species from the community data:
if(sp %in% comm.data$sp==T) comm.data<-comm.data[-which(comm.data$sp==sp),]

# add continuous trait data to the community data:
comm.data<-merge(comm.data,spinf,by="sp",all.x=T,all.y=F)

# the spinf data used in the prev line has been subsetted on species that have continuous traits. For categorical trait data for the community, use the ft data set
comm.data<-merge(comm.data,ft,by="sp",all.x=T,all.y=F)

# get continuous trait data for the current species:
tr.data<-spinf[spinf$sp==sp,]

# add categorical trait data for the current species:
tr.data<-merge(tr.data,ft,by="sp",all.x=T,all.y=F)

colonise[(i-1):(i+1),]
comm.data

# Calculate values of categorical traits relative to the community:
comm.cover<-sum(comm.data$std_cover)

colonise$prop_peren[i]<-sum(comm.data$std_cover[comm.data$lifespan=="P"])/comm.cover
colonise$prop_grass[i]<-sum(comm.data$std_cover[comm.data$grass==1])/comm.cover
colonise$prop_legume[i]<-sum(comm.data$std_cover[comm.data$legume==1])/comm.cover

no.incomm<-length(comm.data[,1])

colonise[(i-1):(i+1),]

comm.data
tr.data

# calculate FD measures only when there is trait data for the current species:

if(length(tr.data[,1])>0){

# COMMUNITY WEIGHTED MEAN:
colonise$hgt_cwm[i]<-weighted.mean(comm.data$hgt_mean,comm.data$std_cover, na.rm=T)
colonise$ldmc_cwm[i]<-weighted.mean(comm.data$ldmc_mean,comm.data$std_cover, na.rm=T)
colonise$sla_cwm[i]<-weighted.mean(comm.data$sla_mean,comm.data$std_cover, na.rm=T)

# RAW:
colonise$hgt_raw[i]<-tr.data$hgt_mean
colonise$seed_raw[i]<-tr.data$seed_mean
colonise$ldmc_raw[i]<-tr.data$ldmc_mean
colonise$sla_raw[i]<-tr.data$sla_mean

# relative FD:
colonise$hgt_relFD[i]<-tr.data$hgt_mean-weighted.mean(comm.data$hgt_mean,comm.data$std_cover, na.rm=T)
colonise$ldmc_relFD[i]<-tr.data$ldmc_mean-weighted.mean(comm.data$ldmc_mean,comm.data$std_cover, na.rm=T)
colonise$sla_relFD[i]<-tr.data$sla_mean-weighted.mean(comm.data$sla_mean,comm.data$std_cover, na.rm=T)

# Use "seed.comm" for ALL seed calculations (same as comm.data but with woodies removed where present):

# If there are woodies present, remove them and tidy seed.comm:
if(length(which(comm.data$sp %in% woodies))>0){
seed.comm<-comm.data
seed.comm<-seed.comm[-which(seed.comm$sp %in% woodies),]
seed.comm<-tidy.df(seed.comm)
colonise$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
colonise$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
}
# If there are no woodies, make seed.comm == comm.data. We will use seed.comm in the NN and gap calcs, so define it here regardless of whether there are woodies or not:
if(length(which(comm.data$sp %in% woodies))==0){
seed.comm<-comm.data
colonise$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
colonise$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
}

# relative NN and GAP, one trait at a time (j):

colonise[(i-1):(i+1),]

# store comm.data in a separate data frame for updating in the j trait loop:
comm.data.store<-comm.data

for (j in 1:length(traits)){

trait.thisrun<-traits[j]

# update comm.data, depending on the trait:

if (trait.thisrun=="seed") comm.data<-seed.comm
if (trait.thisrun!="seed") comm.data<-comm.data.store

# get trait data for the current community and target species and order them by trait value:
nn<-data.frame(tr=comm.data[,which(colnames(comm.data)==paste(trait.thisrun,"mean",sep="_"))],type="comm",stringsAsFactors=F)
nn<-rbind(nn,c(tr.data[,which(colnames(tr.data)==paste(trait.thisrun,"mean",sep="_"))],"target"))
nn$tr<-as.numeric(nn$tr)
nn<-nn[order(nn$tr),]
nn<-tidy.df(nn)
nn

# if an NA belongs to the target species, skip to the next trait j:
if(is.na(nn$tr[nn$type=="target"])) next

# if an NA belongs to someone else in the community, remove NAs:
if(length(which(is.na(nn$tr)))>0) nn<-nn[-which(is.na(nn$tr)),]
nn

# identify the row occupied by the target:
target.ind<-which(nn$type=="target")

# **** SINGLE NEIGHBOUR (smallest & largest):

# If the species is the SMALLEST, take the next smallest as the NN and the full range as the gap:
if(target.ind==1) {

## *** NEAREST NEIGHBOUR:	
# we are not using NN, but I've kept the code in case we need it later:
colonise[i,which(colnames(colonise)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# nn_distance is the distance between the invader who happens to be the smallest for the current trait in the current community and its nearest, next biggest neighbour:
nn_distance<-nn$tr[target.ind+1]-nn$tr[target.ind]

# invader_mean is the mean value for the invader for the current trait:
invader_mean<-tr.data[,grep(trait.thisrun,colnames(tr.data))]

# If NN distance > invader mean trait value, then gap size = NN trait value, otherwise gap = NN distance*2

if (nn_distance > invader_mean) colonise[i,which(colnames(colonise)==paste(trait.thisrun,"gap",sep="_"))]<-nn$tr[target.ind+1] else colonise[i,which(colnames(colonise)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind+1]-nn$tr[target.ind])*2

} # close if smallest

colonise[(i-1):(i+1),]

# If the species is the LARGEST, take the next largest as the NN and the full range as the gap:
if(target.ind==length(nn[,1])) {

## *** NEAREST NEIGHBOUR:	
colonise[i,which(colnames(colonise)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# gap = NN distance*2

colonise[i,which(colnames(colonise)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind]-nn$tr[target.ind-1])*2

} # close if largest

# **** DOUBLE NEIGHBOUR:

colonise[(i-1):(i+1),]

# calculate the distance between the neighbour on the small side and the neighbour on the big side:
if(target.ind>1 & target.ind<length(nn[,1])){
lowNN.dist<-nn$tr[target.ind]-nn$tr[target.ind-1]
highNN.dist<-nn$tr[target.ind+1]-nn$tr[target.ind]

# we want relative NN distance, so each calculation should be the nearest neighbour minus the target species . The code must preserve the sign, rather than just take the smallest distance:

# 21 March 2016 update: switched the sides of the calculation so that, for example, the taller nearest neighbour would be positive, the shorter nearest neighbour would be negative. This should be easier to interpret. Corresponding changes have been made for single (largest and smallest) neighbour above.

if(lowNN.dist < highNN.dist) colonise[i,which(colnames(colonise)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

if(lowNN.dist > highNN.dist) colonise[i,which(colnames(colonise)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

# When low and high neighbour have the same distance, we were getting NAs for relNN values. The next line deals with those cases, and we're using a positive value, from JC email (21 March 2016): Let's use a positive value when the upper and lower NN are the same.  
if(lowNN.dist == highNN.dist) colonise[i,which(colnames(colonise)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

colonise[(i-1):(i+1),]

# gap distance:
colonise[i,which(colnames(colonise)==paste(trait.thisrun,"gap",sep="_"))]<-lowNN.dist + highNN.dist

} # close if species has two neighbours

} # close j traits

} # close if species has trait data

} # close i rows of data frame

# Add absolute FD:
colonise$hgt_absFD<-abs(colonise$hgt_relFD)
colonise$seed_absFD<-abs(colonise$seed_relFD)
colonise$ldmc_absFD<-abs(colonise$ldmc_relFD)
colonise$sla_absFD<-abs(colonise$sla_relFD)

# check:
colonise[sample(1:length(colonise[,1]),4),]
# head(colonise[which(colonise$sp=="Nep_cat"),])

} # close functional distance

#### ------ ADD ENVIRO VARIABLES:
{
## ~~ ADD: no. planted (subplot level). np1 contains the number of species seeded for each seeded subplot, excluding those two species with unviable seed. It was defined in the "Tidy subplot & species data" section of "Data setup":
head(np1)
head(colonise,3)

colonise<-merge(colonise,np1,by="plsub",all.x=T, all.y=F)

# minus one for the target invader:
colonise$no_seeded<-colonise$no_seeded-1
head(colonise,3)

## ~~ ADD: elevation (plot level):
tp<-topog[,c("plot","masl")]
colonise<-merge(colonise,tp,by="plot",all.x=T, all.y=F)

## ~~ ADD: NO3 and MOISTURE (subplot level):
no3<-soil.mean[,c("plsub","NO3soil","moisture")]
colonise<-merge(colonise, no3, by="plsub", all.x=T, all.y=F)

## ~~ ADD: BURN and PRECIP (YEAR level):
head(burn_precip)
colonise<-merge(colonise, burn_precip, by="year", all.x=T, all.y=F)
head(colonise,4)

## ~~ ADD: SOIL DISTURBANCE (subplot x year):
# From JC: Lets just use  bare ground by itself, and all three of the other temporal vars ("dsince_fire","rfall_mm","year").  

# Tidy disturbance data:
soil_dist$plsub<-paste(soil_dist$plot,soil_dist$subplot,sep="")
plsubs<-c("W","X","Y","Z")
soil_dist<-soil_dist[soil_dist$subplot %in% plsubs,]
soil_dist<-tidy.df(soil_dist)
head(soil_dist)

# Summarise bare ground disturbance:
soil.bg<-soil_dist[soil_dist$disturb=="bare",]
soil.bg<-tidy.df(soil.bg)
head(soil.bg)
dist1<-aggregate(std_cover~plsub+year,sum,data=soil.bg)
head(dist1)

# Check
yr.now<-sample(soil_dist$year,1)
plsub.now<-sample(soil_dist$plsub,1)
soil_dist[soil_dist$year==yr.now & soil_dist$plsub==plsub.now,]
dist1[dist1$year==yr.now & dist1$plsub==plsub.now,]

# Change the colname to soil_dist:
colnames(dist1)[which(colnames(dist1)=="std_cover")]<-"soil_dist"

# Merge with main data:
colonise<-merge(colonise,dist1,by=c("year","plsub"),all.x=T,all.y=F)
colonise[sample(1:length(colonise[,1]),size=4),]

# Make soil disturbance NAs zero (nothing was recorded):
colonise$soil_dist[which(is.na(colonise$soil_dist))]<-0

} # close enviro variables

#### ------ TIDY AND SCALE:
{
## ~~ Order by plot, subplot, year:
colonise$plot<-as.numeric(colonise$plot)
colonise<-colonise[order(colonise$plot, colonise$subplot, colonise$year),]
colonise<-tidy.df(colonise)

## ~~ Factorise:
colonise$plot<-as.factor(colonise$plot)
colonise$subplot<-as.factor(colonise$subplot)

# Make year numeric, beginning at zero:
colonise$yr<-colonise$year-min(colonise$year)

## ~~ Scale ALL numeric variables (exclude response, raw year, numeric year, factors and functional dummies):
to_scale<-sapply(colonise,is.numeric)
to_scale [names(to_scale) %in% c("pr","year","yr","grass","legume")]<-FALSE
to_scale

colScale<-colonise
colScale[,to_scale]<-scale(colScale[,to_scale])
sapply(colScale[,to_scale],range,na.rm=T)
head(colScale,4)

# Remove 'woody' invaders
coloniseSc<-colScale[-which(colScale$func_grp=="W"),]
coloniseSc<-tidy.df(coloniseSc)
check.rows(coloniseSc,2)

# Make a new data frame from the colonise (unscaled) base data, with woodies removed and functional groups added as dummies. This will be directly comparable to coloniseSc, and can be used to plot the x-axis using unscaled variables:
head(colonise,2)

col.unscaled<-colonise

# Remove 'woody' invaders
col.unscaled<-col.unscaled[-which(col.unscaled$func_grp=="W"),]
col.unscaled<-tidy.df(col.unscaled)
check.rows(col.unscaled,2)

# Double check that none of the data sets have NAs:
head(colonise,3)
head(colScale,3)
head(coloniseSc,3)

na.colonise<-data.frame(pred=names(apply(colonise,2,function(x) length(which(is.na(x))))),no_NAs=apply(colonise,2,function(x) length(which(is.na(x)))))
na.colonise<-tidy.df(na.colonise)

na.colScale<-data.frame(pred=names(apply(colScale,2,function(x) length(which(is.na(x))))),no_NAs=apply(colScale,2,function(x) length(which(is.na(x)))))
na.colScale<-tidy.df(na.colScale)

na.coloniseSc<-data.frame(pred=names(apply(coloniseSc,2,function(x) length(which(is.na(x))))),no_NAs=apply(coloniseSc,2,function(x) length(which(is.na(x)))))
na.coloniseSc<-tidy.df(na.coloniseSc)

# These should all be zero:
length(which(na.colonise$no_NAs>0))
length(which(na.colScale$no_NAs>0))
length(which(na.coloniseSc$no_NAs>0))

# Re-level lifespan to P in ALL data sets
colonise$lifespan<-relevel(colonise$lifespan,ref="P")
colScale$lifespan<-relevel(colScale$lifespan,ref="P")
coloniseSc$lifespan<-relevel(coloniseSc$lifespan,ref="P")
col.unscaled$lifespan<-relevel(col.unscaled$lifespan,ref="P")

} # close tidy and scale

} # close probability of colonisation

######################################
# 	   PROBABILITY OF SPREAD		 #
######################################
{

#### ------ BASE SPREAD DATA:
{
# Calculate probability of spread for all seeded species, removing species that were anywhere in the plot in 1991. 

# unseeded subplots (note there are 66 subplots, some W and Y plots are included that were not seeded, there were only 54 seeded subplots in total out of 120)
unsd.sub

# seeded subplots
sd.sub

# survey years
syr<-unique(sp_dat$year)

# Remove 1998 and 2004 as we have missing unseeded subplot data from those years:
spr.syr<-syr[-which(syr %in% c(1998,2004))]

# The SPREAD phase will be analysed at the PLOT level.  

# The data set 's' contains all seeded species, with those that were present in 1991 removed (see 'probability of colonisation' code for derivation). This will be the 'base' data for each year, i.e., in each year of interest (1993 onwards), have each of the seeded species that were not present in 1991 anywhere on the plot, spread to the unseeded subplots? It's a probability, so will take into account failed spreads (0) as well as successes (1)
head(s,10)

# If the data are correct, Ast_nov and Mim_rin should not appear and there should be only 52 seeded species:
which(s$sp=="Ast_nov"); which(s$sp=="Mim_rin"); length(levels(s$sp))

# Put PLOT back into 's' and remove plsub:
s.plot<-s
s.plot$plot<-substr(s.plot$plsub,1,unlist(gregexpr("[A-Z]",s.plot$plsub))-1)
s.plot[sample(1:length(s.plot[,1]),10),]
s.plot$plsub<-NULL

# Remove duplicate species x plot combinations:
s.plot$sppl<-paste(s.plot$sp,s.plot$plot,sep="")
s.plot<-s.plot[-which(duplicated(s.plot$sppl)),]
s.plot$sppl<-NULL
s.plot<-tidy.df(s.plot)
head(s.plot)

# The data set 'abund' contains species abundances, with satellite plots removed and the data columns simplified for summarising. This was used for colonisation, but had its raw abundance value (value) chopped. Make a new data frame in the same way but without removing the value column:
head(abund)

# Remove satellite plots
spr.abund<-sp_dat[sp_dat$subplot %in% c("W","X","Y","Z"),]
spr.abund<-spr.abund[,-which(colnames(spr.abund) %in% c("species"))]
spr.abund<-tidy.df(spr.abund)

# remove rows that are not seeded species from abundance data:
spr.abund<-spr.abund[spr.abund$sp %in% s.plot$sp,]
spr.abund<-tidy.df(spr.abund)
if(length(which(!spr.abund$sp %in% s.plot$sp))>0) print("data subset is wrong")

# in the loop we need abundances for the SEEDED subplots and presences/absence for the UNSEEDED subplots

unsd.sub # unseeded subplots
sd.sub # seeded subplots

# So make a SEEDED subplot abundance data frame and one UNSEEDED subplot abundance data frame
head(spr.abund)

sd.abund<-spr.abund[spr.abund$plsub %in% sd.sub,]
unsd.abund<-spr.abund[spr.abund$plsub %in% unsd.sub,]
sd.abund<-tidy.df(sd.abund)
unsd.abund<-tidy.df(unsd.abund)

# These contain the abundances for SEEDED SPECIES on seeded and unseeded subplots:
head(sd.abund)
head(unsd.abund)

# Add a "sp_plot" variable for merging to the three data frames of interest. Species by plot is the observation level.
sd.abund$sp_plot<-paste(sd.abund$sp,sd.abund$plot,sep="_")
unsd.abund$sp_plot<-paste(unsd.abund$sp,unsd.abund$plot,sep="_")
s.plot$sp_plot<-paste(s.plot$sp,s.plot$plot,sep="_")

# list for the output:
outlist.spr<-list()

# This will begin in 1993. No seeded species on the plot in 1991 are included.

# The loop is done one year at a time, for all seeded species in year i
# It replicates the base spread data (s.plot), and gives a probability of spread for each year
# The result should be dim(s.plot)[1] * length(3:length(spr.syr)), that is the full base spread data, times the number of years of interest (five years)
head(s.plot)

for (i in 3:length(spr.syr)){

# SPREAD rules depending on sampling time/frequency:

# EARLY samples: 1993 and 1994, the species MUST occur on SEEDED SUBPLOTS at sometime before the current year, not including the current year: 
##--- For 1993, the species must have occurred in 1992
##--- For 1994, the species must have occurred in 1992 and/or 1993

# LATE samples: 1997, 2004, 2008 and 2012, the species MUST occur somewhere on SEEDED SUBPLOTS in the year 1992 onwards, including the current year

# Determine if the species was present on SEEDED SUBPLOTS in the previous years, following the above rules:

this.year<-spr.syr[i]

if (this.year < 1995) years.before<-spr.syr[which(spr.syr=="1992"):(which(spr.syr==this.year-1))]

if (this.year > 1995) years.before<-spr.syr[which(spr.syr=="1992"):(which(spr.syr==this.year))]

splot.thisrun<-s.plot

# For each species x plot combination in splot.thisrun, determine whether it was present in the "years.before" on SEEDED subplots:
seeded.yearsbefore<-sd.abund[sd.abund$year %in% years.before,]
seeded.yearsbefore<-tidy.df(seeded.yearsbefore)
head(seeded.yearsbefore)
head(splot.thisrun)

splot.thisrun$prev_pres<-splot.thisrun$sp_plot %in% seeded.yearsbefore$sp_plot 

# For each species x plot combination in splot.thisrun, determine if it has spread to the UNSEEDED subplots in the current year:
unsd.thisyear<-unsd.abund[unsd.abund$year %in% this.year,]
unsd.thisyear<-tidy.df(unsd.thisyear)
splot.thisrun$pres_thisyr<-splot.thisrun$sp_plot %in% unsd.thisyear$sp_plot 
head(unsd.thisyear)
head(splot.thisrun)

splot.thisrun$spread<-NA

# Code spread only for species that received T for prev_pres:
splot.thisrun$spread[splot.thisrun$prev_pres==T & splot.thisrun$pres_thisyr==F]<-0
splot.thisrun$spread[splot.thisrun$prev_pres==T & splot.thisrun$pres_thisyr==T]<-1

check.rows(splot.thisrun)
head(splot.thisrun)

# Remove rows for species that received F for prev_pres:
splot.thisrun<-splot.thisrun[-which(is.na(splot.thisrun$spread)),]
splot.thisrun<-tidy.df(splot.thisrun)
splot.thisrun$prev_pres<-NULL
splot.thisrun$pres_thisyr<-NULL

head(splot.thisrun)
check.rows(splot.thisrun)

# Code spread RAW ABUNDANCE only for species that received 1 for spread:

# First, summarise abundance across the two subplots:
spab.thisyear<-unsd.thisyear[,c("sp_plot","value")]
spab2<-aggregate(value~sp_plot,sum,data=spab.thisyear)

# then add to df:
splot.thisrun<-merge(splot.thisrun,spab2,by="sp_plot",all.x=T,all.y=F)

# change the name of value so it doesn't conflict with cumulative abundance below:

colnames(splot.thisrun)[colnames(splot.thisrun)=="value"]<-"spr_abund_raw"

# There are no positive values for abundance on any of the spread==0 rows because they were removed in the previous step. This should be zero:
if (length(which(is.na(splot.thisrun$value[splot.thisrun$spread==0])==F))!=0) stop(print("something is wrong"))

# There are positive values on all of the spread==1 rows. Check this is zero:
if (length(which(is.na(splot.thisrun$value[splot.thisrun$spread==1])))!=0) stop(print("something is wrong"))

# Code spread STANDARDISED ABUNDANCE only for species that received 1 for spread:

# First, summarise abundance across the two subplots:
stdab.thisyear<-unsd.thisyear[,c("sp_plot","std_cover")]
spab3<-aggregate(std_cover~sp_plot,sum,data=stdab.thisyear)

# then add to df:
splot.thisrun<-merge(splot.thisrun,spab3,by="sp_plot",all.x=T,all.y=F)

# There are no positive values for abundance on any of the spread==0 rows because they were removed in the previous step. This should be zero:
if (length(which(is.na(splot.thisrun$std_cover[splot.thisrun$spread==0])==F))!=0) stop(print("something is wrong"))

# There are positive values on all of the spread==1 rows. Check this is zero:
if (length(which(is.na(splot.thisrun$std_cover[splot.thisrun$spread==1])))!=0) stop(print("something is wrong"))

# Create a data frame for calculating the cumulative abundance in seeded subplots for years 1992-year i, inclusive of year i. This will be identical to the "seeded.yearsbefore" data above for survey years beyond 1995, but different for early surveys before 1995. So create a new df to use on all surveys 
cumul.yrs<-spr.syr[which(spr.syr=="1992"):(which(spr.syr==this.year))]
cumul.abund<-sd.abund[sd.abund$year %in% cumul.yrs,]
cumul.abund<-tidy.df(cumul.abund)

head(cumul.abund)
head(splot.thisrun)

# Sum abundances over all years, including year i:
abund.agg<-aggregate(value~sp_plot,data=cumul.abund, FUN=sum)

# Include only sp_plot combinations present in the data set for this year:
abund.agg<-abund.agg[abund.agg$sp_plot %in% splot.thisrun$sp_plot,]
abund.agg<-tidy.df(abund.agg)
head(abund.agg)

# check:
spp<-sample(abund.agg$sp_plot,1)
abund.agg[which(abund.agg$sp_plot==spp),]
cumul.abund[which(cumul.abund$sp_plot==spp),]

# Merge abundance into the main data frame:
splot.thisrun<-merge(splot.thisrun,abund.agg,by="sp_plot",all.x=T, all.y=F)
head(splot.thisrun)
colnames(splot.thisrun)[which(colnames(splot.thisrun)=="value")]<-"sdsub_abund"

# Add a year column and tidy up:
splot.thisrun$year<-this.year
splot.thisrun<-splot.thisrun[,c("year","sp","plot","spread","spr_abund_raw","std_cover","sdsub_abund")]
splot.thisrun<-tidy.df(splot.thisrun)

# add to list:
outlist.spr[[i]]<-splot.thisrun

} # close i survey year

spread<-data.frame(do.call(rbind, outlist.spr))
head(spread,3)
head(spr.abund,3)

# random checks (do many of these):
# X and Z are the unseeded
row.now<-spread[sample(1:length(spread[,1]),1),]
sp.now<-as.character(row.now$sp)
pl.now<-row.now$plot
yr.now<-row.now$year
spread[spread$plot==pl.now & spread$year==yr.now & spread$sp==sp.now,]
spr.abund[spr.abund$plot==pl.now & spr.abund$sp==sp.now,]

} # close base spread

#### ------ FUNCTIONAL DISTANCE:
{
# spinf is the species info, reduced to species with trait data and relevant cols: 
head(spinf)
head(spread)

# Calculate continuous trait variables, for each seeded species (s), in each year, in each plot:

# 1. RAW = actual trait value

# 2. GAP = gap size; nearest neighbour high - nearest neighbour low

# 3. relFD = relative functional difference; trait(s) - trait(community weighted mean):

# 4. relNN = relative neighbour; trait(s) - trait(nearest neighbour):

# Add categorical trait data to spread data. This uses the ft data frame which was created for the colonisation functional distance:
head(ft)
head(spread)

spread<-merge(spread, ft, by="sp", all.x=T, all.y=F)
head(spread)

# check:
sp.now<-sample(levels(spread$sp),1)
ftraits[ftraits$sp==sp.now,]
head(spread[spread$sp==sp.now, ])

# New variables:

# The prop_native and prop_intro variables both relate to the proportion native, the difference is where the unknown species were grouped (into native in prop_native and introduced in prop_intro):
spread$prop_peren<-NA
spread$prop_grass<-NA
spread$prop_legume<-NA

spread$hgt_cwm<-NA
spread$seed_cwm<-NA
spread$ldmc_cwm<-NA
spread$sla_cwm<-NA

spread$hgt_raw<-NA
spread$seed_raw<-NA
spread$ldmc_raw<-NA
spread$sla_raw<-NA

spread$hgt_gap<-NA
spread$seed_gap<-NA
spread$ldmc_gap<-NA
spread$sla_gap<-NA

spread$hgt_relFD<-NA
spread$seed_relFD<-NA
spread$ldmc_relFD<-NA
spread$sla_relFD<-NA

spread$hgt_relNN<-NA
spread$seed_relNN<-NA
spread$ldmc_relNN<-NA
spread$sla_relNN<-NA

traits<-c("hgt","seed","ldmc","sla")

# make a vector of sp for separating woodies:
woodies<-as.character(sp_info$sp[sp_info$func_grp=="W"])

# Use spr.abund for the community data because it is formatted for the plot level:

# We will use plot level cover for community weighted means, so standardise the value column at the plot level.

# Remove the standardised cover column as this was standardised to subplot and will not be used:
spabund<-spr.abund
spabund$std_cover<-NULL

# Aggregate the value cover to the sp x plot level, removing duplicate records for sp x plot within a year:
spabund<-aggregate(value~sp+year+plot,data=spabund,FUN=sum)
spabund$yr_plot<-paste(spabund$year, spabund$plot, sep="_")
head(spabund)

# Standardise to plot level for each sp x plot:
yrpls<-unique(spabund$yr_plot)
outl1<-list()

for (i in 1:length(yrpls)){
yrpl.thisrun<-yrpls[i]
data.thisrun<-spabund[spabund$yr_plot==yrpl.thisrun,]
data.thisrun$stdcov_plot<-data.thisrun$value/sum(data.thisrun$value)
outl1[[i]]<-data.thisrun
}
spab<-data.frame(do.call(rbind,outl1))
head(spab)
head(spabund)

# check:
rowx<-spab[sample(1:length(spab[,1]),1),]
sum(spab[spab$year==rowx$year & spab$plot==rowx$plot,]$stdcov_plot)
# Should be perfectly linear within sp x plot:
# plot(spab[spab$year==rowx$year & spab$plot==rowx$plot,]$value,spab[spab$year==rowx$year & spab$plot==rowx$plot,]$stdcov_plot)

### --- START CALCULATIONS:

for (i in 1:length(spread[,1])){

# define species, subplot and year for each row:
sp<-as.character(spread$sp[i])
plot<-as.character(spread$plot[i])
year<-as.character(spread$year[i])

# get community abundance data for the current subplot and year:
comm.data<-spab[spab$plot==plot & spab$year==year,]

# remove the target species from the community data:
if(sp %in% comm.data$sp==T) comm.data<-comm.data[-which(comm.data$sp==sp),]

# add continuous trait data to the community data:
comm.data<-merge(comm.data,spinf,by="sp",all.x=T,all.y=F)

# the spinf data used in the prev line has been subsetted on species that have continuous traits. For categorical trait data for the community, use the ft data set
comm.data<-merge(comm.data,ft,by="sp",all.x=T,all.y=F)

# get continuous trait data for the current species:
tr.data<-spinf[spinf$sp==sp,]

# add categorical trait data for the current species:
tr.data<-merge(tr.data,ft,by="sp",all.x=T,all.y=F)

head(spread,3)
comm.data

# Calculate values of categorical traits relative to the community:
comm.cover<-sum(comm.data$stdcov_plot)

spread$prop_peren[i]<-sum(comm.data$stdcov_plot[comm.data$lifespan=="P"])/comm.cover
spread$prop_grass[i]<-sum(comm.data$stdcov_plot[comm.data$grass==1])/comm.cover
spread$prop_legume[i]<-sum(comm.data$stdcov_plot[comm.data$legume==1])/comm.cover

no.incomm<-length(comm.data[,1])

spread[(i-1):(i+1),]

comm.data
tr.data

# calculate FD measures only when there is trait data for the current species:

if(length(tr.data[,1])>0){

# COMMUNITY WEIGHTED MEAN:
spread$hgt_cwm[i]<-weighted.mean(comm.data$hgt_mean,comm.data$stdcov_plot, na.rm=T)
spread$ldmc_cwm[i]<-weighted.mean(comm.data$ldmc_mean,comm.data$stdcov_plot, na.rm=T)
spread$sla_cwm[i]<-weighted.mean(comm.data$sla_mean,comm.data$stdcov_plot, na.rm=T)

# RAW:
spread$hgt_raw[i]<-tr.data$hgt_mean
spread$seed_raw[i]<-tr.data$seed_mean
spread$ldmc_raw[i]<-tr.data$ldmc_mean
spread$sla_raw[i]<-tr.data$sla_mean

# relative FD:
spread$hgt_relFD[i]<-tr.data$hgt_mean-weighted.mean(comm.data$hgt_mean,comm.data$stdcov_plot, na.rm=T)
spread$ldmc_relFD[i]<-tr.data$ldmc_mean-weighted.mean(comm.data$ldmc_mean,comm.data$stdcov_plot, na.rm=T)
spread$sla_relFD[i]<-tr.data$sla_mean-weighted.mean(comm.data$sla_mean,comm.data$stdcov_plot, na.rm=T)

# Use "seed.comm" for ALL seed calculations (same as comm.data but with woodies removed where present):

# If there are woodies present, remove them and tidy seed.comm:
if(length(which(comm.data$sp %in% woodies))>0){
seed.comm<-comm.data
seed.comm<-seed.comm[-which(seed.comm$sp %in% woodies),]
seed.comm<-tidy.df(seed.comm)
spread$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$stdcov_plot, na.rm=T)
spread$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$stdcov_plot, na.rm=T)
}
# If there are no woodies, make seed.comm == comm.data. We will use seed.comm in the NN and gap calcs, so define it here regardless of whether there are woodies or not:
if(length(which(comm.data$sp %in% woodies))==0){
seed.comm<-comm.data
spread$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$stdcov_plot, na.rm=T)
spread$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$stdcov_plot, na.rm=T)
}

# relative NN and GAP, one trait at a time (j):

spread[(i-1):(i+1),]

# store comm.data in a separate data frame for updating in the j trait loop:
comm.data.store<-comm.data

for (j in 1:length(traits)){

trait.thisrun<-traits[j]

# update comm.data, depending on the trait:

if (trait.thisrun=="seed") comm.data<-seed.comm
if (trait.thisrun!="seed") comm.data<-comm.data.store

# get trait data for the current community and target species and order them by trait value:
nn<-data.frame(tr=comm.data[,which(colnames(comm.data)==paste(trait.thisrun,"mean",sep="_"))],type="comm",stringsAsFactors=F)
nn<-rbind(nn,c(tr.data[,which(colnames(tr.data)==paste(trait.thisrun,"mean",sep="_"))],"target"))
nn$tr<-as.numeric(nn$tr)
nn<-nn[order(nn$tr),]
nn<-tidy.df(nn)
nn

# if an NA belongs to the target species, skip to the next trait j:
if(is.na(nn$tr[nn$type=="target"])) next

# if an NA belongs to someone else in the community, remove NAs:
if(length(which(is.na(nn$tr)))>0) nn<-nn[-which(is.na(nn$tr)),]
nn

# identify the row occupied by the target:
target.ind<-which(nn$type=="target")

# **** SINGLE NEIGHBOUR (smallest & largest):

# If the species is the SMALLEST, take the next smallest as the NN and the full range as the gap:
if(target.ind==1) {

## *** NEAREST NEIGHBOUR:	

spread[i,which(colnames(spread)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# nn_distance is the distance between the invader who happens to be the smallest for the current trait in the current community and its nearest, next biggest neighbour:
nn_distance<-nn$tr[target.ind+1]-nn$tr[target.ind]

# invader_mean is the mean value for the invader for the current trait:
invader_mean<-tr.data[,grep(trait.thisrun,colnames(tr.data))]

# If NN distance > invader mean trait value, then gap size = NN trait value, otherwise gap = NN distance*2

if (nn_distance > invader_mean) spread[i,which(colnames(spread)==paste(trait.thisrun,"gap",sep="_"))]<-nn$tr[target.ind+1] else spread[i,which(colnames(spread)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind+1]-nn$tr[target.ind])*2

} # close if smallest

spread[(i-1):(i+1),]

# If the species is the LARGEST, take the next largest as the NN and the full range as the gap:
if(target.ind==length(nn[,1])) {
	
## *** NEAREST NEIGHBOUR:	
spread[i,which(colnames(spread)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# gap = NN distance*2

spread[i,which(colnames(spread)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind]-nn$tr[target.ind-1])*2

} # close if largest

# **** DOUBLE NEIGHBOUR:

spread[(i-1):(i+1),]

# calculate the distance between the neighbour on the small side and the neighbour on the big side:
if(target.ind>1 & target.ind<length(nn[,1])){
lowNN.dist<-nn$tr[target.ind]-nn$tr[target.ind-1]
highNN.dist<-nn$tr[target.ind+1]-nn$tr[target.ind]

# we want relative NN distance, so each calculation should be the nearest neighbour minus the target species . The code must preserve the sign, rather than just take the smallest distance:

# 21 March 2016 update: switched the sides of the calculation so that, for example, the taller nearest neighbour would be positive, the shorter nearest neighbour would be negative. This should be easier to interpret. Corresponding changes have been made for single (largest and smallest) neighbour above.

if(lowNN.dist < highNN.dist) spread[i,which(colnames(spread)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

if(lowNN.dist > highNN.dist) spread[i,which(colnames(spread)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

# When low and high neighbour have the same distance, we were getting NAs for relNN values. The next line deals with those cases, and we're using a positive value, from JC email (21 March 2016): Let's use a positive value when the upper and lower NN are the same.  
if(lowNN.dist == highNN.dist) spread[i,which(colnames(spread)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

spread[(i-1):(i+1),]

# gap distance:
spread[i,which(colnames(spread)==paste(trait.thisrun,"gap",sep="_"))]<-lowNN.dist + highNN.dist

} # close if species has two neighbours

} # close j traits

} # close if species has trait data

} # close i rows of data frame

# Add absolute FD:
spread$hgt_absFD<-abs(spread$hgt_relFD)
spread$seed_absFD<-abs(spread$seed_relFD)
spread$ldmc_absFD<-abs(spread$ldmc_relFD)
spread$sla_absFD<-abs(spread$sla_relFD)

# check:
check.rows(spread,3)

} # close functional distance

#### ------ ADD ENVIRO VARIABLES:
{

# All variables will need to be added at the plot level. 
head(spread,3)

## ~~ ADD: no. planted (PLOT level). np1 contains the number of species seeded for each seeded subplot, excluding those two species with unviable seed. It was defined in the "Tidy subplot & species data" section of "Data setup". Modify this so it works at the plot level (Total number seeded in the whole plot (i.e. across both of the seeded subplots):

np4<-np1

# Add plot and subplot to np4
np4$plot<-substr(np4$plsub,1,unlist(gregexpr("[A-Z]",np4$plsub))-1)
np4$subplot<-substr(np4$plsub,unlist(gregexpr("[A-Z]",np4$plsub)),unlist(gregexpr("[A-Z]",np4$plsub)))
head(np4)

# Get total no_seeded across both subplots:
np5<-data.frame(aggregate(no_seeded~plot,sum,data=np4))

head(np5)
head(spread,3)

spread<-merge(spread,np5,by="plot",all.x=T, all.y=F)

# minus one for the target invader:
spread$no_seeded<-spread$no_seeded-1
head(spread,3)

## ~~ ADD: NO3 (take mean for each plot across subplots):
no3m<-soil.mean[,c("plot","NO3soil","moisture")]
no3.ag<-aggregate(NO3soil~plot,FUN=mean,data=no3m)
head(no3.ag)
spread<-merge(spread, no3.ag, by="plot", all.x=T, all.y=F)
head(spread,3)

} # close enviro variables

#### ------ TIDY AND SCALE:
{

head(spread,3)

## ~~ Order by plot, subplot, year:
spread$plot<-as.numeric(spread$plot)
spread<-spread[order(spread$plot, spread$year),]
spread<-tidy.df(spread)

## ~~ Factorise:
spread$plot<-as.factor(spread$plot)

# Make year numeric, beginning at zero:
spread$yr<-spread$year-min(spread$year)

## ~~ Scale ALL numeric variables (exclude responses, raw year, numeric year, factors and functional dummies):

to_scale<-sapply(spread,is.numeric)
to_scale [names(to_scale) %in% c("spread","spr_abund_raw","std_cover","year","yr","grass","legume")]<-FALSE
to_scale

sprScale<-spread
sprScale[,to_scale]<-scale(sprScale[,to_scale])
sapply(sprScale[,to_scale],range,na.rm=T)
head(sprScale,4)

# Remove 'woody' invaders
spreadSc<-sprScale[-which(sprScale$func_grp=="W"),]
spreadSc<-tidy.df(spreadSc)
check.rows(spreadSc,2)

# Make a new data frame from the colonise (unscaled) base data, with woodies removed and functional groups added as dummies. This will be directly comparable to coloniseSc, and can be used to plot the x-axis using unscaled variables:
head(spread,2)

spr.unscaled<-spread

# Remove 'woody' invaders
spr.unscaled<-spr.unscaled[-which(spr.unscaled$func_grp=="W"),]
spr.unscaled<-tidy.df(spr.unscaled)
check.rows(spr.unscaled,2)

# Double check that none of the data sets have NAs:
# For the new spread abundance analysis, I'm keeping the NAs for now
head(spread,3)
head(sprScale,3)
head(spreadSc,3)
head(spr.unscaled,3)

na.spread<-data.frame(pred=names(apply(spread,2,function(x) length(which(is.na(x))))),no_NAs=apply(spread,2,function(x) length(which(is.na(x)))))
na.spread<-tidy.df(na.spread)

na.sprScale<-data.frame(pred=names(apply(sprScale,2,function(x) length(which(is.na(x))))),no_NAs=apply(sprScale,2,function(x) length(which(is.na(x)))))
na.sprScale<-tidy.df(na.sprScale)

na.spreadSc<-data.frame(pred=names(apply(spreadSc,2,function(x) length(which(is.na(x))))),no_NAs=apply(spreadSc,2,function(x) length(which(is.na(x)))))
na.spreadSc<-tidy.df(na.spreadSc)

na.spr.unscaled<-data.frame(pred=names(apply(spr.unscaled,2,function(x) length(which(is.na(x))))),no_NAs=apply(spr.unscaled,2,function(x) length(which(is.na(x)))))
na.spr.unscaled<-tidy.df(na.spr.unscaled)

# These should all be zero:
length(which(na.spread$no_NAs>0))
length(which(na.sprScale$no_NAs>0))
length(which(na.spreadSc$no_NAs>0))
length(which(na.spr.unscaled$no_NAs>0))

# Re-level lifespan to P in ALL data sets
spread$lifespan<-relevel(spread$lifespan,ref="P")
sprScale$lifespan<-relevel(sprScale$lifespan,ref="P")
spreadSc$lifespan<-relevel(spreadSc$lifespan,ref="P")
spr.unscaled$lifespan<-relevel(spr.unscaled$lifespan,ref="P")

} # close tidy and scale

} # close probability of spread

######################################
#  	  COLONISATION ABUNDANCE		 #
######################################
{

#### ------ BASE COLONISE ABUND DATA:
{	
# Colonisation abundance for all seeded species, removing species that were anywhere in the plot in 1991:

# seeded species (the two species with non-viable seed have already been removed):
sd.sp
head(s)

# 54 seeded subplots
sd.sub

# species data
head(sp_dat)

# survey years
syr<-unique(sp_dat$year)

# updated 7 Feb 2017: we were previously using 'seeded' as the base data, containing every seeded subplot and the species that were seeded including those present in 1991. However, we are now updating this to remove species x subplot combinations where that species was anywhere on the plot in 1991. That is, we are now using 'newcomers only' for the colonisation abundance analysis. Thus, use 's' as the base data for colab.

# If the data are correct, Ast_nov and Mim_rin should not appear and there should be only 52 seeded species:
which(s$sp=="Ast_nov"); which(s$sp=="Mim_rin"); length(levels(s$sp))

# Remove satellite plots
colab.abund<-sp_dat[sp_dat$subplot %in% c("W","X","Y","Z"),]
colab.abund<-tidy.df(colab.abund)

# remove species name:
colab.abund<-colab.abund[,-which(colnames(colab.abund) %in% c("species"))]

# remove rows that are not seeded species from abundance data:
colab.abund<-colab.abund[colab.abund$sp %in% seeded$sp,]
colab.abund<-tidy.df(colab.abund)
if(length(which(!colab.abund$sp %in% seeded$sp))>0) print("data subset is wrong")

# remove value column:
colab.abund<-colab.abund[,-which(colnames(colab.abund)=="value")]

# 's' already has plot and subplot removed for for merging, so just update its name here:
colab.seeded<-s
colab.abund<-colab.abund[,c("sp","year","std_cover","plsub")]

head(colab.seeded)
head(colab.abund)

# Check it:
sp.now<-sample(levels(colab.seeded$sp),1)
ab91<-colab.abund[colab.abund$year==1991,]
ab91<-tidy.df(ab91)
plsub.abnow<-as.character(ab91$plsub[ab91$sp==sp.now])
plsub.sdnow<-as.character(colab.seeded$plsub[colab.seeded$sp==sp.now])
# There should never be plsub.sdnow in plsub.abnow in 1991:
which(plsub.sdnow %in% plsub.abnow==T)

# if we run this from syr==1 we get zeros for all species x subplot combinations, so just run from syr==2 as with the colonisation probability analysis.

# list for the output:
outlist.colab<-list()

for (i in 2:length(syr)){

# do each year separately:
data.thisrun<-colab.abund[colab.abund$year==syr[i],]

# remove rows that are not seeded subplots:
data.thisrun<-data.thisrun[data.thisrun$plsub %in% sd.sub,]
data.thisrun<-tidy.df(data.thisrun)
head(data.thisrun)
head(colab.seeded)

# add std_cover to base data:
mg.dat<-merge(colab.seeded,data.thisrun,by=c("sp","plsub"),all.x=T,all.y=F)
head(mg.dat)

# update year
mg.dat$year<-syr[i]

# Make NAs zero:
mg.dat$std_cover[which(is.na(mg.dat$std_cover))]<-0

# add to list:
outlist.colab[[i]]<-mg.dat

} # close i survey year

colab<-data.frame(do.call(rbind, outlist.colab))

# random checks:
sp.now<-sample(sd.sp$sp,1)
sub.now<-sample(sd.sub,1)
yr.now<-sample(syr,1)
colab[colab$plsub==sub.now & colab$year==yr.now,]
colab.abund[colab.abund$plsub==sub.now & colab.abund$year==yr.now,]
colab.seeded[colab.seeded$plsub==sub.now,]

# add plot and subplot back in (removed in loop for merging):
colab$plot<-substr(colab$plsub,1,unlist(gregexpr("[A-Z]",colab$plsub))-1)
colab$subplot<-substr(colab$plsub,unlist(gregexpr("[A-Z]",colab$plsub)),unlist(gregexpr("[A-Z]",colab$plsub)))
head(colab)

# Remove zeros:
colab<-colab[-which(colab$std_cover==0),]
colab<-tidy.df(colab)

# Transform. No need to do this now that we're using beta models
# colab$ab_logit<-log(colab$std_cover /(1-colab$std_cover))

head(colab)

} # close base colonise abundance

#### ------ FUNCTIONAL DISTANCE:
{
# spinf is the species info, reduced to species with trait data and relevant cols: 
head(spinf)
head(colab)

# Calculate continuous trait variables, for each seeded species (s), in each year, in each plot:

# 1. RAW = actual trait value

# 2. GAP = gap size; nearest neighbour high - nearest neighbour low

# 3. relFD = relative functional difference; trait(s) - trait(community weighted mean):

# 4. relNN = relative neighbour; trait(s) - trait(nearest neighbour):

# Add categorical trait data to spread data. This uses the ft data frame which was created for the colonisation functional distance:
head(ft)
head(colab)	

colab<-merge(colab, ft, by="sp", all.x=T, all.y=F)
head(colab)

# check:
sp.now<-sample(levels(colab$sp),1)
ftraits[ftraits$sp==sp.now,]
head(colab[colab$sp==sp.now, ])

# New variables:

# The prop_native and prop_intro variables both relate to the proportion native, the difference is where the unknown species were grouped (into native in prop_native and introduced in prop_intro):
colab$prop_peren<-NA
colab$prop_grass<-NA
colab$prop_legume<-NA

colab$hgt_cwm<-NA
colab$seed_cwm<-NA
colab$ldmc_cwm<-NA
colab$sla_cwm<-NA

colab$hgt_raw<-NA
colab$seed_raw<-NA
colab$ldmc_raw<-NA
colab$sla_raw<-NA

colab$hgt_gap<-NA
colab$seed_gap<-NA
colab$ldmc_gap<-NA
colab$sla_gap<-NA

colab$hgt_relFD<-NA
colab$seed_relFD<-NA
colab$ldmc_relFD<-NA
colab$sla_relFD<-NA

colab$hgt_relNN<-NA
colab$seed_relNN<-NA
colab$ldmc_relNN<-NA
colab$sla_relNN<-NA

traits<-c("hgt","seed","ldmc","sla")

# make a vector of sp for separating woodies:
woodies<-as.character(sp_info$sp[sp_info$func_grp=="W"])

# Use abund for the community data, formatted for the subplot level:
head(abund)

### --- START CALCULATIONS:

for (i in 1:length(colab[,1])){

# define species, subplot and year for each row:
sp<-as.character(colab$sp[i])
plsub<-as.character(colab$plsub[i])
year<-as.character(colab$year[i])

# get community abundance data for the current subplot and year:
comm.data<-abund[abund$plsub==plsub & abund$year==year,]

# remove the target species from the community data:
if(sp %in% comm.data$sp==T) comm.data<-comm.data[-which(comm.data$sp==sp),]

# add continuous trait data to the community data:
comm.data<-merge(comm.data,spinf,by="sp",all.x=T,all.y=F)

# the spinf data used in the prev line has been subsetted on species that have continuous traits. For categorical trait data for the community, use the ft data set
comm.data<-merge(comm.data,ft,by="sp",all.x=T,all.y=F)

# get continuous trait data for the current species:
tr.data<-spinf[spinf$sp==sp,]

# add categorical trait data for the current species:
tr.data<-merge(tr.data,ft,by="sp",all.x=T,all.y=F)

head(colab,3)
comm.data

# Calculate values of categorical traits relative to the community:
comm.cover<-sum(comm.data$std_cover)

colab$prop_peren[i]<-sum(comm.data$std_cover[comm.data$lifespan=="P"])/comm.cover
colab$prop_grass[i]<-sum(comm.data$std_cover[comm.data$grass==1])/comm.cover
colab$prop_legume[i]<-sum(comm.data$std_cover[comm.data$legume==1])/comm.cover

no.incomm<-length(comm.data[,1])

colab[(i-1):(i+1),]

comm.data
tr.data

# calculate FD measures only when there is trait data for the current species:

if(length(tr.data[,1])>0){

# COMMUNITY WEIGHTED MEAN:
colab$hgt_cwm[i]<-weighted.mean(comm.data$hgt_mean,comm.data$std_cover, na.rm=T)
colab$ldmc_cwm[i]<-weighted.mean(comm.data$ldmc_mean,comm.data$std_cover, na.rm=T)
colab$sla_cwm[i]<-weighted.mean(comm.data$sla_mean,comm.data$std_cover, na.rm=T)

# RAW:
colab$hgt_raw[i]<-tr.data$hgt_mean
colab$seed_raw[i]<-tr.data$seed_mean
colab$ldmc_raw[i]<-tr.data$ldmc_mean
colab$sla_raw[i]<-tr.data$sla_mean

# relative FD:
colab$hgt_relFD[i]<-tr.data$hgt_mean-weighted.mean(comm.data$hgt_mean,comm.data$std_cover, na.rm=T)
colab$ldmc_relFD[i]<-tr.data$ldmc_mean-weighted.mean(comm.data$ldmc_mean,comm.data$std_cover, na.rm=T)
colab$sla_relFD[i]<-tr.data$sla_mean-weighted.mean(comm.data$sla_mean,comm.data$std_cover, na.rm=T)

# Use "seed.comm" for ALL seed calculations (same as comm.data but with woodies removed where present):

# If there are woodies present, remove them and tidy seed.comm:
if(length(which(comm.data$sp %in% woodies))>0){
seed.comm<-comm.data
seed.comm<-seed.comm[-which(seed.comm$sp %in% woodies),]
seed.comm<-tidy.df(seed.comm)
colab$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
colab$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
}
# If there are no woodies, make seed.comm == comm.data. We will use seed.comm in the NN and gap calcs, so define it here regardless of whether there are woodies or not:
if(length(which(comm.data$sp %in% woodies))==0){
seed.comm<-comm.data
colab$seed_relFD[i]<-tr.data$seed_mean-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
colab$seed_cwm[i]<-weighted.mean(seed.comm$seed_mean,seed.comm$std_cover, na.rm=T)
}

# relative NN and GAP, one trait at a time (j):

colab[(i-1):(i+1),]

# store comm.data in a separate data frame for updating in the j trait loop:
comm.data.store<-comm.data

for (j in 1:length(traits)){

trait.thisrun<-traits[j]

# update comm.data, depending on the trait:

if (trait.thisrun=="seed") comm.data<-seed.comm
if (trait.thisrun!="seed") comm.data<-comm.data.store

# get trait data for the current community and target species and order them by trait value:
nn<-data.frame(tr=comm.data[,which(colnames(comm.data)==paste(trait.thisrun,"mean",sep="_"))],type="comm",stringsAsFactors=F)
nn<-rbind(nn,c(tr.data[,which(colnames(tr.data)==paste(trait.thisrun,"mean",sep="_"))],"target"))
nn$tr<-as.numeric(nn$tr)
nn<-nn[order(nn$tr),]
nn<-tidy.df(nn)
nn

# if an NA belongs to the target species, skip to the next trait j:
if(is.na(nn$tr[nn$type=="target"])) next

# if an NA belongs to someone else in the community, remove NAs:
if(length(which(is.na(nn$tr)))>0) nn<-nn[-which(is.na(nn$tr)),]
nn

# identify the row occupied by the target:
target.ind<-which(nn$type=="target")

# **** SINGLE NEIGHBOUR (smallest & largest):

# If the species is the SMALLEST, take the next smallest as the NN and the full range as the gap:
if(target.ind==1) {

## *** NEAREST NEIGHBOUR:	
colab[i,which(colnames(colab)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# nn_distance is the distance between the invader who happens to be the smallest for the current trait in the current community and its nearest, next biggest neighbour:
nn_distance<-nn$tr[target.ind+1]-nn$tr[target.ind]

# invader_mean is the mean value for the invader for the current trait:
invader_mean<-tr.data[,grep(trait.thisrun,colnames(tr.data))]

# If NN distance > invader mean trait value, then gap size = NN trait value, otherwise gap = NN distance*2

if (nn_distance > invader_mean) colab[i,which(colnames(colab)==paste(trait.thisrun,"gap",sep="_"))]<-nn$tr[target.ind+1] else colab[i,which(colnames(colab)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind+1]-nn$tr[target.ind])*2

} # close if smallest

colab[(i-1):(i+1),]

# If the species is the LARGEST, take the next largest as the NN and the full range as the gap:
if(target.ind==length(nn[,1])) {

## *** NEAREST NEIGHBOUR:	

colab[i,which(colnames(colab)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

## *** GAP:	

# This is the method using 2 x the gap between the nearest neighbour:

# gap = NN distance*2

colab[i,which(colnames(colab)==paste(trait.thisrun,"gap",sep="_"))]<-(nn$tr[target.ind]-nn$tr[target.ind-1])*2

} # close if largest

# **** DOUBLE NEIGHBOUR:

colab[(i-1):(i+1),]

# calculate the distance between the neighbour on the small side and the neighbour on the big side:
if(target.ind>1 & target.ind<length(nn[,1])){
lowNN.dist<-nn$tr[target.ind]-nn$tr[target.ind-1]
highNN.dist<-nn$tr[target.ind+1]-nn$tr[target.ind]

# we want relative NN distance, so each calculation should be the nearest neighbour minus the target species . The code must preserve the sign, rather than just take the smallest distance:

# 21 March 2016 update: switched the sides of the calculation so that, for example, the taller nearest neighbour would be positive, the shorter nearest neighbour would be negative. This should be easier to interpret. Corresponding changes have been made for single (largest and smallest) neighbour above.

if(lowNN.dist < highNN.dist) colab[i,which(colnames(colab)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind-1]-nn$tr[target.ind]

if(lowNN.dist > highNN.dist) colab[i,which(colnames(colab)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

# When low and high neighbour have the same distance, we were getting NAs for relNN values. The next line deals with those cases, and we're using a positive value, from JC email (21 March 2016): Let's use a positive value when the upper and lower NN are the same.  
if(lowNN.dist == highNN.dist) colab[i,which(colnames(colab)==paste(trait.thisrun,"relNN",sep="_"))]<-nn$tr[target.ind+1]-nn$tr[target.ind]

colab[(i-1):(i+1),]

# gap distance:
colab[i,which(colnames(colab)==paste(trait.thisrun,"gap",sep="_"))]<-lowNN.dist + highNN.dist

} # close if species has two neighbours

} # close j traits

} # close if species has trait data

} # close i rows of data frame

# Add absolute FD:
colab$hgt_absFD<-abs(colab$hgt_relFD)
colab$seed_absFD<-abs(colab$seed_relFD)
colab$ldmc_absFD<-abs(colab$ldmc_relFD)
colab$sla_absFD<-abs(colab$sla_relFD)

# check:
check.rows(colab,3)

} # close functional distance

#### ------ ADD ENVIRO VARIABLES:
{

check.rows(colab,3)

## ~~ ADD: no. planted (subplot level). np1 contains the number of species seeded for each seeded subplot, excluding those two species with unviable seed. It was defined in the "Tidy subplot & species data" section of "Data setup":
head(np1)
head(colab,3)

colab<-merge(colab,np1,by="plsub",all.x=T, all.y=F)

# minus one for the target invader:
colab$no_seeded<-colab$no_seeded-1
head(colab,3)

## ~~ ADD: elevation (plot level):
tp<-topog[,c("plot","masl")]
colab<-merge(colab,tp,by="plot",all.x=T, all.y=F)

## ~~ ADD: NO3 and MOISTURE (subplot level):
no3<-soil.mean[,c("plsub","NO3soil","moisture")]
colab<-merge(colab, no3, by="plsub", all.x=T, all.y=F)

## ~~ ADD: BURN and PRECIP (YEAR level):
head(burn_precip)
colab<-merge(colab, burn_precip, by="year", all.x=T, all.y=F)
head(colab,3)

## ~~ ADD: SOIL DISTURBANCE (subplot x year):
# From JC: Lets just use  bare ground by itself, and all three of the other temporal vars ("dsince_fire","rfall_mm","year").  

# Tidy disturbance data:
soil_dist$plsub<-paste(soil_dist$plot,soil_dist$subplot,sep="")
plsubs<-c("W","X","Y","Z")
soil_dist<-soil_dist[soil_dist$subplot %in% plsubs,]
soil_dist<-tidy.df(soil_dist)
head(soil_dist)

# Summarise bare ground disturbance:
soil.bg<-soil_dist[soil_dist$disturb=="bare",]
soil.bg<-tidy.df(soil.bg)
head(soil.bg)
dist1<-aggregate(std_cover~plsub+year,sum,data=soil.bg)
head(dist1)

# Check
yr.now<-sample(soil_dist$year,1)
plsub.now<-sample(soil_dist$plsub,1)
soil_dist[soil_dist$year==yr.now & soil_dist$plsub==plsub.now,]
dist1[dist1$year==yr.now & dist1$plsub==plsub.now,]

# Change the colname to soil_dist:
colnames(dist1)[which(colnames(dist1)=="std_cover")]<-"soil_dist"

# Merge with main data:
colab<-merge(colab,dist1,by=c("year","plsub"),all.x=T,all.y=F)
colab[sample(1:length(colab[,1]),size=4),]

# Make soil disturbance NAs zero (nothing was recorded):
colab$soil_dist[which(is.na(colab$soil_dist))]<-0

check.rows(colab,3)

} # close enviro variables

#### ------ TIDY AND SCALE:
{

check.rows(colab,3)
head(colab,3)

## ~~ Order by plot, subplot, year:
colab$plot<-as.numeric(colab$plot)
colab<-colab[order(colab$plot, colab$subplot, colab$year),]
colab<-tidy.df(colab)

## ~~ Factorise:
colab$plot<-as.factor(colab$plot)
colab$subplot<-as.factor(colab$subplot)

# Make year numeric, beginning at zero:
colab$yr<-colab$year-min(colab$year)

## ~~ Scale ALL numeric variables (exclude response, raw year, numeric year, factors and functional dummies):
to_scale<-sapply(colab,is.numeric)
to_scale [names(to_scale) %in% c("std_cover","year","yr","grass","legume")]<-FALSE
to_scale

cabScale<-colab
cabScale[,to_scale]<-scale(cabScale[,to_scale])
sapply(cabScale[,to_scale],range,na.rm=T)
head(cabScale,3)

# Remove 'woody' invaders
colabSc<-cabScale[-which(cabScale$func_grp=="W"),]
colabSc<-tidy.df(colabSc)
check.rows(colabSc,2)

# Make a new data frame from the colab (unscaled) base data, with woodies removed and functional groups added as dummies. This will be directly comparable to colabSc, and can be used to plot the x-axis using unscaled variables:
head(colab,2)

cab.unscaled<-colab

# Remove 'woody' invaders
cab.unscaled<-cab.unscaled[-which(cab.unscaled$func_grp=="W"),]
cab.unscaled<-tidy.df(cab.unscaled)
check.rows(cab.unscaled,2)

# Double check that none of the data sets have NAs:
head(colab,3)
head(cabScale,3)
head(colabSc,3)

na.colab<-data.frame(pred=names(apply(colab,2,function(x) length(which(is.na(x))))),no_NAs=apply(colab,2,function(x) length(which(is.na(x)))))
na.colab<-tidy.df(na.colab)

na.cabScale<-data.frame(pred=names(apply(cabScale,2,function(x) length(which(is.na(x))))),no_NAs=apply(cabScale,2,function(x) length(which(is.na(x)))))
na.cabScale<-tidy.df(na.cabScale)

na.colabSc<-data.frame(pred=names(apply(colabSc,2,function(x) length(which(is.na(x))))),no_NAs=apply(colabSc,2,function(x) length(which(is.na(x)))))
na.colabSc<-tidy.df(na.colabSc)

# These should all be zero:
length(which(na.colab$no_NAs>0))
length(which(na.cabScale$no_NAs>0))
length(which(na.colabSc$no_NAs>0))

# Re-level lifespan to P in ALL data sets
colab$lifespan<-relevel(colab$lifespan,ref="P")
cabScale$lifespan<-relevel(cabScale$lifespan,ref="P")
colabSc$lifespan<-relevel(colabSc$lifespan,ref="P")
cab.unscaled$lifespan<-relevel(cab.unscaled$lifespan,ref="P")

} # close tidy and scale

} # close colonisation abundance



