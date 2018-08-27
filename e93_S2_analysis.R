#------------------------------------#
#------~~~~~ SECTION 2 ~~~~~---------#
#-----~~~~~ E93 ANALYSIS ~~~~~-------#
#------------------------------------#

# Authors: Annabel Smith & Jane Catford

# *** LOAD WORKSPACE:
load("e93_wksp_aug18")

# load functions and libraries:
source("/Users/annabelsmith/Documents/01_Current/PROJECTS/03_CATFORD_invasion_traits/DATA/AS_data_and_analysis/ANALYSIS_files/e93_script/e93_S0_function_library.R")

library("lme4")
library("AICcmodavg")
library("MuMIn")
library(arm)
library("aods3")
library(glmmADMB)
library(piecewiseSEM)

# We are using piecewiseSEM to calculate the r2m and r2c for all models. It is much quicker than the function in MuMIn which seems to refit the model before doing the calculation. I have checked both functions against hand-calculations and they agree. I have also checked them both using the examples in Johnson 2014 MEE which explains the situation specific to random slopes models and they give the same result.

######################################
# 		   DATA EXPLANATIONS		 #
######################################
{
#### ------ COLONISATION & SPREAD:

# full data set, including woodies, variables unscaled
head(colonise,3); dim(colonise) 
head(spread,3); dim(spread)
head(colab,3); dim(colab)
# head(impact,3); dim(impact)

# full data set, including woodies, variables scaled
head(colScale,3); dim(colScale) 
head(sprScale,3); dim(sprScale) 
head(cabScale,3); dim(cabScale)
# head(impScale,3); dim(impScale)

# tidy data set, woodies removed, variables scaled. These are the data to model:
head(coloniseSc,3); dim(coloniseSc) 
head(spreadSc,3); dim(spreadSc) 
head(colabSc,3); dim(colabSc)

# tidy data set, woodies removed, variables UNSCALED. This is identical to coloniseSc/spreadSc, but unscaled so can be used for plotting numbers on the x-axis:
head(col.unscaled,3); dim(col.unscaled) 
head(spr.unscaled,3); dim(spr.unscaled) 
head(cab.unscaled,3); dim(cab.unscaled)

} # close data explanations

## ~~~~ ---- ** ANALYSIS ** ---- ~~~~ ##

## ---- STEP 0: LOAD MODEL TABLES AND REMOVE TERMS WITH VIF > 3:
{

mod.dir<-"../e93_model_tables"

col.mods1.allterms<-read.table(paste(mod.dir,"model_table_STEP1_col.txt",sep="/"),header=T)
spr.mods1.allterms<-read.table(paste(mod.dir,"model_table_STEP1_spr.txt",sep="/"),header=T)
cab.mods1.allterms<-read.table(paste(mod.dir,"model_table_STEP1_colab.txt",sep="/"),header=T)
sab.mods1.allterms<-read.table(paste(mod.dir,"model_table_STEP1_sprab.txt",sep="/"),header=T)


col.mods2.allterms<-read.table(paste(mod.dir,"model_table_STEP2_col.txt",sep="/"),header=T)
spr.mods2.allterms<-read.table(paste(mod.dir,"model_table_STEP2_spr.txt",sep="/"),header=T)
cab.mods2.allterms<-read.table(paste(mod.dir,"model_table_STEP2_colab.txt",sep="/"),header=T)
sab.mods2.allterms<-read.table(paste(mod.dir,"model_table_STEP2_sprab.txt",sep="/"),header=T)

# Update the following and run separately for each data set:
model.table<-sab.mods2.allterms
dataset<-sprabSc

head(model.table[,1:7],10)

updated<-list()
offenders.out<-list(model.name=NA,offenders=NA)

for (i in 1:length(model.table)){

model.thisrun<-model.table[,i]
orig.model.thisrun<-model.table[,i]

# remove NA:
if(length(which(is.na(model.thisrun)))>0) model.thisrun<-as.character(model.thisrun[-which(is.na(model.thisrun))])

# get interaction identifier:
if(length(grep(":",model.thisrun))>0) interactions.pres<-"yes" else interactions.pres<-"no"

# remove interactions:
if(length(grep(":",model.thisrun))>0) model.thisrun<-model.thisrun[-grep(":",model.thisrun)]

# get variables for VIF analysis:
vars.thisrun<-dataset[,colnames(dataset) %in% model.thisrun]
head(vars.thisrun)

# make lifespan numeric:
if(length(which(colnames(vars.thisrun)=="lifespan"))>0){
vars.thisrun$lifespan<-as.character(vars.thisrun$lifespan)
vars.thisrun$lifespan[vars.thisrun$lifespan=="P"]<-0
vars.thisrun$lifespan[vars.thisrun$lifespan=="A"]<-1
vars.thisrun$lifespan<-as.numeric(vars.thisrun$lifespan)
} # close if lifespan present

# do VIF analysis and order large to small:
mcl.thisrun<-mcl("vars.thisrun")
mcl.thisrun<-mcl.thisrun[order(-mcl.thisrun$VIF),]

# if VIF in the first term is > 3, remove it and re-run the VIF analysis until there are no terms with VIF > 3:

offenders<-list()

while(mcl.thisrun$VIF[1]>3){

offending.term<-as.character(mcl.thisrun$response[1])

print(offending.term)

model.thisrun<-model.thisrun[-which(model.thisrun==offending.term)]
vars.thisrun<-dataset[,colnames(dataset) %in% model.thisrun]
head(vars.thisrun)

# make lifespan numeric:
if(length(which(colnames(vars.thisrun)=="lifespan"))>0){
vars.thisrun$lifespan<-as.character(vars.thisrun$lifespan)
vars.thisrun$lifespan[vars.thisrun$lifespan=="P"]<-0
vars.thisrun$lifespan[vars.thisrun$lifespan=="A"]<-1
vars.thisrun$lifespan<-as.numeric(vars.thisrun$lifespan)
} # close if lifespan present

mcl.thisrun<-mcl("vars.thisrun")
mcl.thisrun<-mcl.thisrun[order(-mcl.thisrun$VIF),]

offenders[[length(offenders)+1]]<-offending.term

} # close while VIF>3

offenders.vec<-as.character(unlist(offenders))

# assign the offenders and their corresponding model to a list:
if(length(offenders.vec)>0){
offenders.out$model.name[[i]]<-names(model.table)[i]
offenders.out$offenders[[i]]<-offenders.vec
} # close assign to list

# if any variables have been removed, update the model table:
if(length(offenders.vec)>0){
print("there are offenders")

updated.model.thisrun<-orig.model.thisrun
updated.model.thisrun[grep(offenders.vec, orig.model.thisrun)]<-NA
updated[[i]]<-as.character(updated.model.thisrun)

} # close if offenders present

# if no variables have been removed, keep the model the same:
if(length(offenders.vec)==0){
print("no offenders")
updated.model.thisrun<-orig.model.thisrun
updated[[i]]<-as.character(updated.model.thisrun)
} # close if no offenders

} # close i loop

# ***** WARNING: these need to be run separately for each model table after the loop if any offenders were identified:

# col.mods1<-data.frame(do.call(cbind,updated))
# colnames(col.mods1)<-colnames(col.mods1.allterms)
# col.removed<-data.frame(dataset="colonisation prob",model=offenders.out$model.name,term=offenders.out$offenders)

# cab.mods1<-data.frame(do.call(cbind,updated))
# colnames(cab.mods1)<-colnames(cab.mods1.allterms)
# cab.removed<-data.frame(dataset="colonisation abund",model=offenders.out$model.name,term=offenders.out$offenders)

# cab.mods2<-data.frame(do.call(cbind,updated))
# colnames(cab.mods2)<-colnames(cab.mods2.allterms)
# cab.removed<-data.frame(dataset="colonisation abund",model=offenders.out$model.name,term=offenders.out$offenders)

# removed.terms<-data.frame(rbind(col.removed,spr.removed,cab.removed,imp.removed))

# write.table(col.removed,"rm.txt",sep="\t",row.names=F, quote=F)

# original model tables:
head(col.mods1.allterms[,1:5],10)
head(cab.mods1.allterms[,1:5],10)
head(cab.mods2.allterms[,1:5],10)

# updated model tables with VIF terms removed:
head(col.mods1[1:5],10)
head(cab.mods1[1:5],10)
head(cab.mods2[1:5],10)

# When there are no offending terms removed, just re-name the model table:

# _____ STEP 1:

spr.mods1<-spr.mods1.allterms
sab.mods1<-sab.mods1.allterms

# _____ STEP 2:

col.mods2<-col.mods2.allterms
spr.mods2<-spr.mods2.allterms
sab.mods2<-sab.mods2.allterms


} # close update model tables

## ---- COLONISATION ANALYSIS
{

# Modelled data:
head(coloniseSc,2)
# Unscaled data for plotting:
head(col.unscaled,2)

## ---- STEP 1: Determine the best trait type:
{

# Model table imported and adjusted above:
head(col.mods1[1:10],10)
head(col.mods1[,1:5],10)

# Create basic model META DATA:
mod2_meta<-data.frame(model=colnames(col.mods1))
mod2_meta$trait.type<-sapply(mod2_meta$model,function(x) substr(x,1,gregexpr("_",x)[[1]][1]-1))

# Add the FORMULAE to a separate data frame (they are long so make the df hard to read):
mod2_data<-mod2_meta
mod2_data$formula<-NA
for (i in 1:length(mod2_data[,1])){

if (length(which(is.na(col.mods1[,i])))>0) mod2_data$formula[i]<-paste("pr ~ ",paste(col.mods1[,i][-which(is.na(col.mods1[,i]))],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")

if (length(which(is.na(col.mods1[,i])))==0) mod2_data$formula[i]<-paste("pr ~ ",paste(col.mods1[,i],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")
}

# RUN MODELS in a loop and assign to model name in mod2_meta:
mod2_meta$BIC<-NA
mod2_meta$r2m<-NA
mod2_meta$r2c<-NA
mod2_data
for (i in 1:length(mod2_meta[,1])){

print(paste("beggining model ",mod2_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(mod2_meta$model[i],"_mod",sep="")
formula.thisrun<-mod2_data$formula[i]

mod.thisrun<-glmer(formula=formula.thisrun,data=coloniseSc,family=binomial)

mod2_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(coloniseSc$plsub))))

assign(modname.thisrun,mod.thisrun)

r2.thisrun<-rsquared(mod.thisrun)

mod2_meta$r2m[i]<-r2.thisrun$Marginal
mod2_meta$r2c[i]<-r2.thisrun$Conditional

}

# Use model names in mod2_meta$model, plus "_mod" to get model outputs:
mod2_meta
summary(RelFD_5_mod)

# Add a null model for comparison:
null_form<-paste("pr ~ 1"," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")
null_mod<-glmer(formula=null_form,data=coloniseSc,family=binomial)
summary(null_mod)

null_BIC<-AIC(null_mod,k=log(length(levels(coloniseSc$plsub))))

# Summarise results table:
mod2_res<-mod2_meta[,c("model","trait.type","BIC","r2m","r2c")]
mod2_res$model<-as.character(mod2_res$model)
# UPDATE trait type for three.way:
mod2_res<-rbind(mod2_res,c("Null","NA",null_BIC,"NA","NA"))
mod2_res$BIC<-as.numeric(mod2_res$BIC)
mod2_res<-mod2_res[order(mod2_res$BIC),]
mod2_res<-tidy.df(mod2_res)
mod2_res$rank<-row.names(mod2_res)
mod2_res

# write.table(mod2_res,file="mod.results1.txt",sep="\t",quote=F,row.names=F)

} # close STEP 1

## --- STEP 1 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in mod2_res$model, plus "_mod" (e.g. RelFD_5_mod)
mod2_res

# Extract ALL coefficient tables:
for (i in 1:length(mod2_res[,1])){
coefname.thisrun<-paste(mod2_res$model[i],"_coef",sep="")
modname.thisrun<-paste(mod2_res$model[i],"_mod",sep="")
if(modname.thisrun=="Null_mod") next
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in mod2_res$model, plus "_coef" (e.g. RelFD_5_coef). 
get(paste(mod2_res$model[2],"_coef",sep=""))
summary(get(top.step1))

# Extract TOP MODEL coefficient tables:
{
# remove three way:
tts<-unique(mod2_res$trait.type)
tts<-tts[-which(tts=="NA")]
tm.names<-mod2_res$model[mod2_res$BIC %in% tapply(mod2_res$BIC,mod2_res$trait.type,min)]
tm.names<-tm.names[-which(tm.names=="Null")]
old.names<-c(paste(tm.names,"_coef",sep=""))
new.names<-paste(old.names,"_top",sep="")
coef.list<-lapply(old.names,get)
for (i in 1:length(new.names)){
assign(new.names[i], coef.list[[i]][,c("term","est","se","P","lci.link","uci.link")],)
}

new.names

# MANUAL UPDATE:
RelFD_6_coef_top$model<-"RelFD_6"
AbsFD_8_coef_top$model<-"AbsFD_8"
Gap_6_coef_top$model<-"Gap_6"
Raw_2_coef_top$model<-"Raw_2"
Comm_2_coef_top$model<-"Comm_2"

coefs.top<-rbind(RelFD_6_coef_top,AbsFD_8_coef_top,Gap_6_coef_top,Raw_2_coef_top,Comm_2_coef_top)
head(coefs.top)
tail(coefs.top)
coefs.top$signif<-sign(coefs.top$lci.link)==sign(coefs.top$uci.link)
coefs.top$signif[which(coefs.top$signif==T)]<-"yes"
coefs.top$signif[which(coefs.top$signif==F)]<-"no"
# write.table(coefs.top, "coefstop.txt",sep="\t",row.names=F,quote=F)

} # close top mod coefs

# Plot EFFECT SIZES for all models:
for (i in 1:length(mod2_res[,1])){
coef.thisrun<-paste(mod2_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
if(coef.thisrun=="Null_coef") next
heading.thisrun<-paste("Model = ",mod2_res$model[i],", Rank = ",mod2_res$rank[i],"/",length(mod2_res[,1]),", BIC = ", round(mod2_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 1 coefficients

## ---- STEP 2: Determine the best continuous trait type based on results from STEP 1:
{

# Model table imported and adjusted above:
head(col.mods2,8)

# Create basic model META DATA:
step3_meta<-data.frame(model=colnames(col.mods2))
step3_meta$model<-sapply(as.character(step3_meta$model),function(x) substr(x,gregexpr("_",x)[[1]][1]+1,nchar(x)))

# Add the FORMULAE to a separate data frame:
step3_data<-step3_meta
step3_data$formula<-NA
for (i in 1:length(step3_meta[,1])){
if (length(which(is.na(col.mods2[,i])))>0) step3_data$formula[i]<-paste("pr ~ ",paste(col.mods2[,i][-which(is.na(col.mods2[,i]))],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")

if (length(which(is.na(col.mods2[,i])))==0) step3_data$formula[i]<-paste("pr ~ ",paste(col.mods2[,i],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")
}

# RUN MODELS in a loop and assign to model name in step3_meta:
step3_meta$BIC<-NA
step3_meta$r2m<-NA
step3_meta$r2c<-NA
for (i in 1:length(step3_meta[,1])){

print(paste("beggining model ",step3_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(step3_meta$model[i],"_mod",sep="")
formula.thisrun<-step3_data$formula[i]

mod.thisrun<-glmer(formula=formula.thisrun,data=coloniseSc,family=binomial)

step3_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(coloniseSc$plsub))))

assign(modname.thisrun,mod.thisrun)

r2.thisrun<-rsquared(mod.thisrun)

step3_meta$r2m[i]<-r2.thisrun$Marginal
step3_meta$r2c[i]<-r2.thisrun$Conditional

}

# Summarise results table:
step3_res<-step3_meta
step3_res$BIC<-as.numeric(step3_res$BIC)
step3_res<-step3_res[order(step3_res$BIC),]
step3_res<-tidy.df(step3_res)
step3_res$rank<-1:length(step3_res[,1])
step3_res

# write.table(step3_res,file="step3_res.txt",sep="\t",quote=F,row.names=F)

} # close STEP 2

## --- STEP 2 CONVERGENCE:
{
step3_res

top.step2<-get(paste(step3_res$model[step3_res$rank==1],"_mod",sep=""))
summary(top.step2)

# Re-start the top model (50,000 and 300,000 give the same convergence on re-start, so just do 50,000):
top2.ss<-getME(top.step2,c("theta","fixef"))
top2.restart50<-update(top.step2,start=top2.ss,control=glmerControl(optCtrl=list(maxfun=50000)))
# top2.restart300<-update(top.step2,start=top2.ss,control=glmerControl(optCtrl=list(maxfun=300000)))

# relative convergence criterion:
relgrad.top2.restart50<-with(top2.restart50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.top2.restart50))

relgrad.top2.orig<-with(top.step2@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.top2.orig))

# relgrad.top2.restart300<-with(top2.restart300@optinfo$derivs,solve(Hessian,gradient))
# max(abs(relgrad.top2.restart300))

logLik(top.step2)
logLik(top2.restart50)
AIC(top.step2)
AIC(top2.restart50)

# Converged model coefficient table:
assign(paste("top2.restart50","_coef",sep=""),coef.sum(top2.restart50, est.type="logit"))
top2.restart50_coef

# Plot effect sizes of converged model to ensure the coefficients are the same as the original model
quartz(file=paste("Re_start_top",".pdf",sep=""),type="pdf")
efs.plot(top2.restart50_coef,"Re-started (converged model)")
dev.off()

} # close step 2 convergence

## --- STEP 2 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in step3_res$model, plus "_mod" (e.g. comb_6_mod). Add rank to model results table:
step3_res

# Extract ALL coefficient tables:
for (i in 1:length(step3_res[,1])){
coefname.thisrun<-paste(step3_res$model[i],"_coef",sep="")
modname.thisrun<-paste(step3_res$model[i],"_mod",sep="")
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Extract TOP MODEL coefficient table:
step2_top_coef<-get(paste(step3_res$model[step3_res$rank==1],"_coef",sep=""))

# write.table(step2_top_coef, "step2_top_coef.txt",sep="\t",row.names=F,quote=F)

# Plot EFFECT SIZES for all models:
for (i in 1:length(step3_res[,1])){
coef.thisrun<-paste(step3_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
heading.thisrun<-paste("Model = ",step3_res$model[i],", Rank = ",step3_res$rank[i],"/",length(step3_res[,1]),", BIC = ", round(step3_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 2 coefficients

## --- STEP 2 TOP MODEL ESTIMATES:
{
# Modelled data:
head(coloniseSc,2)
# Unscaled data for plotting:
head(col.unscaled,2)

# TOP MODEL coefficient table:
step2_top_coef
step3_res

# The model to plot:
summary(top.step2)

### ---- NEW DATA FOR ESTIMATES:
{

head(step2_top_coef)

step2_top_coef$signif<-sign(step2_top_coef$lci.link)==sign(step2_top_coef$uci.link)

# Reduce to significant terms:
nd.s3<-data.frame(term=step2_top_coef$term[step2_top_coef$signif==T])
# Remove intercept
nd.s3<-data.frame(term=nd.s3[-which(nd.s3$term=="(Intercept)"),])
# Add an interaction identifier:
nd.s3$interaction<-NA
nd.s3$interaction[grep(":",nd.s3$term)]<-"int"
nd.s3$interaction[-grep(":",nd.s3$term)]<-"main"
# Remove main effects that appear as a significant interaction:
nd.s3<-nd.s3[-which(nd.s3$term  %in% nd.s3$term[nd.s3$interaction=="main"][nd.s3$term[nd.s3$interaction=="main"] %in% unlist(strsplit(as.character(nd.s3$term[nd.s3$interaction=="int"]),":"))]),]
nd.s3<-tidy.df(nd.s3)
head(nd.s3)

# Add an x.axis term column for the estimates.nd function:
nd.s3$xaxis.term<-c("NO3soil","sla_relFD","hgt_absFD","prop_peren","prop_grass","prop_grass","prop_legume","yr","yr","yr","yr","yr")

# Add an interaction term column that will be used by the estimates.nd function:
nd.s3$int.term<-c(rep(NA,3),"lifespan","grass","legume","legume","seed_relFD","ldmc_cwm","lifespan","grass","legume")

# For each term, add a name to be used as the data frame name:
nd.s3$df.name<-c("no3.top2","sla.top2","hgt.top2","perenlfsp.top2","propgrassGR.top2","propgrassLEG.top2","proplegLEG.top2","yrseed.top2","yrldmc.top2","yrlife.top2","yrgrass.top2","yrlegume.top2")

# Get estimates with estimates.nd function for each data frame:
for (i in 1:length(nd.s3[,1])){

if (nd.s3$interaction[i]=="main"){
assign(nd.s3$df.name[i],estimates.nd("top.step2",xaxis.term=nd.s3$xaxis.term[i],modelled.data="coloniseSc",unscaled.data="col.unscaled"))
} # close if main

if (nd.s3$interaction[i]=="int"){
assign(nd.s3$df.name[i],estimates.nd("top.step2", xaxis.term=nd.s3$xaxis.term[i],int.term=nd.s3$int.term[i],modelled.data="coloniseSc",unscaled.data="col.unscaled"))
} # close if interaction

} # close i for data frame loop

} # close new data 

} # close top model estimates

## --- STEP 2 TOP MODEL PLOTS:
{

# Names for the estimate data frames are in nd.s3$df.name:
nd.s3

# Convert nd.s3 into a plot table:
step2.plot<-data.frame(data=nd.s3$df.name,x.label=c("NO3 soil","Relative FD SLA","Absolute FD height","Proportion perennial","Proportion Grass","Proportion Grass","Proportion Legume","Year","Year","Year","Year","Year"),x.variable=paste(nd.s3$xaxis.term,"_unsc",sep=""),interaction=nd.s3$interaction, x.subset=nd.s3$int.term,x1.label=c(rep(NA,3),"P","Grass=0","Legume=0","Legume=0","seed RelFD min","LDMC CWM min","P","Grass=0","Legume=0"),x2.label=c(rep(NA,3),"A","Grass=1","Legume=1","Legume=1","seed RelFD max","LDMC CWM max","A","Grass=1","Legume=1"))

# TOP MODEL plots (moved to separate sections below):

} # close top model plots

} # close colonisation analysis

## ---- SPREAD ANALYSIS
{
head(spreadSc,3)

# Modelled data:
head(spreadSc,2)
# Unscaled data for plotting:
head(spr.unscaled,2)

## ---- STEP 1: Determine the best trait type:
{

# Model table imported and adjusted above:
head(spr.mods1[,1:5],10)
head(spr.mods1.allterms[,1:5],10)

# Create basic model META DATA:
spr1_meta<-data.frame(model=colnames(spr.mods1))
spr1_meta$trait.type<-sapply(spr1_meta$model,function(x) substr(x,1,gregexpr("_",x)[[1]][1]-1))

# Add the FORMULAE to a separate data frame (they are long so make the df hard to read):
spr1_data<-spr1_meta
spr1_data$formula<-NA
for (i in 1:length(spr1_data[,1])){

if (length(which(is.na(spr.mods1[,i])))>0) spr1_data$formula[i]<-paste("spread ~ ",paste(spr.mods1[,i][-which(is.na(spr.mods1[,i]))],collapse=" + ")," + (1 | sp) + (1  | plot )",sep="")

if (length(which(is.na(spr.mods1[,i])))==0) spr1_data$formula[i]<-paste("spread ~ ",paste(spr.mods1[,i],collapse=" + ")," + (1 | sp) + (1  | plot)",sep="")
}

# RUN MODELS in a loop and assign to model name in spr1_meta:
spr1_meta$BIC<-NA
spr1_meta$r2m<-NA
spr1_meta$r2c<-NA
spr1_data
for (i in 1:length(spr1_meta[,1])){

print(paste("beggining model ",spr1_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(spr1_meta$model[i],"_mod",sep="")
formula.thisrun<-spr1_data$formula[i]

mod.thisrun<-glmer(formula=formula.thisrun,data=spreadSc,family=binomial)

spr1_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(spreadSc$plot))))

assign(modname.thisrun,mod.thisrun)

r2.thisrun<-rsquared(mod.thisrun)

spr1_meta$r2m[i]<-r2.thisrun$Marginal
spr1_meta$r2c[i]<-r2.thisrun$Conditional

}
save.image("e93_wksp")

# Use model names in spr1_meta$model, plus "_mod" to get model outputs:
spr1_meta

# Add a null model for comparison:
null_form_spr<-paste("spread ~ 1"," + (1 | sp) + (1  | plot)",sep="")
null_mod_spr<-glmer(formula=null_form_spr,data=spreadSc,family=binomial)
summary(null_mod_spr)

null_BIC_spr<-AIC(null_mod_spr,k=log(length(levels(spreadSc$plot))))

# Summarise results table:
spr1_res<-spr1_meta[,c("model","trait.type","BIC","r2m","r2c")]
spr1_res$model<-as.character(spr1_res$model)
spr1_res<-rbind(spr1_res,c("Null","NA",null_BIC_spr,"NA","NA"))
spr1_res$BIC<-as.numeric(spr1_res$BIC)
spr1_res<-spr1_res[order(spr1_res$BIC),]
spr1_res<-tidy.df(spr1_res)
spr1_res$rank<-row.names(spr1_res)
spr1_res

top.spread1<-get(paste(spr1_res$model[spr1_res$rank==1],"_mod",sep=""))
summary(top.spread1)

# write.table(spr1_res,file="spr1_result.txt",sep="\t",quote=F,row.names=F)

} # close STEP 1

## --- STEP 1 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in spr1_res$model, plus "_mod" (e.g. RelFD_2_mod)
spr1_res

# Extract ALL coefficient tables:
for (i in 1:length(spr1_res[,1])){
coefname.thisrun<-paste(spr1_res$model[i],"_coef",sep="")
modname.thisrun<-paste(spr1_res$model[i],"_mod",sep="")
if(modname.thisrun=="Null_mod") next
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in mod2_res$model, plus "_coef" (e.g. RelFD_5_coef). 
get(paste(spr1_res$model[1],"_coef",sep=""))
summary(top.spread1)

# Extract TOP MODEL coefficient tables:
{
tts<-unique(spr1_res$trait.type)
tts<-tts[-which(tts=="NA")]
tm.names<-spr1_res$model[spr1_res$BIC %in% tapply(spr1_res$BIC,spr1_res$trait.type,min)]
tm.names<-tm.names[-which(tm.names=="Null")]
old.names<-c(paste(tm.names,"_coef",sep=""))
new.names<-paste(old.names,"_top",sep="")
coef.list<-lapply(old.names,get)
for (i in 1:length(new.names)){
assign(new.names[i], coef.list[[i]][,c("term","est","se","P","lci.link","uci.link")],)
}

new.names

# MANUAL UPDATE:
RelFD_2_coef_top$model<-"RelFD_2"
AbsFD_2_coef_top$model<-"AbsFD_2"
Gap_1_coef_top$model<-"Gap_1"
Raw_2_coef_top$model<-"Raw_2"
Comm_1_coef_top$model<-"Comm_1"

coefs.top<-rbind(RelFD_2_coef_top,AbsFD_2_coef_top,Gap_1_coef_top,Raw_2_coef_top,Comm_1_coef_top)
head(coefs.top)
tail(coefs.top)
coefs.top$signif<-sign(coefs.top$lci.link)==sign(coefs.top$uci.link)
coefs.top$signif[which(coefs.top$signif==T)]<-"yes"
coefs.top$signif[which(coefs.top$signif==F)]<-"no"
# write.table(coefs.top, "coefstop.txt",sep="\t",row.names=F,quote=F)

} # close top mod coefs

# Plot EFFECT SIZES for all models:
for (i in 1:length(spr1_res[,1])){
coef.thisrun<-paste(spr1_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
if(coef.thisrun=="Null_coef") next
heading.thisrun<-paste("Model = ",spr1_res$model[i],", Rank = ",spr1_res$rank[i],"/",length(spr1_res[,1]),", BIC = ", round(spr1_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 1 coefficients

## --- STEP 1 CONVERGENCE:
{
spr1_res
summary(top.spread1)

# Re-start the top model:
topspr1.ss<-getME(top.spread1,c("theta","fixef"))
topspr1.restart50<-update(top.spread1,start=topspr1.ss,control=glmerControl(optCtrl=list(maxfun=50000)))

# relative convergence criterion:
relgrad.topspr1.restart50<-with(topspr1.restart50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.topspr1.restart50))
format(max(abs(relgrad.topspr1.restart50)),scientific=FALSE)

logLik(top.spread1)
logLik(topspr1.restart50)
AIC(top.spread1)
AIC(topspr1.restart50)

# Converged model coefficient table:
assign(paste("topspr1.restart50","_coef",sep=""),coef.sum(topspr1.restart50, est.type="logit"))
topspr1.restart50_coef

# Plot effect sizes of converged model to ensure the coefficients are the same as the original model
efs.plot(topspr1.restart50_coef,"Re-started (converged model)")

# Re-start best AbsFD model:

summary(AbsFD_2_mod)$coefficients
absspr1.ss<-getME(AbsFD_2_mod,c("theta","fixef"))
absspr1.restart50<-update(AbsFD_2_mod,start=absspr1.ss,control=glmerControl(optCtrl=list(maxfun=1000000)))
summary(absspr1.restart50)$coefficients

# This model has sla_absFD and sla_absFD:yr removed to make it converge:
ab2r<-glmer(spread ~ yr + sdsub_abund +  hgt_absFD + seed_absFD +   ldmc_absFD + lifespan + grass + legume + hgt_absFD:yr +   seed_absFD:yr + ldmc_absFD:yr + lifespan:yr + grass:yr +  legume:yr + (1 | sp) + (1 | plot), data=spreadSc, family=binomial,glmerControl(optCtrl=list(maxfun=1000000)))
summary(ab2r)
summary(AbsFD_1_mod)

# This other simplified model (no hgt_absFD and ldmc_absFD and their interactions) shows that sla_absFD and sla_absFD:yr are not significant:
ab2r2<-glmer(spread ~ yr + sdsub_abund +  sla_absFD + seed_absFD + lifespan + grass + legume + sla_absFD:yr + seed_absFD:yr + lifespan:yr + grass:yr +  legume:yr + (1 | sp) + (1 | plot), data=spreadSc, family=binomial,glmerControl(optCtrl=list(maxfun=1000000)))
summary(ab2r2)

logLik(AbsFD_2_mod)
logLik(absspr1.restart50)
AIC(AbsFD_2_mod)
AIC(absspr1.restart50)

# relative convergence criterion:
relgrad.absspr1.restart50<-with(absspr1.restart50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.absspr1.restart50))
format(max(abs(relgrad.absspr1.restart50)),scientific=FALSE)

coef.absspr1.restart<-coef.sum(absspr1.restart50, est.type="logit")
coef.ab2r<-coef.sum(ab2r,est.type="logit")
coef.ab2r

coef.ab2r<-coef.sum(ab2r,est.type="logit")
coef.ab2r

coef.ab2r2<-coef.sum(ab2r2,est.type="logit")
coef.ab2r2

summary(AbsFD_1_mod)
coef.ab1<-coef.sum(AbsFD_1_mod,est.type="logit")
coef.ab1

efs.plot(coef.absspr1.restart,"Best AbsFD model (AbsFD_2) restarted")

coef.ab1$signif<-sign(coef.ab1$lci.link)==sign(coef.ab1$uci.link)
coef.ab1$signif[which(coef.ab1$signif==T)]<-"yes"
coef.ab1$signif[which(coef.ab1$signif==F)]<-"no"
# write.table(coef.ab1, "coefstop.txt",sep="\t",row.names=F,quote=F)

} # close step 2 convergence

## ---- STEP 2: Determine the best continuous trait type based on results from STEP 1:
{

# Model table imported and adjusted above:
head(spr.mods2[,1:4],10)
head(spr.mods2.allterms[,1:4],10)

# Create basic model META DATA:
spr2_meta<-data.frame(model=colnames(spr.mods2))

# Add the FORMULAE to a separate data frame:
spr2_data<-spr2_meta
spr2_data$formula<-NA
for (i in 1:length(spr2_meta[,1])){

if (length(which(is.na(spr.mods2[,i])))>0) spr2_data$formula[i]<-paste("spread ~ ",paste(spr.mods2[,i][-which(is.na(spr.mods2[,i]))],collapse=" + ")," + (1 | sp) + (1  | plot )",sep="")

if (length(which(is.na(spr.mods2[,i])))==0) spr2_data$formula[i]<-paste("spread ~ ",paste(spr.mods2[,i],collapse=" + ")," + (1 | sp) + (1  | plot)",sep="")

}

# RUN MODELS in a loop and assign to model name in step3_meta:
spr2_meta$BIC<-NA
spr2_meta$r2m<-NA
spr2_meta$r2c<-NA
spr2_data
for (i in 1:length(spr2_meta[,1])){

print(paste("beggining model ",spr2_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(spr2_meta$model[i],"_mod",sep="")
formula.thisrun<-spr2_data$formula[i]

mod.thisrun<-glmer(formula=formula.thisrun,data=spreadSc,family=binomial)

spr2_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(spreadSc$plot))))

assign(modname.thisrun,mod.thisrun)

r2.thisrun<-rsquared(mod.thisrun)

spr2_meta$r2m[i]<-r2.thisrun$Marginal
spr2_meta$r2c[i]<-r2.thisrun$Conditional

}	
save.image("e93_wksp")

# Summarise results table:
spr2_res<-spr2_meta
spr2_res$BIC<-as.numeric(spr2_res$BIC)
spr2_res<-spr2_res[order(spr2_res$BIC),]
spr2_res<-tidy.df(spr2_res)
spr2_res$rank<-1:length(spr2_res[,1])
spr2_res

# write.table(spr2_res,file="spr2_res.txt",sep="\t",quote=F,row.names=F)

top.spread2<-get(paste(spr2_res$model[spr2_res$rank==1],"_mod",sep=""))
summary(top.spread2)

} # close STEP 2

## --- STEP 2 CONVERGENCE:
{
spr2_res

summary(top.spread2)

# Re-start the top model:
topspr2.ss<-getME(top.spread2,c("theta","fixef"))
topspr2.restart50<-update(top.spread2,start=topspr2.ss,control=glmerControl(optCtrl=list(maxfun=50000)))

# relative convergence criterion:
relgrad.topspr2.restart50<-with(topspr2.restart50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.topspr2.restart50))
format(max(abs(relgrad.topspr2.restart50)),scientific=FALSE)

# original model:
relgrad.topspr2<-with(top.spread2@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.topspr2))
format(max(abs(relgrad.topspr2)),scientific=FALSE)

logLik(top.spread2)
logLik(topspr2.restart50)
AIC(top.spread2)
AIC(topspr2.restart50)

# Converged model coefficient table:
assign(paste("topspr2.restart50","_coef",sep=""),coef.sum(topspr2.restart50, est.type="logit"))
topspr2.restart50_coef

# Plot effect sizes of converged model to ensure the coefficients are the same as the original model
quartz(file=paste("Re_start_top",".pdf",sep=""),type="pdf")
efs.plot(topspr2.restart50_coef,"Re-started (converged model)")
dev.off()

} # close step 2 convergence

## --- STEP 2 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in spr2_res$model, plus "_mod" (e.g. Mdl2_1_mod)
spr2_res

# Extract ALL coefficient tables:
for (i in 1:length(spr2_res[,1])){
coefname.thisrun<-paste(spr2_res$model[i],"_coef",sep="")
modname.thisrun<-paste(spr2_res$model[i],"_mod",sep="")
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in spr2_res$model, plus "_coef" (e.g. Mdl2_1_coef). 
get(paste(spr2_res$model[1],"_coef",sep=""))
summary(top.spread2)

# Extract TOP MODEL coefficient table:
step2_topspr_coef<-get(paste(spr2_res$model[spr2_res$rank==1],"_coef",sep=""))

# write.table(step2_topspr_coef, "step2_topspr_coef.txt",sep="\t",row.names=F,quote=F)

# Plot EFFECT SIZES for all models:
for (i in 1:length(spr2_res[,1])){
coef.thisrun<-paste(spr2_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
heading.thisrun<-paste("Model = ",spr2_res$model[i],", Rank = ",spr2_res$rank[i],"/",length(spr2_res[,1]),", BIC = ", round(spr2_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 2 coefficients

## --- STEP 2 TOP MODEL ESTIMATES:
{
# Modelled data:
head(spreadSc,2)
# Unscaled data for plotting:
head(spr.unscaled,2)

# TOP MODEL coefficient table:
step2_topspr_coef
spr2_res

# The model to plot:
summary(top.spread2)

} # close estimates

### ---- STEP 2 TOP MODEL PLOTS:
{

# NEW DATA ESTIMATES
{

step2b_topspr_coef<-step2_topspr_coef
head(step2b_topspr_coef)

step2b_topspr_coef$signif<-FALSE
step2b_topspr_coef$signif[step2b_topspr_coef$P<0.06]<-TRUE
step2b_topspr_coef

# Reduce to significant terms:
nd.spr2b<-data.frame(term=step2b_topspr_coef$term[step2b_topspr_coef$signif==T])
# Remove intercept
nd.spr2b<-data.frame(term=nd.spr2b[-which(nd.spr2b$term=="(Intercept)"),])
# Add an interaction identifier:
nd.spr2b$interaction<-NA
nd.spr2b$interaction[grep(":",nd.spr2b$term)]<-"int"
nd.spr2b$interaction[-grep(":",nd.spr2b$term)]<-"main"
nd.spr2b<-tidy.df(nd.spr2b)
nd.spr2b

# Add an x.axis term column for the estimates.nd function:
nd.spr2b$xaxis.term<-c("sdsub_abund","hgt_absFD","grass","yr","yr","yr","yr")

# Add an interaction term column that will be used by the estimates.nd function:
nd.spr2b$int.term<-c(rep(NA,3),"sla_relFD","seed_relFD","ldmc_relFD","lifespan")

# For each term, add a name to be used as the data frame name:
nd.spr2b$df.name<-c("sdabund.topspr2","hgt.topspr2","grass.topspr2","yrsla.topspr2","yrseed.topspr2","yrldmc.topspr2","lifesp.topspr2")

# Get estimates with estimates.nd function for each data frame:
for (i in 1:length(nd.spr2b[,1])){

if (nd.spr2b$interaction[i]=="main"){
assign(nd.spr2b$df.name[i],estimates.nd("top.spread2",xaxis.term=nd.spr2b$xaxis.term[i],modelled.data="spreadSc",unscaled.data="spr.unscaled"))
} # close if main

if (nd.spr2b$interaction[i]=="int"){
assign(nd.spr2b$df.name[i],estimates.nd("top.spread2", xaxis.term=nd.spr2b$xaxis.term[i],int.term=nd.spr2b$int.term[i],modelled.data="spreadSc",unscaled.data="spr.unscaled"))
} # close if interaction

} # close for data frame loop

} # close new data

## --- STEP 2 TOP MODEL PLOTS:
{

# Names for the estimate data frames are in nd.spr2b$df.name:
nd.spr2b

# Convert nd.spr2b into a plot table:
step2bspr.plot<-data.frame(data=nd.spr2b$df.name,x.label=c("Seeded subplot abundance","Absolute FD height","Grass","Year","Year","Year","Year"),x.variable=paste(nd.spr2b$xaxis.term,"_unsc",sep=""),interaction=nd.spr2b$interaction, x.subset=nd.spr2b$int.term,x1.label=c(rep(NA,2),"Grass=0","SLA RelFD min","seed RelFD min","LDMC RelFD min","Lifespan=P"),x2.label=c(rep(NA,2),"Grass=1","SLA RelFD max","seed RelFD max","LDMC RelFD max","Lifespan=A"))

# TOP MODEL plots (moved to separate section below):

} # close top model plots

} # close plots

} # close spread analysis

## ---- COLONISATION ABUNDANCE ANALYSIS
{
# Modelled data:
head(colabSc,2)
# Unscaled data for plotting:
head(cab.unscaled,2)

## ---- STEP 1: Determine the best trait type:
{

# Model table imported and adjusted above:
head(cab.mods1[,1:5],10)
head(cab.mods1.allterms[,1:5],10)

# Create basic model META DATA:
cab1_meta<-data.frame(model=colnames(cab.mods1))
cab1_meta$trait.type<-sapply(cab1_meta$model,function(x) substr(x,1,gregexpr("_",x)[[1]][1]-1))

# Add the FORMULAE to a separate data frame (they are long so make the df hard to read):
cab1_data<-cab1_meta
cab1_data$formula<-NA
cab1_data
for (i in 1:length(cab1_data[,1])){

if (length(which(is.na(cab.mods1[,i])))>0) cab1_data$formula[i]<-paste("std_cover ~ ",paste(cab.mods1[,i][-which(is.na(cab.mods1[,i]))],collapse=" + ")," + (NO3soil | sp) + (1  | plot / plsub)",sep="")

if (length(which(is.na(cab.mods1[,i])))==0) cab1_data$formula[i]<-paste("std_cover ~ ",paste(cab.mods1[,i],collapse=" + ")," + (NO3soil | sp) + (1  | plot / plsub)",sep="")

}

# RUN MODELS in a loop and assign to model name in cab1_meta:
######## WARNING 4.5 HOURS
cab1_meta$BIC<-NA
# cab1_meta$r2m<-NA
# cab1_meta$r2c<-NA
cab1_data

# The model names (e.g. "Comm_1") have been overwritten by the spread models, so change the names here and run again:

cab1_meta$model<-paste(cab1_meta$model,"_ab",sep="")

for (i in 1:length(cab1_meta[,1])){

print(paste("beggining model ",cab1_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(cab1_meta$model[i],"_mod",sep="")
formula.thisrun<-cab1_data$formula[i]

mod.thisrun<-glmmadmb(formula=as.formula(formula.thisrun),data=colabSc,family="beta")

cab1_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(colabSc$plsub))))

assign(modname.thisrun,mod.thisrun)

}
save.image("e93_wksp")

# Use model names in mod2_meta$model, plus "_mod" to get model outputs:
cab1_meta
summary(RelFD_1_mod)
lmer.pval(summary(RelFD_1_mod)$coefficients[8,3])

# Add a null model for comparison:
null_form_cab<-paste("std_cover ~ 1"," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")
null_mod_cab<-glmmadmb(formula=as.formula(null_form_cab),data=colabSc,family="beta")
summary(null_mod_cab)

null_BIC_cab<-AIC(null_mod_cab,k=log(length(levels(colabSc$plsub))))

# Summarise results table:
cab1_res<-cab1_meta[,c("model","trait.type","BIC")]
cab1_res$model<-as.character(cab1_res$model)
cab1_res<-rbind(cab1_res,c("Null","NA",null_BIC_cab))
cab1_res$BIC<-as.numeric(cab1_res$BIC)
cab1_res<-cab1_res[order(cab1_res$BIC),]
cab1_res<-tidy.df(cab1_res)
cab1_res$rank<-row.names(cab1_res)
cab1_res	

head(colabSc,2)

top.cab1<-get(paste(cab1_res$model[cab1_res$rank==1],"_mod",sep=""))
summary(top.cab1)

# write.table(cab1_res,file="cab1_res.txt",sep="\t",quote=F,row.names=F)

} # close STEP 1

## --- STEP 1 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in cab1_res$model, plus "_mod" (e.g. RelFD_1_mod)
cab1_res

# Extract ALL coefficient tables:
for (i in 1:length(cab1_res[,1])){
coefname.thisrun<-paste(cab1_res$model[i],"_coef",sep="")
modname.thisrun<-paste(cab1_res$model[i],"_mod",sep="")
if(modname.thisrun=="Null_mod") next
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in cab1_res$model, plus "_coef" (e.g. RelFD_1_coef). 
get(paste(cab1_res$model[1],"_coef",sep=""))
summary(get(top.step1))

# Extract TOP MODEL coefficient tables:
{
# remove three way:
cab1_res2<-cab1_res[-which(cab1_res$model %in% "Null"),]
tts<-unique(cab1_res2$trait.type)
tm.names<-cab1_res2$model[cab1_res2$BIC %in% tapply(cab1_res2$BIC,cab1_res2$trait.type,min)]
old.names<-c(paste(tm.names,"_coef",sep=""))
new.names<-paste(old.names,"_top",sep="")
coef.list<-lapply(old.names,get)
for (i in 1:length(new.names)){
assign(new.names[i], coef.list[[i]][,c("term","est","se","P","lci.link","uci.link")],)
}

new.names

# MANUAL UPDATE:
Comm_1_coef_top$model<-"Comm_1"
RelFD_2_coef_top$model<-"RelFD_2"
AbsFD_5_coef_top$model<-"AbsFD_5"
Gap_1_coef_top$model<-"Gap_1"
Raw_sla_seed_2_coef_top$model<-"Raw_sla_seed_2"

coefs.top<-rbind(Comm_1_coef_top,RelFD_2_coef_top,AbsFD_5_coef_top,Gap_1_coef_top,Raw_sla_seed_2_coef_top)
head(coefs.top)
tail(coefs.top)
coefs.top$signif<-sign(coefs.top$lci.link)==sign(coefs.top$uci.link)
coefs.top$signif[which(coefs.top$signif==T)]<-"yes"
coefs.top$signif[which(coefs.top$signif==F)]<-"no"
# write.table(coefs.top, "coefstop.txt",sep="\t",row.names=F,quote=F)

} # close top mod coefs

# Plot EFFECT SIZES for all models:
for (i in 1:length(cab1_res[,1])){
coef.thisrun<-paste(cab1_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
if(coef.thisrun=="Null_coef") next
heading.thisrun<-paste("Model = ",cab1_res$model[i],", Rank = ",cab1_res$rank[i],"/",length(cab1_res[,1]),", BIC = ", round(cab1_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 1 coefficients

## ---- STEP 2: Determine the best continuous trait type based on results from STEP 1:
{

# Model table imported and adjusted above:
head(cab.mods2[,1:4],10)
head(cab.mods2.allterms[,1:4],10)

# Create basic model META DATA:
cab2_meta<-data.frame(model=colnames(cab.mods2))

# Add the FORMULAE to a separate data frame:
cab2_data<-cab2_meta
cab2_data$formula<-NA
for (i in 1:length(cab2_meta[,1])){

if (length(which(is.na(cab.mods2[,i])))>0) cab2_data$formula[i]<-paste("std_cover ~ ",paste(cab.mods2[,i][-which(is.na(cab.mods2[,i]))],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")

if (length(which(is.na(cab.mods2[,i])))==0) cab2_data$formula[i]<-paste("std_cover ~ ",paste(cab.mods2[,i],collapse=" + ")," + (1 + NO3soil | sp) + (1  | plot / plsub)",sep="")

} # close add formulae

# RUN MODELS in a loop and assign to model name in step3_meta:
cab2_meta$BIC<-NA
cab2_data
for (i in 1:length(cab2_meta[,1])){

print(paste("beggining model ",cab2_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(cab2_meta$model[i],"_mod",sep="")
formula.thisrun<-cab2_data$formula[i]

mod.thisrun<-glmmadmb(formula=as.formula(formula.thisrun),data=colabSc,family="beta")

cab2_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(colabSc$plsub))))

assign(modname.thisrun,mod.thisrun)

}
save.image("e93_wksp")

# Summarise results table:
cab2_res<-cab2_meta
cab2_res$BIC<-as.numeric(cab2_res$BIC)
cab2_res<-cab2_res[order(cab2_res$BIC),]
cab2_res<-tidy.df(cab2_res)
cab2_res$rank<-1:length(cab2_res[,1])
cab2_res

# write.table(cab2_res,file="cab2_res.txt",sep="\t",quote=F,row.names=F)

top.cab2<-get(paste(cab2_res$model[cab2_res$rank==1],"_mod",sep=""))
summary(top.cab2)

} # close STEP 2

## --- STEP 2 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in spr2_res$model, plus "_mod" (e.g. Model1_13)
cab2_res

# Extract ALL coefficient tables:
for (i in 1:length(cab2_res[,1])){
coefname.thisrun<-paste(cab2_res$model[i],"_coef",sep="")
modname.thisrun<-paste(cab2_res$model[i],"_mod",sep="")
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in spr2_res$model, plus "_coef" (e.g. Model1_13). 
get(paste(cab2_res$model[1],"_coef",sep=""))
summary(top.cab2)

# Extract TOP MODEL coefficient table:
step2_topcab_coef<-get(paste(cab2_res$model[cab2_res$rank==1],"_coef",sep=""))

# write.table(step2_topcab_coef, "step2_topcab_coef.txt",sep="\t",row.names=F,quote=F)

# Plot EFFECT SIZES for all models:
for (i in 1:length(cab2_res[,1])){
coef.thisrun<-paste(cab2_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
heading.thisrun<-paste("Model = ",cab2_res$model[i],", Rank = ",cab2_res$rank[i],"/",length(cab2_res[,1]),", BIC = ", round(cab2_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 2 coefficients

### ---- STEP 2 TOP MODEL PLOTS:
{

# NEW DATA ESTIMATES
{
step2b_topcab_coef<-step2_topcab_coef
head(step2b_topcab_coef)

step2b_topcab_coef$signif<-FALSE
step2b_topcab_coef$signif[step2b_topcab_coef$P<0.06]<-TRUE
step2b_topcab_coef

# Reduce to significant terms:
nd.cab2b<-data.frame(term=step2b_topcab_coef$term[step2b_topcab_coef$signif==T])
# Remove intercept
nd.cab2b<-data.frame(term=nd.cab2b[-which(nd.cab2b$term=="(Intercept)"),])
# Add an interaction identifier:
nd.cab2b$interaction<-NA
if(length(grep(":",nd.cab2b$term))>0) nd.cab2b$interaction[grep(":",nd.cab2b$term)]<-"int"
if(length(grep(":",nd.cab2b$term))>0) nd.cab2b$interaction[-grep(":",nd.cab2b$term)]<-"main"
if(length(grep(":",nd.cab2b$term))==0) nd.cab2b$interaction<-"main"
nd.cab2b<-tidy.df(nd.cab2b)

# Add an x.axis term column for the estimates.nd function:
nd.cab2b$xaxis.term<-c("yr","NO3soil","prop_peren","prop_legume","hgt_gap","seed_cwm","ldmc_cwm")

# For each term, add a name to be used as the data frame name:
nd.cab2b$df.name<-c("yr.topcab2","no3.topcab2","propperen.topcab2","propleg.topcab2","hgt.topcab2","seed.topcab2","ldmc.topcab2")

# Get estimates with estimates.nd function for each data frame:
for (i in 1:length(nd.cab2b[,1])){

if (nd.cab2b$interaction[i]=="main"){
assign(nd.cab2b$df.name[i],estimates.nd("top.cab2",xaxis.term=nd.cab2b$xaxis.term[i],modelled.data="colabSc",unscaled.data="cab.unscaled"))
} # close if main

} # close for data frame loop

} # close new data 

## --- STEP 2 TOP MODEL PLOTS:
{

# Names for the estimate data frames are in nd.cab2b$df.name:
nd.cab2b

# Convert nd.cab2b into a plot table:
cab2b.plot<-data.frame(data=nd.cab2b$df.name,x.label=c("Year","NO3 soil","Propn perennial","Propn legume","Height gap","Seed community wt mean","LDMC community wt mean"),x.variable=paste(nd.cab2b$xaxis.term,"_unsc",sep=""),interaction=nd.cab2b$interaction,x1.label=c(rep(NA,7)),x2.label=c(rep(NA,7)))

# TOP MODEL plots (moved to separate section below):

} # close top model plots

} # close plots 

} # close col abund analysis

## ---- SPREAD ABUNDANCE ANALYSIS
{
# spread abundance followed the same rules as spread probability, e.g. the species had to have been previously present for it to be coded as spread, etc. (see details in data set-up script)

# Make NAs zero because the data frame includes all species that were seeded in that plot (and not present in 1991) that haven't spread (spread=0). Not sure what this step was for - there are no NAs in the spread data. 
sprabSc<-spreadSc
sab.unscaled<-spr.unscaled
sprabSc$spr_abund_raw[which(is.na(sprabSc$spr_abund_raw))]<-0
sab.unscaled$spr_abund_raw[which(is.na(sab.unscaled$spr_abund_raw))]<-0
sprabSc$std_cover[which(is.na(sprabSc$std_cover))]<-0
sab.unscaled$std_cover[which(is.na(sab.unscaled$std_cover))]<-0

head(sprabSc,3)
head(sab.unscaled,3)

sab.unscaled$plot_yr<-as.factor(paste(sab.unscaled$plot,sab.unscaled$year,sep="_"))
sprabSc$plot_yr<-as.factor(paste(sprabSc$plot,sprabSc$year,sep="_"))

# Remove zeros:
sprabSc<-sprabSc[-which(sprabSc$std_cover==0),]
sab.unscaled<-sab.unscaled[-which(sab.unscaled$std_cover==0),]
sprabSc<-tidy.df(sprabSc)
sab.unscaled<-tidy.df(sab.unscaled)

# Divide std_cover by 2:
sprabSc$std_cover<-sprabSc$std_cover/2
sab.unscaled$std_cover<-sab.unscaled$std_cover/2

# Modelled data:
head(sprabSc,3); dim(sprabSc)

# Unscaled data for plotting:
head(sab.unscaled,3); dim(sab.unscaled)

## ---- STEP 1: Determine the best trait type:
{

# Model table imported and adjusted above:
head(sab.mods1[,1:5],10)
head(sab.mods1.allterms[,1:5],10)

# Create basic model META DATA:
sab1_meta<-data.frame(model=colnames(sab.mods1))
sab1_meta$trait.type<-sapply(sab1_meta$model,function(x) substr(x,1,gregexpr("_",x)[[1]][1]-1))

# Add the FORMULAE to a separate data frame (they are long so make the df hard to read):
sab1_data<-sab1_meta
sab1_data$formula<-NA
for (i in 1:length(sab1_data[,1])){

if (length(which(is.na(sab.mods1[,i])))>0) sab1_data$formula[i]<-paste("std_cover ~ ",paste(sab.mods1[,i][-which(is.na(sab.mods1[,i]))],collapse=" + ")," + (1 | sp) + (1  | plot )",sep="")

if (length(which(is.na(sab.mods1[,i])))==0) sab1_data$formula[i]<-paste("std_cover ~ ",paste(sab.mods1[,i],collapse=" + ")," + (1 | sp) + (1  | plot)",sep="")
}

# RUN MODELS in a loop and assign to model name in spr1_meta:
sab1_meta$BIC<-NA
sab1_data

for (i in 1:length(sab1_meta[,1])){

print(paste("beggining model ",sab1_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(sab1_meta$model[i],"_mod",sep="")
formula.thisrun<-sab1_data$formula[i]

mod.thisrun<-glmmadmb(formula=as.formula(formula.thisrun),data=sprabSc,family="beta")

sab1_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(sprabSc$plot))))

assign(modname.thisrun,mod.thisrun)

}

save.image("e93_wksp_sprab")

# Use model names in spr1_meta$model, plus "_mod" to get model outputs:
sab1_meta

# Add a null model for comparison:
null_form_sab<-paste("std_cover ~ 1"," + (1 | sp) + (1  | plot)",sep="")
null_mod_sab<-glmmadmb(formula=as.formula(null_form_sab),data=sprabSc,family="beta")
summary(null_mod_sab)

null_BIC_sab<-AIC(null_mod_sab,k=log(length(levels(sprabSc$plot))))

# Summarise results table:
sab1_res<-sab1_meta[,c("model","trait.type","BIC")]
sab1_res$model<-as.character(sab1_res$model)
sab1_res<-rbind(sab1_res,c("Null","NA",null_BIC_sab))
sab1_res$BIC<-as.numeric(sab1_res$BIC)
sab1_res<-sab1_res[order(sab1_res$BIC),]
sab1_res<-tidy.df(sab1_res)
sab1_res$rank<-row.names(sab1_res)
sab1_res

top.sab1<-get(paste(sab1_res$model[sab1_res$rank==1],"_mod",sep=""))
summary(top.sab1)

# write.table(sab1_res,file="sab1_result.txt",sep="\t",quote=F,row.names=F)

} # close STEP 1

## --- STEP 1 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in sab1_res$model, plus "_mod" (e.g. Comm_2_mod)
get(paste(sab1_res$model[1],"_mod",sep=""))

# Extract ALL coefficient tables:
for (i in 1:length(sab1_res[,1])){
coefname.thisrun<-paste(sab1_res$model[i],"_coef",sep="")
modname.thisrun<-paste(sab1_res$model[i],"_mod",sep="")
if(modname.thisrun=="Null_mod") next
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in sab1_res$model, plus "_coef" (e.g. RelFD_1_coef). 
get(paste(sab1_res$model[1],"_coef",sep=""))
summary(top.sab1)

# Extract TOP MODEL coefficient tables:
{
# remove three way:
sab1_res2<-sab1_res[-which(sab1_res$model %in% "Null"),]
tts<-unique(sab1_res2$trait.type)
tm.names<-sab1_res2$model[sab1_res2$BIC %in% tapply(sab1_res2$BIC,sab1_res2$trait.type,min)]
old.names<-c(paste(tm.names,"_coef",sep=""))
new.names<-paste(old.names,"_top",sep="")
coef.list<-lapply(old.names,get)
for (i in 1:length(new.names)){
assign(new.names[i], coef.list[[i]][,c("term","est","se","P","lci.link","uci.link")],)
}

new.names

# MANUAL UPDATE:
Comm_2_coef_top$model<-"Comm_2"
RelFD_1_coef_top$model<-"RelFD_1"
AbsFD_1_coef_top$model<-"AbsFD_1"
Gap_1_coef_top$model<-"Gap_1"
Raw_1_coef_top$model<-"Raw_1"

coefs.top<-rbind(Comm_2_coef_top,RelFD_1_coef_top,AbsFD_1_coef_top,Gap_1_coef_top,Raw_1_coef_top)
head(coefs.top)
tail(coefs.top)
coefs.top$signif<-sign(coefs.top$lci.link)==sign(coefs.top$uci.link)
coefs.top$signif[which(coefs.top$signif==T)]<-"yes"
coefs.top$signif[which(coefs.top$signif==F)]<-"no"
# write.table(coefs.top, "coefstop.txt",sep="\t",row.names=F,quote=F)

} # close top mod coefs

# Plot EFFECT SIZES for all models:
for (i in 1:length(sab1_res[,1])){
coef.thisrun<-paste(sab1_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
if(coef.thisrun=="Null_coef") next
heading.thisrun<-paste("Model = ",sab1_res$model[i],", Rank = ",sab1_res$rank[i],"/",length(sab1_res[,1]),", BIC = ", round(sab1_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 1 coefficients

## ---- STEP 2: Determine the best continuous trait type based on results from STEP 1:
{

# Model table imported and adjusted above:
head(sab.mods2[,1:4],10)
head(sab.mods2.allterms[,1:4],10)

# Create basic model META DATA:
sab2_meta<-data.frame(model=colnames(sab.mods2))

# Add the FORMULAE to a separate data frame:
sab2_data<-sab2_meta
sab2_data$formula<-NA
for (i in 1:length(sab2_meta[,1])){

if (length(which(is.na(sab.mods2[,i])))>0) sab2_data$formula[i]<-paste("std_cover ~ ",paste(sab.mods2[,i][-which(is.na(sab.mods2[,i]))],collapse=" + ")," + (1 | sp) + (1  | plot )",sep="")

if (length(which(is.na(sab.mods2[,i])))==0) sab2_data$formula[i]<-paste("std_cover ~ ",paste(sab.mods2[,i],collapse=" + ")," + (1 | sp) + (1  | plot)",sep="")

} # close add formulae

# RUN MODELS in a loop and assign to model name in step3_meta:
sab2_meta$BIC<-NA
sab2_data
for (i in 1:length(sab2_meta[,1])){

print(paste("beggining model ",sab2_meta$model[i],", i = ",i,", time = ",Sys.time(),sep=""))

modname.thisrun<-paste(sab2_meta$model[i],"_mod",sep="")
formula.thisrun<-sab2_data$formula[i]

mod.thisrun<-glmmadmb(formula=as.formula(formula.thisrun),data=sprabSc,family="beta")

sab2_meta$BIC[i]<-AIC(mod.thisrun,k=log(length(levels(sprabSc$plot))))

assign(modname.thisrun,mod.thisrun)

}
save.image("e93_wksp_sprab")

# Summarise results table:
sab2_res<-sab2_meta
sab2_res$BIC<-as.numeric(sab2_res$BIC)
sab2_res<-sab2_res[order(sab2_res$BIC),]
sab2_res<-tidy.df(sab2_res)
sab2_res$rank<-1:length(sab2_res[,1])
sab2_res

# write.table(sab2_res,file="sab2_res.txt",sep="\t",quote=F,row.names=F)

top.sab2<-get(paste(sab2_res$model[sab2_res$rank==1],"_mod",sep=""))
summary(top.sab2)

} # close STEP 2

## --- STEP 2 COEFFICIENTS & EFFECT SIZES:
{
# MODEL OUTPUTS are stored in spr2_res$model, plus "_mod" (e.g. Model1_13)
sab2_res

# Extract ALL coefficient tables:
for (i in 1:length(sab2_res[,1])){
coefname.thisrun<-paste(sab2_res$model[i],"_coef",sep="")
modname.thisrun<-paste(sab2_res$model[i],"_mod",sep="")
coef.thisrun<-coef.sum(get(modname.thisrun), est.type="logit")
assign(coefname.thisrun,coef.thisrun)
}

# Coef tables are stored in spr2_res$model, plus "_coef" (e.g. Model1_13). 
get(paste(sab2_res$model[1],"_coef",sep=""))
summary(top.sab2)

# Extract TOP MODEL coefficient table:
step2_topsab_coef<-get(paste(sab2_res$model[sab2_res$rank==1],"_coef",sep=""))

# write.table(step2_topsab_coef, "step2_topsab_coef.txt",sep="\t",row.names=F,quote=F)

# Plot EFFECT SIZES for all models:
for (i in 1:length(sab2_res[,1])){
coef.thisrun<-paste(sab2_res$model[i],"_coef",sep="")
quartz(file=paste(coef.thisrun,".pdf",sep=""),type="pdf")
heading.thisrun<-paste("Model = ",sab2_res$model[i],", Rank = ",sab2_res$rank[i],"/",length(sab2_res[,1]),", BIC = ", round(sab2_res$BIC[i],0),sep="")
efs.plot(get(coef.thisrun),heading.thisrun)
dev.off()
}

} # close STEP 2 coefficients

## --- STEP 2 TOP MODEL ESTIMATES:
{
# Modelled data:
head(sprabSc,2)
# Unscaled data for plotting:
head(sab.unscaled,2)

# TOP MODEL coefficient table:
step2_topsab_coef
sab2_res

# The model to plot:
summary(top.sab2)

### ---- NEW DATA FOR ESTIMATES:
{

head(step2_topsab_coef)

step2_topsab_coef$signif<-sign(step2_topsab_coef$lci.link)==sign(step2_topsab_coef$uci.link)

# Reduce to significant terms:
nd.sab2<-data.frame(term=step2_topsab_coef$term[step2_topsab_coef$signif==T])
# Remove intercept
nd.sab2<-data.frame(term=nd.sab2[-which(nd.sab2$term=="(Intercept)"),])
# Add an interaction identifier:
nd.sab2$interaction<-NA
if(length(grep(":",nd.sab2$term))>0) nd.sab2$interaction[grep(":",nd.sab2$term)]<-"int"
if(length(grep(":",nd.sab2$term))>0) nd.sab2$interaction[-grep(":",nd.sab2$term)]<-"main"
if(length(grep(":",nd.sab2$term))==0) nd.sab2$interaction<-"main"
nd.sab2<-tidy.df(nd.sab2)

# Add an x.axis term column for the estimates.nd function:
nd.sab2$xaxis.term<-c("yr","sdsub_abund","hgt_gap","prop_legume",rep("yr",4))

# Add an interaction term column that will be used by the estimates.nd function:
nd.sab2$int.term<-c(rep(NA,4),"sla_cwm","hgt_gap","prop_peren","prop_legume")

# For each term, add a name to be used as the data frame name:
nd.sab2$df.name<-c("yr.topsab2","sdabund.topsab2","hgtgap.topsab2","propleg.topsab2","yrsla.topsab2","yrhgtgap.topsab2","yrpropperen.topsab2","yrpropleg.topsab2")

# Get estimates with estimates.nd function for each data frame:
# Use first and third quartile as the range in the predictors is not normal
for (i in 1:length(nd.sab2[,1])){

if (nd.sab2$interaction[i]=="main"){
assign(nd.sab2$df.name[i],estimates.nd.V2("top.sab2",xaxis.term=nd.sab2$xaxis.term[i],modelled.data="sprabSc",unscaled.data="sab.unscaled"))
} # close if main

if (nd.sab2$interaction[i]=="int"){
assign(nd.sab2$df.name[i],estimates.nd.V2("top.sab2", xaxis.term=nd.sab2$xaxis.term[i],int.term=nd.sab2$int.term[i],modelled.data="sprabSc",unscaled.data="sab.unscaled",minmax=F))
} # close if interaction

} # close for data frame loop

save.image("e93_wksp_sprab")

} # close new data 

} # close top model estimates

## --- STEP 2 TOP MODEL PLOTS:
{

# Names for the estimate data frames are in nd.sab2$df.name:
nd.sab2$df.name[1]

# Convert nd.sab2 into a plot table:
# Use first and third quartile as the range in the predictors is not normal
step2sab.plot<-data.frame(data=nd.sab2$df.name,x.label=c("Year","Seeded subplot abundance","Height gap","Proportion legume","Year","Year","Year","Year"),x.variable=paste(nd.sab2$xaxis.term,"_unsc",sep=""),interaction=nd.sab2$interaction, x.subset=nd.sab2$int.term,x1.label=c(rep(NA,4),"SLA CWM 1st Q","Height gap 1st Q","Prop. perennial 1st Q","Prop. legume 1st Q"),x2.label=c(rep(NA,4),"SLA CWM 3rd Q","Height gap 3rd Q","Prop. perennial 3rd Q","Prop. legume 3rd Q"))

# Height gap and proportion legume are included in interactions, so remove them:
step2sab.plot<-step2sab.plot[-which(step2sab.plot$data %in% c("hgtgap.topsab2","propleg.topsab2")),]
step2sab.plot<-tidy.df(step2sab.plot)

# plots moved to separate section below

} # close top model plots

} # close spread abundance

## ~~~~ ---- ** COEF PLOTS ** ---- ~~~~ ##
# Plot multiple competing models on the same effect size plot, for SI
{

# Set up tables for each response separately:
{
# MODEL OUTPUTS are stored in step3_res$model, plus "_mod" (e.g. comb_6_mod). Add rank to model results table:
# Coef tables are stored in spr2_res$model, plus "_coef" (e.g. Model1_13). 

step3_res # colonisation prob step 2
colprob_BIC5<-step3_res[which(abs(step3_res$BIC[1]-step3_res$BIC)<5),]
get("6_4_coef")

spr2_res # spread step 2
spread_BIC5<-spr2_res[which(abs(spr2_res$BIC[1]-spr2_res$BIC)<5),]

cab2_res # colonisation abundance step 2
colab_BIC5<-cab2_res[which(abs(cab2_res$BIC[1]-cab2_res$BIC)<5),]
Comm_1_22_coef

sab2_res # spread abundance step 2
sab_BIC5<-sab2_res[which(abs(sab2_res$BIC[1]-sab2_res$BIC)<5),]

colprob_BIC5
spread_BIC5
colab_BIC5
sab_BIC5

# COLONISATION PROB. top models coef plot:

df_colpr<-do.call(rbind,lapply(paste(colprob_BIC5$model,"_coef",sep=""),get))

b_colpr=paste(colprob_BIC5$model,"_coef",sep="") # a vector with the coef tables to plot
# df1=do.call(rbind,lapply(b,function(x) rbind (get(x))))

yaxislab_colpr=unique(df_colpr$term)
coefall_colpr=data.frame(term= yaxislab_colpr)
# Add interaction identifyer:
coefall_colpr$int<-ifelse(rownames(coefall_colpr) %in% grep(":",coefall_colpr$term),"int","main")

# add term type:
coefall_colpr$type<-ifelse(rownames(coefall_colpr) %in% unlist(lapply(c("sla","hgt","seed","ldmc"),function(x) grep(x,coefall_colpr$term))),"trait","other")

# fix no_seeded:
coefall_colpr$type[coefall_colpr$term=="no_seeded"]<-"other"

coefall_colpr<-coefall_colpr[order(coefall_colpr$int,coefall_colpr$type,decreasing=c(TRUE,FALSE),method="radix"),]
coefall_colpr$term<-as.character(coefall_colpr$term)
coefall_colpr$term[coefall_colpr$type=="trait"]<-coefall_colpr$term[coefall_colpr$type=="trait"][order(coefall_colpr$term[coefall_colpr$type=="trait"])]
coefall_colpr$index<-1:nrow(coefall_colpr)

# Manually add labels
coefall_colpr$label<-c("Intercept","Year","Number seeded","Moisture","Soil NO3","Proportion perennial","Proportion grass","Proportion legume","Lifespan","Grass","Legume","Height absolute FD","LDMC community wt. mean","Seed relative FD","SLA relative FD","Proportion perennial x lifespan","Proportion grass x grass","Proportion grass x legume","Proportion legume x grass","Proportion legume x legume","Year x lifespan","Year x grass","Year x legume","Year x height absolute FD","Year x LDMC community wt. mean","Year x seed relative FD","Year x SLA relative FD")

# SPREAD top models coef plot:

df_sprpr<-do.call(rbind,lapply(paste(spread_BIC5$model,"_coef",sep=""),get))

b_sprpr=paste(spread_BIC5$model,"_coef",sep="") # a vector with the coef tables to plot
# df1=do.call(rbind,lapply(b,function(x) rbind (get(x))))

yaxislab_sprpr=unique(df_sprpr$term)
coefall_sprpr=data.frame(term=yaxislab_sprpr)
# Add interaction identifyer:
coefall_sprpr$int<-ifelse(rownames(coefall_sprpr) %in% grep(":",coefall_sprpr$term),"int","main")

# add term type:
coefall_sprpr$type<-ifelse(rownames(coefall_sprpr) %in% unlist(lapply(c("sla","hgt","seed","ldmc"),function(x) grep(x,coefall_sprpr$term))),"trait","other")

coefall_sprpr<-coefall_sprpr[order(coefall_sprpr$int,coefall_sprpr$type,decreasing=c(TRUE,FALSE),method="radix"),]
coefall_sprpr$term<-as.character(coefall_sprpr$term)
coefall_sprpr$term[coefall_sprpr$type=="trait"]<-coefall_sprpr$term[coefall_sprpr$type=="trait"][order(coefall_sprpr$term[coefall_sprpr$type=="trait"])]
coefall_sprpr$index<-1:nrow(coefall_sprpr)

# Manually add labels
coefall_sprpr$label<-c("Intercept","Year","Seeded subplot abundance","Lifespan","Grass","Legume","Height absolute FD","LDMC community wt. mean","LDMC raw","LDMC relative FD","Seed relative FD","SLA raw","SLA relative FD","Year x lifespan","Year x grass","Year x legume","Year x height absolute FD","Year x LDMC community wt. mean","Year x LDMC raw","Year x LDMC relative FD","Year x seed relative FD","Year x SLA raw","Year x SLA relative FD")

# COLONISATION ABUND. top models coef plot:

df_colab<-do.call(rbind,lapply(paste(colab_BIC5$model,"_coef",sep=""),get))

b_colab=paste(colab_BIC5$model,"_coef",sep="") # a vector with the coef tables to plot
# df1=do.call(rbind,lapply(b,function(x) rbind (get(x))))

yaxislab_colab=unique(df_colab$term)
coefall_colab=data.frame(term=yaxislab_colab)
# Add interaction identifyer:
coefall_colab$int<-ifelse(rownames(coefall_colab) %in% grep(":",coefall_colab$term),"int","main")

# add term type:
coefall_colab$type<-ifelse(rownames(coefall_colab) %in% unlist(lapply(c("sla","hgt","seed","ldmc"),function(x) grep(x,coefall_colab$term))),"trait","other")

# fix no_seeded:
coefall_colab$type[coefall_colab$term=="no_seeded"]<-"other"

coefall_colab<-coefall_colab[order(coefall_colab$int,coefall_colab$type,decreasing=c(TRUE,FALSE),method="radix"),]
coefall_colab$term<-as.character(coefall_colab$term)
coefall_colab$term[coefall_colab$type=="trait"]<-coefall_colab$term[coefall_colab$type=="trait"][order(coefall_colab$term[coefall_colab$type=="trait"])]
coefall_colab$index<-1:nrow(coefall_colab)

# Manually add labels
coefall_colab$label<-c("Intercept","Year","Number seeded","Moisture","Soil NO3","Proportion perennial","Proportion grass","Proportion legume","Height gap","LDMC community wt. mean","LDMC gap","LDMC raw","LDMC relative FD","Seed community wt. mean","SLA absolute FD","SLA relative FD")

# SPREAD ABUND. top models coef plot:

df_sprab<-do.call(rbind,lapply(paste(sab_BIC5$model,"_coef",sep=""),get))

b_sprab=paste(sab_BIC5$model,"_coef",sep="") # a vector with the coef tables to plot
# df1=do.call(rbind,lapply(b,function(x) rbind (get(x))))

yaxislab_sprab=unique(df_sprab$term)
coefall_sprab=data.frame(term=yaxislab_sprab)
# Add interaction identifyer:
coefall_sprab$int<-ifelse(rownames(coefall_sprab) %in% grep(":",coefall_sprab$term),"int","main")

# add term type:
coefall_sprab$type<-ifelse(rownames(coefall_sprab) %in% unlist(lapply(c("sla","hgt","seed","ldmc"),function(x) grep(x,coefall_sprab$term))),"trait","other")

# fix no_seeded:
coefall_sprab$type[coefall_sprab$term=="no_seeded"]<-"other"

coefall_sprab<-coefall_sprab[order(coefall_sprab$int,coefall_sprab$type,decreasing=c(TRUE,FALSE),method="radix"),]
coefall_sprab$term<-as.character(coefall_sprab$term)
coefall_sprab$term[coefall_sprab$type=="trait"]<-coefall_sprab$term[coefall_sprab$type=="trait"][order(coefall_sprab$term[coefall_sprab$type=="trait"])]
coefall_sprab$index<-1:nrow(coefall_sprab)

# Manually add labels
coefall_sprab$label<-c("Intercept","Year","Seeded subplot abundance","Proportion perennial","Proportion legume","Height gap","LDMC raw","Seed absolute FD","Seed raw","Seed relative FD","SLA community wt. mean","Year x proportion perennial","Year x proportion legume","Year x height gap","Year x LDMC raw","Year x seed absolute FD","Year x seed raw","Year x seed relative FD","Year x SLA community wt. mean")

} # close tables

# Put them together:

coefall_colpr; b_colpr
coefall_colab; b_colab
coefall_sprpr; b_sprpr
coefall_sprab; b_sprab

coefdf<-c("coefall_colpr","coefall_colab","coefall_sprpr","coefall_sprab")
bs<-c("b_colpr","b_colab","b_sprpr","b_sprab")
labs<-c("Establishment occupancy","Establishment abundance","Spread occupancy","Spread abundance")
dfs<-c("df_colpr","df_colab","df_sprpr","df_sprab")

# width for a4: 8.27
quartz(title="",width=10,height=11.69, dpi=64)

par(mfrow=c(2,2),mar=c(6,14,2,1), mgp=c(2.6,1,0),oma=c(0,0,0,0))

yincs=c(0.2,0.12,0.2,0.2)
greyincs=c(0.2,0.1,0.2,0.2)
legx<-c(-100,-3.2,0.5,-4)

for (j in 1:length(bs)){

yoffset=0
greyoffset<-1
coef.thisrun<-get(coefdf[j])
bs.thirun<-get(bs[j])
df.thisrun<-get(dfs[j])
lab.thisrun<-labs[j]

nm=length(bs.thirun)
ylaboffset=(yincs[j]/2)*(nm-1)

plot(1:nrow(coef.thisrun),1:nrow(coef.thisrun),type="n",xlim=c(min(df.thisrun$lci.link),max(df.thisrun$uci.link)), ylim=c(0,nrow(coef.thisrun)), xlab="Effect size", ylab="", yaxt="n", bty="l")
title(main=bquote(.("(")*.(letters[j])*")"~.(lab.thisrun)~Delta~BIC~"<"~5), cex.main=1, font.main=1,adj=0)

axis(side=2, at=(1:nrow(coef.thisrun))-ylaboffset, labels=rev(as.character(coef.thisrun$label)),las=1, cex.axis=0.8)
arrows(0,-1,0,50,code=0)

for (i in 1:length(bs.thirun)){
  coefthisrun=get(bs.thirun[i])
coefmerge=merge(coef.thisrun,coefthisrun,all.x=T,all.y=FALSE)
  coefmerge=coefmerge[order(coefmerge$index),]
  if(i>1) yoffset=yoffset+yincs[j]
   if(i>1) greyoffset=greyoffset-greyincs[j]
points(rev(coefmerge$est),(1:nrow(coefmerge)-yoffset),pch=20,col=i,cex=1) 
  arrows(rev(coefmerge$lci.link),(1:length(coefmerge[,1])-yoffset),rev(coefmerge$uci.link),(1:length(coefmerge[,1])-yoffset),code=0)
}
par(xpd=NA)
legend(legx[j],6,legend=paste("Rank ",seq(1,length(bs.thirun))," ", bs.thirun,sep=""),bty="n",pch=20,col=seq(1,length(bs.thirun),length.out=length(bs.thirun)),cex=0.8)
par(xpd=F)

} # close j models

} # close coef plots

## ~~~ -- ** PLOTS FOR PAPER ** -- ~~~ ##
## ~~~~ ---- ** AUG 2018 ** ---- ~~~~ ##
{
# COLONISATION PROBABILITY:

step2.plot2<-data.frame(data=nd.s3$df.name,x.label=c("Soil NO3","Relative SLA distance","Absolute height distance","Proportion perennial","Proportion Grass","Proportion Grass","Proportion Legume","Year","Year","Year","Year","Year"),x.variable=paste(nd.s3$xaxis.term,"_unsc",sep=""),interaction=nd.s3$interaction, x.subset=nd.s3$int.term,x1.label=c(rep(NA,3),"Perennial","Forb","Non-legume","Non-legume","Seed RelDist min","LDMC CWM min","Perennial","Forb","Non-legume"),x2.label=c(rep(NA,3),"Annual","Grass","Legume","Legume","Seed RelDist max","LDMC CWM max","Annual","Grass","Legume"))

quartz(title="",width=11.69,height=8.27, dpi=64, pointsize=19)
par(mfrow=c(3,4),mar=c(4,4,1.5,1), pch=20, las=1, bty="l", mgp=c(2.7,1,0),oma=c(0,0,1,0))

for (i in 1:length(step2.plot2[,1])){

data.thisrun<-get(as.character(step2.plot2$data[i]))
plot.x<-data.thisrun[,as.character(step2.plot2$x.variable[i])]
plot.y<-data.thisrun[,"fit"]

plot(plot.x, plot.y, type="n", ylab="Invader occupancy", xlab="", ylim=c(0,1))
if(step2.plot2$x.label[i]=="Soil NO3") title(xlab=expression(paste("Soil ", NO[3]),mgp=c(2.2,1,0))) else title(xlab=as.character(step2.plot2$x.label[i]),mgp=c(2.2,1,0))

if (step2.plot2$interaction[i]=="main"){
pg.ci(as.character(step2.plot2$x.variable[i]),as.character(step2.plot2$data[i]),colour=rgb(0,0,0,0.1))
lines(plot.x, plot.y)
} # close if main

if (step2.plot2$interaction[i]=="int"){
pg.ci(as.character(step2.plot2$x.variable[i]),as.character(step2.plot2$data[i]),x.subset=as.character(step2.plot2$x.subset[i]),colour=rgb(0,0,0,0.1))

sub1<-unique(data.thisrun[,as.character(step2.plot2$x.subset[i])])[1]
sub2<-unique(data.thisrun[,as.character(step2.plot2$x.subset[i])])[2]

linex.lev1<-data.thisrun[,as.character(step2.plot2$x.variable[i])][data.thisrun[,as.character(step2.plot2$x.subset[i])]==sub1]
liney.lev1<-data.thisrun$fit[data.thisrun[,as.character(step2.plot2$x.subset[i])]==sub1]

lines(linex.lev1, liney.lev1, lty=1)

linex.lev2<-data.thisrun[,as.character(step2.plot2$x.variable[i])][data.thisrun[,as.character(step2.plot2$x.subset[i])]==sub2]
liney.lev2<-data.thisrun$fit[data.thisrun[,as.character(step2.plot2$x.subset[i])]==sub2]

lines(linex.lev2, liney.lev2, lty=2)

par(xpd=NA)
legend(par("usr")[1],1.1,legend=c(as.character(step2.plot2$x1.label[i]),as.character(step2.plot2$x2.label[i])),lty=c(1,2),bty="n")
par(xpd=F)

} # close if interaction

par(xpd=NA)
text(par("usr")[1],par("usr")[4]+0.15,paste("(",letters[i],")",sep=""),cex=1.2)
par(xpd=F)

} # close plot for

# COLONISATION ABUNDANCE:

cab2b.plot2<-data.frame(data=nd.cab2b$df.name,x.label=c("Year","Soil NO3","Proportion perennial","Proportion legume","Height gap","Seed community wt mean","LDMC community wt mean"),x.variable=paste(nd.cab2b$xaxis.term,"_unsc",sep=""),interaction=nd.cab2b$interaction,x1.label=c(rep(NA,7)),x2.label=c(rep(NA,7)))

quartz(title="",width=11.69,height=(8.27/3)*2, dpi=64, pointsize=19)
par(mfrow=c(2,4),mar=c(4,4,1.5,1), pch=20, las=1, bty="l", mgp=c(3,1,0),oma=c(0,0,1,0))

for (i in 1:length(cab2b.plot2[,1])){

data.thisrun<-get(as.character(cab2b.plot2$data[i]))
y.limx<-0.1

plot.x<-data.thisrun[,as.character(cab2b.plot2$x.variable[i])]
if(length(grep(".resp",colnames(data.thisrun)))>0) plot.y<-data.thisrun[,"fit.resp"] else plot.y<-data.thisrun[,"fit"]
if(length(grep(".resp",colnames(data.thisrun)))>0) lci.thisrun<-data.thisrun[,"lci.resp"] else lci.thisrun<-data.thisrun[,"lci"]
if(length(grep(".resp",colnames(data.thisrun)))>0) uci.thisrun<-data.thisrun[,"uci.resp"] else uci.thisrun<-data.thisrun[,"uci"]

if (length(plot.x)==2){

if(is.factor(plot.x)==T){
plot(c(0,1), plot.y, type="n", ylab="Proportion cover", xlab="", ylim=c(0,y.limx), xaxt="n",xlim=c(-0.5,1.5))
points(c(0,1), plot.y)
arrows(c(0,1),lci.thisrun,c(0,1),uci.thisrun, angle=90, code=3, length=0.05)
par(xpd=NA)
axis(side=1, at=c(0,1), labels=c(as.character(cab2b.plot2$x1.label[i]),as.character(cab2b.plot2$x2.label[i])))
par(xpd=F)
} # close if factor

if(is.factor(plot.x)==F){
plot(plot.x, plot.y, type="n", ylab="Invader abundance", xlab="", ylim=c(0,y.limx), xaxt="n",xlim=c(-0.5,1.5))
points(plot.x, plot.y)
arrows(plot.x,lci.thisrun,plot.x,uci.thisrun, angle=90, code=3, length=0.05)
par(xpd=NA)
axis(side=1, at=plot.x, labels=c(as.character(cab2b.plot2$x1.label[i]),as.character(cab2b.plot2$x2.label[i])))
par(xpd=F)
} # close if not factor

} # close if length 2

if (length(plot.x)>2){

plot(plot.x, plot.y, type="n", ylab="Invader abundance", xlab="", ylim=c(0,y.limx))
par(xpd=NA)
if(cab2b.plot2$x.label[i]=="Soil NO3") title(xlab=expression(paste("Soil ", NO[3])),mgp=c(2.35,1,0)) else if(cab2b.plot2$x.label[i]=="Proportion perennial") title(xlab="Proportion perennial\nin community",mgp=c(3.2,1,0)) else if(cab2b.plot2$x.label[i]=="Proportion legume") title(xlab="Proportion legume\nin community",mgp=c(3.2,1,0)) else title(xlab=as.character(cab2b.plot2$x.label[i]),mgp=c(2.2,1,0))
par(xpd=F)

# ; if(cab2b.plot2$x.label[i]=="Proportion perennial") title(xlab=expression(paste("xx ", NO[3])),mgp=c(2.2,1,0))

if (cab2b.plot2$interaction[i]=="main"){
pg.ci(as.character(cab2b.plot2$x.variable[i]),as.character(cab2b.plot2$data[i]),colour=rgb(0,0,0,0.1))
lines(plot.x, plot.y)
} # close if main

if (cab2b.plot2$interaction[i]=="int"){
pg.ci(as.character(cab2b.plot2$x.variable[i]),as.character(cab2b.plot2$data[i]),x.subset=as.character(cab2b.plot2$x.subset[i]),colour=rgb(0,0,0,0.1))

sub1<-unique(data.thisrun[,as.character(cab2b.plot2$x.subset[i])])[1]
sub2<-unique(data.thisrun[,as.character(cab2b.plot2$x.subset[i])])[2]

linex.lev1<-data.thisrun[,as.character(cab2b.plot2$x.variable[i])][data.thisrun[,as.character(cab2b.plot2$x.subset[i])]==sub1]
liney.lev1<-data.thisrun$fit[data.thisrun[,as.character(cab2b.plot2$x.subset[i])]==sub1]

lines(linex.lev1, liney.lev1, lty=1)

linex.lev2<-data.thisrun[,as.character(cab2b.plot2$x.variable[i])][data.thisrun[,as.character(cab2b.plot2$x.subset[i])]==sub2]
liney.lev2<-data.thisrun$fit[data.thisrun[,as.character(cab2b.plot2$x.subset[i])]==sub2]

lines(linex.lev2, liney.lev2, lty=2)

par(xpd=NA)
legend("topleft",legend=c(as.character(cab2b.plot2$x1.label[i]),as.character(cab2b.plot2$x2.label[i])),lty=c(1,2),bty="n")
par(xpd=F)

} # close if interaction

} # close if > 2 levels

par(xpd=NA)
text(par("usr")[1],par("usr")[4]+0.015,paste("(",letters[i+12],")",sep=""),cex=1.2)
par(xpd=F)

} # close plot for

# SPREAD PROBABILITY:

step2bspr.plot2<-data.frame(data=nd.spr2b$df.name,x.label=c("Seeded subplot abundance","Absolute height distance","Grass","Year","Year","Year","Year"),x.variable=paste(nd.spr2b$xaxis.term,"_unsc",sep=""),interaction=nd.spr2b$interaction, x.subset=nd.spr2b$int.term,x1.label=c(rep(NA,2),"Forb","SLA RelDist min","Seed RelDist min","LDMC RelDist min","Perennial"),x2.label=c(rep(NA,2),"Grass","SLA RelDist max","seed RelDist max","LDMC RelDist max","Annual"))

quartz(title="",width=11.69,height=(8.27/3)*2, dpi=64, pointsize=19)
par(mfrow=c(2,4),mar=c(4,4,1.5,1), pch=20, las=1, bty="l", mgp=c(2.7,1,0),oma=c(0,0,1,0))

for (i in 1:length(step2bspr.plot2[,1])){

data.thisrun<-get(as.character(step2bspr.plot2$data[i]))
plot.x<-data.thisrun[,as.character(step2bspr.plot2$x.variable[i])]
plot.y<-data.thisrun[,"fit"]

if (length(plot.x)==2){

plot(plot.x, plot.y, type="n", ylab="Invader occupancy", xlab="", ylim=c(0,1), xaxt="n",xlim=c(-0.5,1.5))
points(plot.x, plot.y)
arrows(plot.x,data.thisrun$lci,plot.x,data.thisrun$uci, angle=90, code=3, length=0.05)
axis(side=1, at=plot.x, labels=c(as.character(step2bspr.plot2$x1.label[i]),as.character(step2bspr.plot2$x2.label[i])))
}

if (length(plot.x)>2){

plot(plot.x, plot.y, type="n", ylab="Invader occupancy", xlab="", ylim=c(0,1))
title(xlab=as.character(step2bspr.plot2$x.label[i]),mgp=c(2.2,1,0))

if (step2bspr.plot2$interaction[i]=="main"){
pg.ci(as.character(step2bspr.plot2$x.variable[i]),as.character(step2bspr.plot2$data[i]),colour=rgb(0,0,0,0.1))
lines(plot.x, plot.y)
} # close if main

if (step2bspr.plot2$interaction[i]=="int"){
pg.ci(as.character(step2bspr.plot2$x.variable[i]),as.character(step2bspr.plot2$data[i]),x.subset=as.character(step2bspr.plot2$x.subset[i]),colour=rgb(0,0,0,0.1))

sub1<-unique(data.thisrun[,as.character(step2bspr.plot2$x.subset[i])])[1]
sub2<-unique(data.thisrun[,as.character(step2bspr.plot2$x.subset[i])])[2]

linex.lev1<-data.thisrun[,as.character(step2bspr.plot2$x.variable[i])][data.thisrun[,as.character(step2bspr.plot2$x.subset[i])]==sub1]
liney.lev1<-data.thisrun$fit[data.thisrun[,as.character(step2bspr.plot2$x.subset[i])]==sub1]

lines(linex.lev1, liney.lev1, lty=1)

linex.lev2<-data.thisrun[,as.character(step2bspr.plot2$x.variable[i])][data.thisrun[,as.character(step2bspr.plot2$x.subset[i])]==sub2]
liney.lev2<-data.thisrun$fit[data.thisrun[,as.character(step2bspr.plot2$x.subset[i])]==sub2]

lines(linex.lev2, liney.lev2, lty=2)

par(xpd=NA)
legend(par("usr")[1],1.1,legend=c(as.character(step2bspr.plot2$x1.label[i]),as.character(step2bspr.plot2$x2.label[i])),lty=c(1,2),bty="n")
par(xpd=F)

} # close if interaction

} # close if > 2 levels

par(xpd=NA)
text(par("usr")[1],par("usr")[4]+0.15,paste("(",letters[i],")",sep=""),cex=1.2)
par(xpd=F)

} # close plot for

# SPREAD ABUNDANCE

quartz(title="",width=11.69,height=(8.27/3)*2, dpi=64, pointsize=19)
par(mfrow=c(2,4),mar=c(4,4,1.5,1), pch=20, las=1, bty="l", mgp=c(2.7,1,0),oma=c(0,0,1,0))
ylim_sprab<-c(0.06,0.2,rep(0.06,4))
let_ofs<-c(0.01,0.025,rep(0.01,4))

for (i in 1:length(step2sab.plot[,1])){

data.thisrun<-get(as.character(step2sab.plot$data[i]))
plot.x<-data.thisrun[,as.character(step2sab.plot$x.variable[i])]
if(length(grep(".resp",colnames(data.thisrun)))>0) plot.y<-data.thisrun[,"fit.resp"] else plot.y<-data.thisrun[,"fit"]
if(length(grep(".resp",colnames(data.thisrun)))>0) lci.thisrun<-data.thisrun[,"lci.resp"] else lci.thisrun<-data.thisrun[,"lci"]
if(length(grep(".resp",colnames(data.thisrun)))>0) uci.thisrun<-data.thisrun[,"uci.resp"] else uci.thisrun<-data.thisrun[,"uci"]

y.limx<-ylim_sprab[i]

if (length(plot.x)==2){

plot(plot.x, plot.y, type="n", ylab="Invader abundance", xlab="", ylim=c(0,y.limx), xaxt="n",xlim=c(-0.5,1.5))
points(plot.x, plot.y)
arrows(plot.x,data.thisrun$lci,plot.x,data.thisrun$uci, angle=90, code=3, length=0.05)
axis(side=1, at=plot.x, labels=c(as.character(step2sab.plot$x1.label[i]),as.character(step2sab.plot$x2.label[i])))
}

if (length(plot.x)>2){

plot(plot.x, plot.y, type="n", ylab="Invader abundance", xlab="", ylim=c(0,y.limx))
title(xlab=as.character(step2sab.plot$x.label[i]),mgp=c(2.2,1,0))

if (step2sab.plot$interaction[i]=="main"){
pg.ci(as.character(step2sab.plot$x.variable[i]),as.character(step2sab.plot$data[i]),colour=rgb(0,0,0,0.1))
lines(plot.x, plot.y)
} # close if main

if (step2sab.plot$interaction[i]=="int"){
pg.ci(as.character(step2sab.plot$x.variable[i]),as.character(step2sab.plot$data[i]),x.subset=as.character(step2sab.plot$x.subset[i]),colour=rgb(0,0,0,0.1))

sub1<-unique(data.thisrun[,as.character(step2sab.plot$x.subset[i])])[1]
sub2<-unique(data.thisrun[,as.character(step2sab.plot$x.subset[i])])[2]

linex.lev1<-data.thisrun[,as.character(step2sab.plot$x.variable[i])][data.thisrun[,as.character(step2sab.plot$x.subset[i])]==sub1]
liney.lev1<-data.thisrun$fit.resp[data.thisrun[,as.character(step2sab.plot$x.subset[i])]==sub1]

lines(linex.lev1, liney.lev1, lty=1)

linex.lev2<-data.thisrun[,as.character(step2sab.plot$x.variable[i])][data.thisrun[,as.character(step2sab.plot$x.subset[i])]==sub2]
liney.lev2<-data.thisrun$fit.resp[data.thisrun[,as.character(step2sab.plot$x.subset[i])]==sub2]

lines(linex.lev2, liney.lev2, lty=2)

par(xpd=NA)
legend(par("usr")[1],0.068,legend=c(as.character(step2sab.plot$x1.label[i]),as.character(step2sab.plot$x2.label[i])),lty=c(1,2),bty="n",cex=0.8)
par(xpd=F)

} # close if interaction

} # close if > 2 levels

par(xpd=NA)
text(par("usr")[1],par("usr")[4]+let_ofs[i],paste("(",letters[i+7],")",sep=""),cex=1.2)
par(xpd=F)

} # close plot for

} # close plots for paper








