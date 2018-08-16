#------------------------------------#
#------~~~~~ SECTION 0 ~~~~~---------#
#------~~~~ E93 FUNCTIONS ~~~~-------#
#------------------------------------#

# Authors: Annabel Smith & Jane Catford, except where other source indicated

normalise<-function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))

# Tidy data frames: drop levels and re-assign rownames to subsetted data frames:
tidy.df<-function(df){
df<-droplevels(df)
rownames(df)<-1:length(df[,1])
return(df)
}

# Randomly check rows of a data frame:
check.rows<-function(df,n=6){df[sample(1:length(df[,1]),n),]}

# Analyse multicollinearity among many terms in a data frame. Version two has been updated so that no two terms of the same type (e.g. sla_cwm and sla_relFD) can be included in the same linear model:
mcl_v2<-function(datasets){

# datasets: vector containing the datasets to be analysed, IN QUOTES

res.list<-list()

for (m in 1:length(datasets)){

data.thisrun<-get(datasets[m])
head(data.thisrun,3)

cnames<-colnames(data.thisrun)
out.mat1<-matrix(data=NA, nrow=length(cnames), ncol=3)

for(i in 1:length(cnames)){
resp.thisrun<-cnames[i]
preds.thisrun<-cnames[-which(cnames==resp.thisrun)]

# remove any like terms:

# first select on terms with an underscore to get the prefix:
if(gregexpr("_",resp.thisrun)[[1]][1]>0){

resp.prefix<-substr(resp.thisrun,1,gregexpr("_",resp.thisrun)[[1]][1])

# is it a trait? Some non-trait terms also have underscores, so determine if it's one to be worried about:
is.trait<-resp.prefix %in% c("hgt_","seed_","ldmc_","sla_")

if(is.trait==T){

preds.thisrun<-preds.thisrun[-which(unlist(gregexpr(resp.prefix,preds.thisrun))>0)]

} # close if term is a trait

} # close if underscore

formula.thisrun<-paste(resp.thisrun,"~",paste(preds.thisrun,collapse="+"), sep="")
r2<-summary(lm(eval(parse(text=formula.thisrun)),data=data.thisrun))$r.squared
VIF<-my_vif(r2)
out.mat1[i,1]<-resp.thisrun
out.mat1[i,2]<-r2
out.mat1[i,3]<-VIF
}
res.r2<-data.frame(out.mat1)
colnames(res.r2)<-c("response","r2","VIF")
res.r2[,2]<-as.numeric(as.character(res.r2[,2]))
res.r2[,3]<-as.numeric(as.character(res.r2[,3]))

res.list[[m]]<-res.r2
}
return(data.frame(do.call(rbind,res.list)))
} # close mcl_v2

# Analyse multicollinearity among many terms in a data frame (original version includes like terms in same model (e.g. sla_cwm and sla_relFD)):
mcl<-function(datasets){
# datasets: vector containing the datasets to be analysed, IN QUOTES

res.list<-list()

for (m in 1:length(datasets)){

data.thisrun<-get(datasets[m])
head(data.thisrun)

cnames<-colnames(data.thisrun)
out.mat1<-matrix(data=NA, nrow=length(cnames), ncol=3)

for(i in 1:length(cnames)){
resp.thisrun<-cnames[i]
preds.thisrun<-cnames[-which(cnames==resp.thisrun)]

formula.thisrun<-paste(resp.thisrun,"~",paste(preds.thisrun,collapse="+"), sep="")
r2<-summary(lm(eval(parse(text=formula.thisrun)),data=data.thisrun))$r.squared
VIF<-my_vif(r2)
out.mat1[i,1]<-resp.thisrun
out.mat1[i,2]<-r2
out.mat1[i,3]<-VIF
}
res.r2<-data.frame(out.mat1)
colnames(res.r2)<-c("response","r2","VIF")
res.r2[,2]<-as.numeric(as.character(res.r2[,2]))
res.r2[,3]<-as.numeric(as.character(res.r2[,3]))

res.list[[m]]<-res.r2
}
return(data.frame(do.call(rbind,res.list)))
} # close mcl function

lmer.pval<-function(x) 2*(1-pnorm(abs(x))) # x is the t-value from the model summary

std.err<-function(x) sd(x)/sqrt(length(x))

# PREDICT function
# extends predictSE to calculate CIs and dataframe-ise predictions with 'new data'
pred<-function(model,new.data,se.fit=T,type="response"){

if(class(model)!="glmmadmb"){
pr1<-predictSE(model,new.data,se.fit=se.fit, type=type)
df1<-data.frame(new.data,fit=pr1$fit,se=pr1$se.fit, lci=pr1$fit-(1.96*pr1$se.fit), uci=pr1$fit+(1.96*pr1$se.fit))
} # close lme4 models

# defined for logit models only
if(class(model)=="glmmadmb"){
pr1<-predict(model,new.data,se.fit=se.fit, type="link")
df1<-data.frame(new.data,fit.link=pr1$fit,se.link=pr1$se.fit, lci.link=pr1$fit-(1.96*pr1$se.fit), uci.link=pr1$fit+(1.96*pr1$se.fit))	
df1$fit.resp<-round(invlogit(df1$fit.link),4)
df1$lci.resp<-round(invlogit(df1$lci.link),4)
df1$uci.resp<-round(invlogit(df1$uci.link),4)
} # close admb models

return(df1)

} # close predict function

# CI transform function
# Calculates CIs from SEs and back transforms them to the response scale. Backtransform after calculating the CIs. For binomial models, use invlogit (arm package), not exp to backtransform. 
CItrans<-function(data,est.type=NULL){

library(arm)

# data: enter data name without quotes
# est.type: the scale on which the estimates were estimated
# estimates should be called est
# SEs should be called se

# defined for the logit scale for binomial models, update for other types when necessary:
if (est.type=="logit"){
data$lci.link<-data$est-(data$se*1.96)
data$uci.link<-data$est+(data$se*1.96)
data$est.resp<-round(invlogit(data$est),4)
data$lci.resp<-round(invlogit(data$lci.link),4)
data$uci.resp<-round(invlogit(data$uci.link),4)
} # close logit

if (est.type=="log"){
data$lci.link<-data$est-(data$se*1.96)
data$uci.link<-data$est+(data$se*1.96)
data$est.resp<-round(exp(data$est),4)
data$lci.resp<-round(exp(data$lci.link),4)
data$uci.resp<-round(exp(data$uci.link),4)
} # close log

if (est.type=="normal"){
data$lci.link<-data$est-(data$se*1.96)
data$uci.link<-data$est+(data$se*1.96)
} # close normal

return(data)
} # close transform CI function

# Polygon CI function
pg.ci<-function(x,data,x.subset=NULL,colour){

# DEFINITION:
# x, data and x.subset must be in quotes
# x: a vector of data on the x-axis
# data: a data frame which includes x 
# currently, the confidence intervals must be specified as $lci and $uci, need to fix this... 
# update 13th Feb 2017: I haven't fixed the $uci and $lci problem, but I've just added in an if to deal with two different cases that distinguish $uci.resp, etc.
# x.subset: a vector of data for subsetting x. If there no subsets, omit this argument or set it at NULL to plot a single polygon. If there are subsets, enter the colname of data frame x to be used as the subset

# For testing:

# i=3
# x<-as.character(step2spr.plot$x.variable[i])
# data<-as.character(step2spr.plot$data[i])
# x.subset<-as.character(step2spr.plot$x.subset[i])
# colour<-rgb(0,0,0,0.1)
# i=1

if(is.null(x.subset)==T){
xx<-paste(data,"$",x,sep="")
# lci<- paste(data,"$lci",sep="")
# uci<- paste(data,"$uci",sep="")
if(length(grep(".resp",colnames(get(data))))>0) lci<-paste(data,"$lci.resp",sep="") else lci<-paste(data,"$lci",sep="")
if(length(grep(".resp",colnames(get(data))))>0) uci<-paste(data,"$uci.resp",sep="") else uci<-paste(data,"$uci",sep="")

xvec <- c(eval(parse(text=xx)), tail(eval(parse(text=xx)), 1), rev(eval(parse(text=xx))), eval(parse(text=xx))[1])
yvec <- c(eval(parse(text=lci)), tail(eval(parse(text=uci)), 1), rev(eval(parse(text=uci))), eval(parse(text=lci))[1])
polygon(xvec, yvec, col=colour, border=NA)
} # close if no subsets

if(is.null(x.subset)==F){

# Get data and vector that is used for subsetting:
data.withsubset<-get(data)
subset.all<-data.withsubset[,x.subset]

# Specify subs.levs: levels for factors, unique numbers for binary variables, and first and third quartiles for continuous variables

if(is.factor(subset.all)) subs.levs<-levels(subset.all)

if(is.factor(subset.all)==F) {

if(length(unique(subset.all))==2) subs.levs<-unique(subset.all) 

} # close if subset is not a factor

for (i in 1:length(subs.levs)){

sub.thisrun<-subs.levs[i]
x.thisrun<-data.withsubset[which(subset.all==sub.thisrun),x]
# lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
# uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]
if(length(grep(".resp",colnames(data.withsubset)))>0) lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci.resp"] else lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
if(length(grep(".resp",colnames(data.withsubset)))>0) uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci.resp"] else uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]

xvec <- c(x.thisrun, tail(x.thisrun, 1), rev(x.thisrun), x.thisrun[1])
yvec <- c(lci.thisrun, tail(uci.thisrun, 1), rev(uci.thisrun), lci.thisrun[1])
polygon(xvec, yvec, col=colour, border=NA)

} # close for sub levels i

} # close if subset present

} # close pg.ci

# Extract coefficient table from model summary
coef.ext<-function(model) {

if(class(model)[1]=="lmerMod"){
mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,3)])[[1]],round(summary(model)$coefficients[,c(1,2,3)],4))
mod.store<-tidy.df(mod.store)
mod.store$P<-lmer.pval(mod.store[,4])
mod.store<-mod.store[,c(1,2,3,5)]
colnames(mod.store)<-c("term","est","se","P")
return(mod.store)
} # close lmer

if(class(model)[1]=="glmmadmb"){
mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,4)])[[1]],round(summary(model)$coefficients[,c(1,2,4)],4))
mod.store<-tidy.df(mod.store)
colnames(mod.store)<-c("term","est","se","P")
return(mod.store)
} # close glmmadmb

if(class(model)[1]=="glmerMod"){
mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,4)])[[1]],round(summary(model)$coefficients[,c(1,2,4)],4))
mod.store<-tidy.df(mod.store)
colnames(mod.store)<-c("term","est","se","P")
return(mod.store)
} # close glmerMod

} # close coef.ext

# Summarise the coefficient table. This combines the CItrans and coef.ext functions:
coef.sum<-function(model,...){
coef.store<-coef.ext(model)
ci.store<-CItrans(coef.store,...)
return(ci.store)
}

# Effect size plots from coef table:
efs.plot<-function(coef.table,heading){
# coef.table is the output from coef.sum
# heading must be in quotes
par(mar=c(4,8,2,1), mgp=c(2.6,1,0))
plot(rev(coef.table$est),1:length(coef.table[,1]), pch=20, xlim=c(min(coef.table$lci.link),max(coef.table$uci.link)), xlab="Effect size", ylab="", yaxt="n", bty="l", main=heading, cex.main=1, font.main=1)
arrows(rev(coef.table$lci.link),1:length(coef.table[,1]),rev(coef.table$uci.link),1:length(coef.table[,1]),code=0)
arrows(0,0,0,50,code=0)
axis(side=2, at=1:1:length(coef.table[,1]), labels=rev(as.character(coef.table$term)),las=1, cex.axis=0.8)
}

# Generate new data and predictions (calls on pred()). V2 includes option to change min and max for numeric x numeric interactions.
estimates.nd.V2<-function(model,xaxis.term, int.term=NULL,modelled.data, unscaled.data,minmax=NULL){

# Set these for testing:
# model<-"top.sab2"
# modelled.data<-"sprabSc"
# unscaled.data<-"sab.unscaled"
# xaxis.term<-"yr"
# int.term<-"prop_peren"
# minmax<-T

# model: the model to plot, in quotes
# xaxis.term: the term which will vary along the x-axis. This can only be a single term as, for now, we will never plot three-way interactions. 
# int.term: the term which will be represented as an interaction. Currently the code is set to do these for all levels for factors and the first and third quartile for numeric interactions
# modelled data: the actual data that were modelled
# unscaled data: a data frame which is identical to the modelled data, only with variables unscaled
# minmax: logical. If TRUE, plot the int.term at the min and max. If FALSE, plot the in.term at the 1st and 3rd quartile

modx.whole<-get(model)
mod.data<-get(modelled.data)
unsc.data<-get(unscaled.data)

# All model terms to appear in nd:
mod.terms<-unique(unlist(strsplit(row.names(coefficients(summary(modx.whole))),":")))
# Remove intercept:
mod.terms<-mod.terms[-which(mod.terms=="(Intercept)")]
# Remove any terms not in the modelled data frame. These will be the factorial variables with levels pasted as a suffix. In the case of this analysis, this is only lifespan, so this is done half manually. 
notin<-mod.terms[-which(mod.terms %in% colnames(mod.data))]
if(length(notin)>0) {
mod.terms<-mod.terms[which(mod.terms %in% colnames(mod.data))]
# Manually add lifespan:
if(length(grep("lifespan",notin))>0) mod.terms<-c(mod.terms,"lifespan")
} # close if notin > 0

# Extract relevant columns from modelled data:
mtx<-mod.terms[-which(mod.terms==xaxis.term)]
ndx<-mod.data[,colnames(mod.data) %in% mtx]
head(ndx)

# Mean of all numeric columns:
# ndx.num<-data.frame(t(apply(ndx[,sapply(ndx,is.numeric)],2,mean)))

# All numeric columns zero:
ndx.num<-data.frame(t(apply(ndx[,sapply(ndx,is.numeric)],2,function(x) x<-0)))

# Manually add categorical variables:
# For lifespan use P
# Use zero for numeric
if("yr" %in% colnames(ndx)) ndx.num$yr<-0
if("c3" %in% colnames(ndx)) ndx.num$c3<-0
if("c4" %in% colnames(ndx)) ndx.num$c4<-0
if("lifespan" %in% colnames(ndx)) ndx.num$lifespan<-factor("P",levels=c("A","P"))
head(ndx.num)

# ----- BASE DATA:

# This includes all terms, with only the x-axis variable varying (i.e. not the interaction yet):

# Put together with term of interest:
tx<-paste(modelled.data,xaxis.term,sep="$")

if(is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])==F){
term.nd<-data.frame(seq(min(eval(parse(text=tx))),max(eval(parse(text=tx))),length.out=50))
colnames(term.nd)<-xaxis.term
term.nd<-cbind(term.nd,ndx.num)
# Add unscaled variable for plotting:
tx.unsc<-paste(unscaled.data,xaxis.term,sep="$")
unsc.df<-data.frame(seq(min(eval(parse(text=tx.unsc))),max(eval(parse(text=tx.unsc))),length.out=50))
colnames(unsc.df)<-paste(xaxis.term,"unsc",sep="_")
term.nd<-cbind(term.nd,unsc.df)
head(term.nd)
} 

if(is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])){
term.nd<-data.frame(factor(c(levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[1],levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[2]),levels=c("P","A")))
colnames(term.nd)<-xaxis.term
term.nd<-cbind(term.nd,ndx.num)
# Add unscaled variable for plotting (although you can't unscale a factor, the plotting code might need this):
tx.unsc<-paste(unscaled.data,xaxis.term,sep="$")
unsc.df<-data.frame(factor(c(levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[1],levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[2]),levels=c("P","A")))
colnames(unsc.df)<-paste(xaxis.term,"unsc",sep="_")
term.nd<-cbind(term.nd,unsc.df)
head(term.nd)
}

# If the term of interest is a main effect, has only two levels and is NOT a factor (e.g. c3 and c4), reduce the data down to those unique levels. 

if(length(unique(eval(parse(text=tx))))==2 & is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])==F){
tx.nd<-paste("term.nd",xaxis.term,sep="$")
term.nd<-term.nd[c(which(eval(parse(text=tx.nd))==min(eval(parse(text=tx.nd)))),which(eval(parse(text=tx.nd))==max(eval(parse(text=tx.nd))))),]
} # close if length 2
head(term.nd)

# If there is an interaction term, rbind the 50 rows together with itself, varying the level of the interaction:

if (is.null(int.term)==F){

term.nd1<-term.nd
term.nd2<-term.nd

if(is.factor(mod.data[,int.term])){

int.lev<-levels(term.nd[,int.term])
lev.toadd<-int.lev[-which(int.lev %in% term.nd1[,int.term])]

term.nd2[,int.term]<-factor(lev.toadd,levels=c("A","P"))

term.nd<-rbind(term.nd1,term.nd2)

} # close if int.term is factor

if(is.numeric(mod.data[,int.term])){

# If the variable is binary (e.g. c3 and c4), make a df with the int.term switched on and one with it switched off:

if (length(unique(mod.data[,int.term]))==2){

term.nd2[,which(colnames(term.nd2)==int.term)]<-1
term.nd<-rbind(term.nd1,term.nd2)

} # close if binary

# If the variable is continuous, plot min and max or 1st and 3rd quartile, depending on the value in minmax:

if (length(unique(mod.data[,int.term]))!=2){

if (minmax==T){

vmin<-summary(mod.data[,int.term])[1]
vmax<-summary(mod.data[,int.term])[6]

term.nd1[,int.term]<-vmin
term.nd2[,int.term]<-vmax

term.nd<-rbind(term.nd1,term.nd2)

} # close minmax == T

if (minmax==F){

vmin<-summary(mod.data[,int.term])[2]
vmax<-summary(mod.data[,int.term])[5]

term.nd1[,int.term]<-vmin
term.nd2[,int.term]<-vmax

term.nd<-rbind(term.nd1,term.nd2)
} # close minmax == F

} # close if continuous numeric

} # close if is.numeric

} # close if int.term !=NULL

term.pr<-pred(modx.whole,term.nd,se.fit=T,type="response")

return(term.pr)

} # close estimtes.nd2

# Generate new data and predictions (calls on pred())
estimates.nd<-function(model,xaxis.term, int.term=NULL,modelled.data, unscaled.data){

# Set these for testing:
# model<-"top.cab2"
# modelled.data<-"colabSc"
# unscaled.data<-"cab.unscaled"
# xaxis.term<-"yr"
# int.term<-NULL

# model: the model to plot, in quotes
# xaxis.term: the term which will vary along the x-axis. This can only be a single term as, for now, we will never plot three-way interactions. 
# int.term: the term which will be represented as an interaction. Currently the code is set to do these for all levels for factors and the first and third quartile for numeric interactions
# modelled data: the actual data that were modelled
# unscaled data: a data frame which is identical to the modelled data, only with variables unscaled

modx.whole<-get(model)
mod.data<-get(modelled.data)
unsc.data<-get(unscaled.data)

# All model terms to appear in nd:
mod.terms<-unique(unlist(strsplit(row.names(coefficients(summary(modx.whole))),":")))
# Remove intercept:
mod.terms<-mod.terms[-which(mod.terms=="(Intercept)")]
# Remove any terms not in the modelled data frame. These will be the factorial variables with levels pasted as a suffix. In the case of this analysis, this is only lifespan, so this is done half manually. 
notin<-mod.terms[-which(mod.terms %in% colnames(mod.data))]
if(length(notin)>0) {
mod.terms<-mod.terms[which(mod.terms %in% colnames(mod.data))]
# Manually add lifespan:
if(length(grep("lifespan",notin))>0) mod.terms<-c(mod.terms,"lifespan")
} # close if notin > 0

# Extract relevant columns from modelled data:
mtx<-mod.terms[-which(mod.terms==xaxis.term)]
ndx<-mod.data[,colnames(mod.data) %in% mtx]
head(ndx)

# Mean of all numeric columns:
# ndx.num<-data.frame(t(apply(ndx[,sapply(ndx,is.numeric)],2,mean)))

# All numeric columns zero:
ndx.num<-data.frame(t(apply(ndx[,sapply(ndx,is.numeric)],2,function(x) x<-0)))

# Manually add categorical variables:
# For lifespan use P
# Use zero for numeric
if("yr" %in% colnames(ndx)) ndx.num$yr<-0
if("c3" %in% colnames(ndx)) ndx.num$c3<-0
if("c4" %in% colnames(ndx)) ndx.num$c4<-0
if("lifespan" %in% colnames(ndx)) ndx.num$lifespan<-factor("P",levels=c("A","P"))
head(ndx.num)

# ----- BASE DATA:

# This includes all terms, with only the x-axis variable varying (i.e. not the interaction yet):

# Put together with term of interest:
tx<-paste(modelled.data,xaxis.term,sep="$")

if(is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])==F){
term.nd<-data.frame(seq(min(eval(parse(text=tx))),max(eval(parse(text=tx))),length.out=50))
colnames(term.nd)<-xaxis.term
term.nd<-cbind(term.nd,ndx.num)
# Add unscaled variable for plotting:
tx.unsc<-paste(unscaled.data,xaxis.term,sep="$")
unsc.df<-data.frame(seq(min(eval(parse(text=tx.unsc))),max(eval(parse(text=tx.unsc))),length.out=50))
colnames(unsc.df)<-paste(xaxis.term,"unsc",sep="_")
term.nd<-cbind(term.nd,unsc.df)
head(term.nd)
} 

if(is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])){
term.nd<-data.frame(factor(c(levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[1],levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[2]),levels=c("P","A")))
colnames(term.nd)<-xaxis.term
term.nd<-cbind(term.nd,ndx.num)
# Add unscaled variable for plotting (although you can't unscale a factor, the plotting code might need this):
tx.unsc<-paste(unscaled.data,xaxis.term,sep="$")
unsc.df<-data.frame(factor(c(levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[1],levels(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])[2]),levels=c("P","A")))
colnames(unsc.df)<-paste(xaxis.term,"unsc",sep="_")
term.nd<-cbind(term.nd,unsc.df)
head(term.nd)
}

# If the term of interest is a main effect, has only two levels and is NOT a factor (e.g. c3 and c4), reduce the data down to those unique levels. 

if(length(unique(eval(parse(text=tx))))==2 & is.factor(get(modelled.data)[,which(colnames(get(modelled.data))==xaxis.term)])==F){
tx.nd<-paste("term.nd",xaxis.term,sep="$")
term.nd<-term.nd[c(which(eval(parse(text=tx.nd))==min(eval(parse(text=tx.nd)))),which(eval(parse(text=tx.nd))==max(eval(parse(text=tx.nd))))),]
} # close if length 2
head(term.nd)

# If there is an interaction term, rbind the 50 rows together with itself, varying the level of the interaction:

if (is.null(int.term)==F){

term.nd1<-term.nd
term.nd2<-term.nd

if(is.factor(mod.data[,int.term])){

int.lev<-levels(term.nd[,int.term])
lev.toadd<-int.lev[-which(int.lev %in% term.nd1[,int.term])]

term.nd2[,int.term]<-factor(lev.toadd,levels=c("A","P"))

term.nd<-rbind(term.nd1,term.nd2)

} # close if int.term is factor

if(is.numeric(mod.data[,int.term])){

# If the variable is binary (e.g. c3 and c4), make a df with the int.term switched on and one with it switched off:

if (length(unique(mod.data[,int.term]))==2){

term.nd2[,which(colnames(term.nd2)==int.term)]<-1
term.nd<-rbind(term.nd1,term.nd2)

} # close if binary

# If the variable is continuous, plot min and max:

if (length(unique(mod.data[,int.term]))!=2){

vmin<-summary(mod.data[,int.term])[1]
vmax<-summary(mod.data[,int.term])[6]

term.nd1[,int.term]<-vmin
term.nd2[,int.term]<-vmax

term.nd<-rbind(term.nd1,term.nd2)

} # close if continuous numeric

} # close if is.numeric

} # close if int.term !=NULL

term.pr<-pred(modx.whole,term.nd,se.fit=T,type="response")

return(term.pr)

} # close estimtes.nd

# R MATRIX FOR WALD TEST. 
R.mat<-function(model,coef.name){
if (class(model)=="glmmadmb") {
	coef.tab<-summary(model)$coefficients
	mod.fr<-"model$frame"}
if (class(model)=="lmerMod") {
	coef.tab<-summary(model)$coefficients
	mod.fr<-"model@frame"}
if (class(model)=="glmerMod") {
	coef.tab<-summary(model)$coefficients
	mod.fr<-"model@frame"}
# STEP 1. Determine what type of variable (main or interaction) coef.name is. Get the names of the variables that make up the variable in question. t.names of main effects will = 1 and of interactions will = 2. t.names should never = 3 unless there are three way interactions.
t.mat<-attr(terms(model),"factors")
t.names<-names(which(t.mat[,colnames(t.mat)==coef.name]==1))

# STEP 2. Get the levels of the factor so that they can be appended to the coef.name:
if (length(t.names)==1) lev<-levels(eval(parse(text=paste(mod.fr,"$",coef.name,sep=""))))

if (length(t.names)==2){
# the classes of the interaction terms:
int.terms<-data.frame(term=t.names,class=c(class(eval(parse(text=paste(mod.fr,"$",t.names[1],sep="")))),class(eval(parse(text=paste(mod.fr,"$",t.names[2],sep=""))))))
factor.term<-int.terms$term[which(int.terms$class=="factor")]
if(length(factor.term)==1) lev<-levels(eval(parse(text=paste(mod.fr,"$",factor.term,sep=""))))
else stop("more than one factor in interaction term")
# This line was added to make the function work on Alice's treatment:yr models. NOTE THAT IT MIGHT NOT WORK FOR ALL INTERACTION MODELS:
lev<-paste(factor.term,lev,":",levels(factor.term)[which(levels(factor.term)!=factor.term)],sep="")
}

# STEP 3: generate the R matrix using the levels and the length of the levels:
if (length(t.names)==1){
R<-matrix(data=0,nrow=length(lev)-1,ncol=length(coef.tab[,1]))
for(k in 1:(length(lev)-1)){
R[k,which(rownames(coef.tab)==paste(coef.name,lev[k+1],sep=""))]<-1
}
}
if(length(t.names)==2){
R<-matrix(data=0,nrow=length(lev)-1,ncol=length(coef.tab[,1]))
for(k in 1:(length(lev)-1)){
R[k,which(rownames(coef.tab)==lev[k+1])]<-1
}
}
R<-R
} # close R.mat

# Blank plot. plot.new() does this but I don't understand the defaults; I defined this so I know what's where.
blankplot<-function()plot(1:10,1:10,bty="n",type="n",xaxt="n",yaxt="n",xlab="",ylab="")

# WALD TEST: to calculate p-vals from a factor with three or more levels:
# NOTE THAT q MUST BE THE LENGTH OF THE MATRIX IN QUESTION.
# Author: Wade Blanchard
Wald<-function(object,R,q) {
if (!is.matrix(R)) stop("Restrictions must be a matrix")
b<-fixef(object) 
vc<-vcov(object)
w<-t(R%*%b-q)%*%solve(R%*%vc%*%t(R))%*%(R%*%b-q)
pw<-1-pchisq(w[1],length(q))
return(invisible(list(chisq=as.vector(w),pvalue=pw)))
}

# From: http://glmm.wikidot.com/faq
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Variance inflation factor for multicollinearity analysis (from https://onlinecourses.science.psu.edu/stat501/node/347):
my_vif<-function(r2){1/(1-r2)}



