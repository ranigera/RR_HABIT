####################################################################################################
# R code for the PILOT DATA of: 
# "Does anxiety moderate training duration effects on habits in humans?  Determining the effects of 
#  trait anxiety on the experimental induction of habits in an instrumental outcome devaluation task"

## Last modified by Eva on NOVEMBER 2018
## Verified by ??
####################################################################################################


# ----------------------------------------- PRELIMINARY STUFF ---------------------------------------------------------------------------

# load libraries
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load(car, lme4, lmerTest, pbkrtest, ggplot2, dplyr, plyr, tidyr, multcomp, mvoutlier, HH, doBy, psych, pastecs, reshape, reshape2, 
               jtools, effects, compute.es, DescTools, MBESS, afex, ez, metafor, influence.ME,GPArotation)

require(lattice)


# Set path
home_path       <- '/Users/evapool/Documents/my_github/TASK_HABITS/RR/PILOT_DATA' # this will need to be made non-specific at the end (source the)
figures_path    <- file.path(home_path,'ANALYSIS','interindividual', 'figures')
utilities_path  <- file.path(home_path,'ANALYSIS','interindividual','R')
setwd (home_path)

# source my utilites
#source (file.path(utilities_path, 'getChangeIndex.R'))
#source (file.path(utilities_path, 'getClassicIndex.R'))
#source (file.path(utilities_path, 'makeIndividualDiffPlot.R'))
source (file.path(utilities_path, 'makeSplitGroupPlotCovariate.R'))


# get database
FULL <- read.delim(file.path(home_path,'DATABASES/PILOT_DATABASE.txt'), header = T, sep ='') # read in dataset

# define factors
FULL$site      <- factor(FULL$site)
FULL$ID        <- factor(FULL$ID)
FULL$session   <- factor(FULL$session)
FULL$run       <- factor(FULL$run)
FULL$trial     <- factor(FULL$trial)
FULL$cue       <- factor(FULL$cue)
FULL$prepost   <- factor(FULL$prepost)
FULL$group     <- factor(FULL$group)

# remove the baseline condition from the data
FULL <- subset(FULL, cue == 'Valued' | cue == 'Devalued')

# get the last run of the last training session and all the runs after satiation
DAY1   <- subset(FULL, group == '1-day')
DAY3   <- subset(FULL, group == '3-day')

DAY1 <- ddply(DAY1, .(ID), transform, averagePress  = mean(pressFreq[prepost=="pre"]))
DAY3 <- ddply(DAY3, .(ID), transform, averagePress  = mean(pressFreq[prepost=="pre"]))

C.DAY1 <- subset(DAY1, run == '2' | run == '3')
DAY3   <- subset(DAY3, session == '3') # we want the last day only
C.DAY3 <- subset(DAY3, run == '4' | run == '5')

CHANGE <- rbind(C.DAY1,C.DAY3)

# get variable of interest
CHANGE <- ddply(CHANGE, .(ID), transform, prePressAverage  = mean(normPressFreq[prepost == "pre"]))
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeBehav  = (mean(normPressFreq[prepost=="post" & cue=='Valued']) - mean(normPressFreq[prepost=="pre" & cue=='Valued'])) - (mean(normPressFreq[prepost=="post" & cue=='Devalued']) - mean(normPressFreq[prepost=="pre" & cue=='Devalued'])))
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeLiking = (mean(normLiking[prepost=="post" & cue=='Valued']) - mean(normLiking[prepost=="pre" & cue=='Valued'])) - (mean(normLiking[prepost=="post" & cue=='Devalued']) - mean(normLiking[prepost=="pre" & cue=='Devalued'])))


# get the test phase only
TEST <- subset(CHANGE, prepost == 'post')
TEST <- ddply(TEST, .(ID), transform, valueDiff = (mean(normPressFreq[cue=='Valued']) - mean(normPressFreq[cue=='Devalued'])))

# ----------------------------------------- EFFECTS OF OVER-TRAINING ON DEVALUATION SENSITIVITY -----------------------------------------

# subset by site
T.CALTECH = subset(TEST, site == 'Caltech1')
T.CALTECH2= subset(TEST, site == 'Caltech2')
T.HAMBURG = subset(TEST, site == 'Hamburg')
T.SYDNEY  = subset(TEST, site == 'Sydney')
T.TELAVIV = subset(TEST, site == 'Tel_Aviv')

############################################### outcome devaluation induced changes

# ------------------------  CALTECH1 
T.CALTECH.mean = ddply(T.CALTECH,.(ID,group,cue),summarise,normPressFreq=mean(normPressFreq), prePressAverage = mean(prePressAverage), valueDiff = mean(valueDiff))

int.caltech = ezANOVA(T.CALTECH.mean, dv = normPressFreq, wid = ID, within = .(cue), within_covariates = prePressAverage, between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue+prePressAverage+ Error (ID/cue), data = T.CALTECH))


# ------------------------  CALTECH2 
T.CALTECH2.mean = ddply(T.CALTECH2,.(ID,group,cue),summarise,normPressFreq=mean(normPressFreq), prePressAverage = mean(prePressAverage), valueDiff = mean(valueDiff))

int.caltech2 = ezANOVA(T.CALTECH2.mean, dv = normPressFreq, wid = ID, within = .(cue), within_covariates = prePressAverage, between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue+prePressAverage+ Error (ID/cue), data = T.CALTECH2))


# ------------------------  HAMBURG
T.HAMBURG.mean = ddply(T.HAMBURG,.(ID,group,cue),summarise,normPressFreq=mean(normPressFreq), prePressAverage = mean(prePressAverage), valueDiff = mean(valueDiff))

int.hamburg = ezANOVA(T.HAMBURG.mean, dv = normPressFreq, wid = ID, within = .(cue), within_covariates = prePressAverage, between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue+prePressAverage+ Error (ID/cue), data = T.HAMBURG))


# ------------------------  SYDNEY
T.SYDNEY.mean = ddply(T.SYDNEY,.(ID,group,cue),summarise,normPressFreq=mean(normPressFreq), prePressAverage = mean(prePressAverage), valueDiff = mean(valueDiff))

int.sydeny = ezANOVA(T.SYDNEY.mean, dv = normPressFreq, wid = ID, within = .(cue), within_covariates = prePressAverage, between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue+prePressAverage+ Error (ID/cue), data = T.SYDNEY))


# ------------------------  TELAVIV
T.TELAVIV.mean = ddply(T.TELAVIV,.(ID,group,cue),summarise,normPressFreq=mean(normPressFreq), prePressAverage = mean(prePressAverage), valueDiff = mean(valueDiff))

int.telaviv = ezANOVA(T.TELAVIV.mean, dv = normPressFreq, wid = ID, within = .(cue), within_covariates = prePressAverage, between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue+prePressAverage+ Error (ID/cue), data = T.TELAVIV))


# --------------------- FIGURE 1 (AND META-ANALYSIS)
# get the mean and the std 
T.CALTECH.mean <- subset(T.CALTECH.mean, cue!='Devalued')
T.CALTECH2.mean <- subset(T.CALTECH2.mean, cue!='Devalued')
T.HAMBURG.mean <- subset(T.HAMBURG.mean, cue!='Devalued')
T.SYDNEY.mean <- subset(T.SYDNEY.mean, cue!='Devalued')
T.TELAVIV.mean <- subset(T.TELAVIV.mean, cue!='Devalued')

estimate.caltech = summaryBy(valueDiff ~ group, data = T.CALTECH.mean,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.hamburg = summaryBy(valueDiff ~ group, data = T.HAMBURG.mean,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.caltech2 = summaryBy(valueDiff ~ group, data = T.CALTECH2.mean,
                              FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.sydney  = summaryBy(valueDiff ~ group, data = T.SYDNEY.mean,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.telaviv  = summaryBy(valueDiff ~ group, data = T.TELAVIV.mean,
                              FUN = function(x) { c(m = mean(x), s = sd(x)) } )

# build database for meta-analysis
site           = c ("Caltech (2017-Sept) "                        , "Hamburg (2018-Jan)"                           ,"Caltech2 (2018-May)"                           ,"Sydeny (2018-May)"                             ,"Tel-Aviv (2018-June)")
year           = c ("2017-sept"                                   , "2018-jan"                                     ,"2018-may"                                      ,"2018-may"                                      , "2018-june")
food_cons      = c ("bysession"                                   , "byrun"                                        , "byrun"                                        ,"byrun"                                         , "byrun")

mean_moderate  = c (estimate.caltech$valueDiff.m[1]               , estimate.hamburg$valueDiff.m[1]                , estimate.caltech2$valueDiff.m[1]                , estimate.sydney$valueDiff.m[1]                , estimate.telaviv$valueDiff.m[1]) # mean difference prepost for moderate training
mean_extensive = c (estimate.caltech$valueDiff.m[2]               , estimate.hamburg$valueDiff.m[2]                , estimate.caltech2$valueDiff.m[2]                , estimate.sydney$valueDiff.m[2]                , estimate.telaviv$valueDiff.m[2]) # mean difference prepost for extinsive trainig
std_moderate   = c (estimate.caltech$valueDiff.s[1]               , estimate.hamburg$valueDiff.s[1]                , estimate.caltech2$valueDiff.s[1]                , estimate.sydney$valueDiff.s[1]                , estimate.telaviv$valueDiff.s[1])
std_extensive  = c (estimate.caltech$valueDiff.s[2]               , estimate.hamburg$valueDiff.s[2]                , estimate.caltech2$valueDiff.s[2]                , estimate.sydney$valueDiff.s[2]                , estimate.telaviv$valueDiff.s[2])
n_moderate     = c (length(which(T.CALTECH.mean$group == '1-day')), length(which(T.HAMBURG.mean$group == '1-day')), length(which(T.CALTECH2.mean$group == '1-day')) , length(which(T.SYDNEY.mean$group == '1-day')) , length(which(T.TELAVIV.mean$group == '1-day')))
n_extensive    = c (length(which(T.CALTECH.mean$group == '3-day')), length(which(T.HAMBURG.mean$group == '3-day')) , length(which(T.CALTECH2.mean$group == '3-day')) , length(which(T.SYDNEY.mean$group == '3-day')) , length(which(T.TELAVIV.mean$group == '3-day')))

metadata = data.frame( site, year, food_cons, mean_moderate, mean_extensive, std_moderate, std_extensive, n_moderate, n_extensive)

# compute effect sizes
meta.data <- escalc(measure="SMD", m1i=mean_moderate, sd1i=std_moderate, n1i=n_moderate,
                    m2i=mean_extensive, sd2i=std_extensive, n2i=n_extensive, data=metadata)

# compute random-effect model
res <- rma.mv(yi, vi, random = ~ 1 | site, data=meta.data)

# plot
par(mar=c(4,4,1,2)) # decrease margins so the full space is used
par(cex=1, font=1)### switch to bold font
forest.plot <-forest(res,slab = (meta.data$site),xlim=c(-3,2),
                     ilab = cbind(meta.data$n_extensive, meta.data$n_moderate),
                     ilab.xpos=c(-1.3,-0.9), cex=1,
                     order=order(meta.data$site))


# add column headings to the plot
par(cex=1, font=4)### switch to bold font
text(-3, 6.2, "STUDY",  pos=4)
text( 2, 6.2, "SMD [95% CI]", pos=2)
par(cex=1, font=3)### switch to bold font
text(c(-1.35,-0.8), 6.2, c("N: day1 ", "N: day3 "))


# save the plot in the figures folder
dev.print(pdf, file.path(figures_path,'Between_Figure_forest.pdf'))
dev.off()



# ---------------------  LINEAR MIXED MODEL 
change.inter = lmer(normPressFreq~ group*cue + prePressAverage + site + (cue|ID) + (1|trial), data = TEST, REML=FALSE)
change.simple = lmer(normPressFreq ~ group+cue + prePressAverage + site + (cue|ID) + (1|trial), data = TEST, REML=FALSE)

anova(change.inter, change.simple)

# ----- assumptions check
plot(fitted(change.inter),residuals(change.inter)) # note heteroscedastisity and the impact of the 0 values
qqnorm(residuals(change.inter))
hist(residuals(change.inter))


# plot distribution of effect of interest to see how the 0 responses affected our tageted effect: 
bg_b = ddply(TEST,.(ID,group),summarise,valueDiff=mean(valueDiff))

behav =data.frame(bg_b$valueDiff)
behav$group = bg_b$group
behav$typeMeasure <- 'changeBehavior'
colnames(behav) [1] <- "valueDiff"


pp <- ggplot(behav, aes(valueDiff, fill = group)) +
  geom_histogram(aes(y=..density..),alpha=0.3,binwidth=0.2)+
  geom_density(alpha = 0.1)+
  facet_grid(~group)+
  theme_bw()+
  labs(
    title = '',
    x = 'Valued - Devalued in pressing Behavior',
    y = "Density"
  ) 

ppp <-  pp + theme_linedraw(base_size = 14, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.justification = c(1,1), legend.position = "right",
        legend.text = element_text(size = 14),
        axis.ticks.x=element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"))  

pdf(file.path(figures_path,'Between_Figure_histograms.pdf'))
print(ppp)
dev.off()

# ----------------------------------------- THE ROLE OF INTERINDIVIDUAL DIFFERENCES -----------------------------------------

# -------------------------------------------------------- STAI 
change.stai = lmer(normPressFreq~ group*cue*ANXIETY + prePressAverage + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+ANXIETY)^2 + prePressAverage + site + (1+cue*prepost|ID) + (1|trial), data = TEST, REML=FALSE)

anova(change.stai, change.basic)

# ----- follow-up simple slope approach

# Anxiety -1 SD
TEST$ANX_pSD <- scale(TEST$ANXIETY, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.pSD = lmer(normPressFreq~ group*cue*ANX_pSD + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)

anova(sslop.pSD)

# Anxiety +1 SD
TEST$ANX_mSD <- scale(TEST$ANXIETY, scale = T) - 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.mSD = lmer(normPressFreq~ group*cue*ANX_mSD + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)

anova(sslop.mSD)

# --------------------------------------------Figure 2

Anxiety.means <- aggregate(TEST$valueDiff, by = list(TEST$ID, TEST$group, TEST$site, TEST$ANX_pSD, TEST$ANX_mSD, TEST$ANXIETY), FUN='mean') # extract means
colnames(Anxiety.means) <- c('ID','group','site', 'ANX_pSD', 'ANX_mSD','ANXIETY', 'valueDiff')

# predicted MEANS
acqC1.aov      <- aov_car(valueDiff  ~ group*ANXIETY + prePressAverage+site +Error(ID), data = Anxiety.means, observed = c("ANXIETY"), factorize = F, anova_table = list(es = "pes"))
acqC1.adjmeans <- emmeans(acqC1.aov, specs = c("group"), by = "ANXIETY", at = list(ANXIETY= c(-1, 1)))
acqC1.adjmeans

# real means
TEST$traitGroup <- ntile(TEST$ANXIETY, 2)
TEST$traitGroup <- factor(TEST$traitGroup)
TEST.means <- aggregate(TEST$valueDiff, by = list(TEST$ID, TEST$group, TEST$site, TEST$BIS_total, TEST$TICS_CSSS, TEST$ANXIETY, TEST$traitGroup), FUN='mean') # extract means
colnames(TEST.means) <- c('ID','group','site','BIS_total', 'TICS_CSSS', 'ANXIETY','traitGroup','valueDiff')

extreamGroupData.means = subset(TEST.means, traitGroup == min(as.numeric(TEST$traitGroup), na.rm = T) | traitGroup == max(as.numeric(TEST$traitGroup), na.rm = T) )
extreamGroupData.means$traitGroup = mapvalues(extreamGroupData.means$traitGroup, from = c("1", "2"), to = c("1: Low-level", "2: High-level"))

makeSplitGroupPlotCovariate(extreamGroupData.means,'Anxiety Groups (STAI questionnaire)', 'Between_Figure_STAI.pdf', figures_path)


# -------------------------------------------------------- TICS
change.inter = lmer(normPressFreq~ group*cue*TICS_CSSS + prePressAverage + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+TICS_CSSS)^2 + prePressAverage + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)

anova(change.inter, change.basic)

# ----- follow up testing
# Stress -1 SD
TEST$STRESS_pSD <- scale(TEST$TICS_CSSS, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low on stress)
sslop.pSD = lmer(normPressFreq~ group*cue*STRESS_pSD + site + (cue|ID) + (1|trial), data = TEST, REML=FALSE)

# Stress +1 SD
TEST$STRESS_mSD <- scale(TEST$TICS_CSSS, scale = T) - 1 # here I'm going to test at - 1SD (so people that are high on stress)
sslop.mSD = lmer(normPressFreq~ group*cue*STRESS_mSD + site + (1+cue|ID) + (1|trial), data = TEST, REML=FALSE)

# -------------------------------------------------------- BIS 
change.inter = lmer(normPressFreq~ group*cue + prePressAverage+BIS_total + site + (cue|ID) + (1|trial), data = TEST, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+BIS_total)^2 + prePressAverage + site + (cue|ID) + (1|trial), data = TEST, REML=FALSE)

anova(change.inter, change.basic)



#----------------------------------------- PRINCIPLE COMPONENT ANALYSIS ---------------------------
  
#------------------------ EXTRACT COMPONENTS
  
# Check if questionnaire data are correlated between them
questionnaires <- aggregate(ANXIETY ~   BIS_total* TICS_CSSS, 
                              data = CHANGE, FUN = mean, na.action = na.pass)
r.questionnaires= cor(questionnaires, use = "pairwise.complete.obs")

# prepare database for the PCA
Q_ACP.means.ID <- aggregate(ANXIETY ~ ID * TICS_SOOV * TICS_PREPE * TICS_WODI * TICS_EXWO * TICS_LACK * TICS_SOTE * TICS_SOIS * TICS_WORY * TICS_WOOV * BIS_motor * BIS_attentional * BIS_nonplanning, 
                            data = CHANGE, FUN = mean, na.action = na.pass) # we do not include the total scales
Q_ACP.means <- Q_ACP.means.ID
Q_ACP.means$ID <- NULL

# quick look at the covarivance structure
r.subscale = cor(Q_ACP.means, use = "pairwise.complete.obs")
cor.plot(Q_ACP.means,numbers=TRUE,main="correlation matrix")
names(Q_ACP.means)[names(Q_ACP.means) == 'V1'] <- 'STAI'

# apply PCA
describe (Q_ACP.means)
pairs.panels(na.omit(Q_ACP.means))

# determine the number of factors
fa.parallel(Q_ACP.means) 
# apply PCA with varimax rotation
quest.1.pca <- psych::principal(Q_ACP.means, rotate="varimax", nfactors=4, scores=TRUE) # "none", "varimax" (Default), "quatimax", "promax", "oblimin", "simplimax", and "cluster"

print(quest.1.pca$loadings,cutoff = 0.3)

# create figure with PCA solution
fa.diagram(quest.1.pca)

#------------------------ USE COMPONENTS AS PREDICTORS

# extract the components
axes <- predict(quest.1.pca, Q_ACP.means, Q_ACP.means)   
# combine it with the participants ID
dat <- cbind(Q_ACP.means.ID, axes)

# label factors:
#names(dat)[names(dat) == 'RC1'] <- 'anxiety_isolation'
#names(dat)[names(dat) == 'RC2'] <- 'impulsivity'
#names(dat)[names(dat) == 'RC3'] <- 'social_pressure'
#names(dat)[names(dat) == 'RC4'] <- 'work_pressure'

# merge with the FULL database
PCA_CHANGE <- join (CHANGE,dat, type = "full")

# run model
change.inter = lmer(normPressFreq~ group*cue*(RC1+RC2+RC3+RC4)+prePressAverage+ site + (cue|ID) + (1|trial), data = PCA_CHANGE, REML=FALSE)
anova(change.inter) 

# ----- assumptions check
plot(fitted(change.inter),residuals(change.inter)) # show this to ben
qqnorm(residuals(change.inter))
hist(residuals(change.inter))

# ******************************** linear approach follow up

# Social Isolation -1 SD
PCA_CHANGE$ANX_pSD <- scale(PCA_CHANGE$RC1, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.pSD = lmer(normPressFreq~ group*cue*ANX_pSD + site + (cue|ID) + (1|trial), data = PCA_CHANGE, REML=FALSE)

anova(sslop.pSD)

# Social Isolation +1 SD
PCA_CHANGE$ANX_mSD <- scale(PCA_CHANGE$RC1, scale = T) - 1 # here I'm going to test at + 1SD (so people that are high in anxiety)
sslop.mSD = lmer(normPressFreq ~ group*cue*ANX_mSD + site + (cue|ID) + (1|trial), data = PCA_CHANGE, REML=FALSE)

anova(sslop.mSD)
