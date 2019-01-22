####################################################################################################
# R code for the PILOT DATA of: 
# "Does anxiety moderate training duration effects on habits in humans?  Determining the effects of 
# trait anxiety on the experimental induction of habits in an instrumental outcome devaluation task"

## Last modified by Eva on NOVEMBER 2018
## Verified by Rani Gera
####################################################################################################



# ----------------------------------------- PRELIMINARY STUFF ---------------------------------------------------------------------------

# load libraries
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load(car, lme4, lmerTest, pbkrtest, ggplot2, dplyr, plyr, tidyr, multcomp, mvoutlier, HH, doBy, psych, pastecs, reshape, reshape2, 
               jtools, effects, compute.es, DescTools, MBESS, afex, ez, metafor, influence.ME, GPArotation)

require(lattice)


# Set path
home_path       <- '/Users/evapool/Documents/my_github/TASK_HABITS/RR/PILOT_DATA' # this will need to be made non-specific at the end (source the)
figures_path    <- file.path(home_path,'ANALYSIS','interindividual', 'figures')
utilities_path  <- file.path(home_path,'ANALYSIS','interindividual','R')
setwd (home_path)

# source my utilites
source (file.path(utilities_path, 'getChangeIndex.R'))
source (file.path(utilities_path, 'getClassicIndex.R'))
source (file.path(utilities_path, 'makeIndividualDiffPlot.R'))
source (file.path(utilities_path, 'makeSplitGroupPlot.R'))
source (file.path(utilities_path, 'countTrialxCondition.R'))


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
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeBehav  = (mean(normPressFreq[prepost=="post" & cue=='Valued']) - mean(normPressFreq[prepost=="pre" & cue=='Valued'])) - (mean(normPressFreq[prepost=="post" & cue=='Devalued']) - mean(normPressFreq[prepost=="pre" & cue=='Devalued'])))
CHANGE <- ddply(CHANGE, .(ID), transform, normChangeLiking = (mean(normLiking[prepost=="post" & cue=='Valued']) - mean(normLiking[prepost=="pre" & cue=='Valued'])) - (mean(normLiking[prepost=="post" & cue=='Devalued']) - mean(normLiking[prepost=="pre" & cue=='Devalued'])))

# code itemxcondition
CHANGE <- ddply(CHANGE, .(ID,prepost), countTrialxCondition)


# ----------------------------------------- EFFECTS OF OVER-TRAINING ON DEVALUATION SENSITIVITY -----------------------------------------

# get total number of participants included 
count(CHANGE$ID) # note that Caltech2 used a slightly different protocol so there are less repeat per condition

# subset by site
C.CALTECH = subset(CHANGE, site == 'Caltech1')
C.CALTECH2= subset(CHANGE, site == 'Caltech2')
C.HAMBURG = subset(CHANGE, site == 'Hamburg')
C.SYDNEY  = subset(CHANGE, site == 'Sydney')
C.TELAVIV = subset(CHANGE, site == 'Tel_Aviv')

############################################### Manipulation check

# hunger
summary(aov(hunger ~ group*prepost + Error (ID/prepost), data = C.CALTECH))
summary(aov(hunger ~ group*prepost + Error (ID/prepost), data = C.CALTECH2))
summary(aov(hunger ~ group*prepost + Error (ID/prepost), data = C.HAMBURG))
summary(aov(hunger ~ group*prepost + Error (ID/prepost), data = C.SYDNEY))
summary(aov(hunger ~ group*prepost + Error (ID/prepost), data = C.TELAVIV))

# liking ratings
summary(aov(normLiking ~ group*cue*prepost + Error (ID/cue*prepost), data = C.CALTECH)) 
summary(aov(normLiking ~ group*cue*prepost + Error (ID/cue*prepost), data = C.CALTECH2)) 
summary(aov(normLiking ~ group*cue*prepost + Error (ID/cue*prepost), data = C.HAMBURG))
summary(aov(normLiking ~ group*cue*prepost + Error (ID/cue*prepost), data = C.SYDNEY)) 
summary(aov(normLiking ~ group*cue*prepost + Error (ID/cue*prepost), data = C.TELAVIV)) 

############################################### outcome devaluation induced changes

# ------------------------  CALTECH1 
C.CALTECH = subset(CHANGE, site == 'Caltech1')
CALTECH.index <- getChangeIndex(C.CALTECH)# aggregate based on pre-post 

int.caltech = ezANOVA(CALTECH.index, dv = pressFreq, wid = ID, within = .(cue), between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue*prepost + Error (ID/cue*prepost), data = C.CALTECH))

# ------------------------  CALTECH2 
C.CALTECH2 = subset(CHANGE, site == 'Caltech2')
CALTECH2.index <- getChangeIndex(C.CALTECH2)# aggregate based on pre-post 

int.caltech2 = ezANOVA(CALTECH2.index, dv = pressFreq, wid = ID, within = .(cue), between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue*prepost + Error (ID/cue*prepost), data = C.CALTECH2))


# ------------------------  HAMBURG
C.HAMBURG = subset(CHANGE, site == 'Hamburg')
HAMBURG.index <- getChangeIndex(C.HAMBURG)# aggregate based on pre-post 

int.hamburg = ezANOVA(HAMBURG.index, dv = pressFreq, wid = ID, within = .(cue), between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue*prepost + Error (ID/cue*prepost), data = C.HAMBURG))


# ------------------------  SYDNEY
C.SYDNEY = subset(CHANGE, site == 'Sydney')
SYDNEY.index <- getChangeIndex(C.SYDNEY)# aggregate based on pre-post 

int.sydeny = ezANOVA(SYDNEY.index, dv = pressFreq, wid = ID, within = .(cue), between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue*prepost + Error (ID/cue*prepost), data = C.SYDNEY))


# ------------------------  TELAVIV
C.TELAVIV = subset(CHANGE, site == 'Tel_Aviv')
TELAVIV.index <- getChangeIndex(C.TELAVIV)# aggregate based on pre-post 

int.telaviv = ezANOVA(TELAVIV.index, dv = pressFreq, wid = ID, within = .(cue), between = group, type = 3, detailed = T, return_aov = T) # quick check because aov uses a type 1 anova
summary(aov(normPressFreq ~ group*cue*prepost + Error (ID/cue*prepost), data = C.TELAVIV))

# --------------------- FIGURE 1 (AND META-ANALYSIS)
CALTECH.index2 <- ddply(CALTECH.index, .(ID), transform, pressFreq = pressFreq-pressFreq[cue=="Devalued"])
CALTECH.index2 <- subset(CALTECH.index2, cue!='Devalued')

CALTECH2.index2 <- ddply(CALTECH2.index, .(ID), transform, pressFreq = pressFreq-pressFreq[cue=="Devalued"])
CALTECH2.index2 <- subset(CALTECH2.index2, cue!='Devalued')

HAMBURG.index2 <- ddply(HAMBURG.index, .(ID), transform, pressFreq = pressFreq-pressFreq[cue=="Devalued"])
HAMBURG.index2 <- subset(HAMBURG.index2, cue!='Devalued')

SYDNEY.index2 <- ddply(SYDNEY.index, .(ID), transform, pressFreq = pressFreq-pressFreq[cue=="Devalued"])
SYDNEY.index2 <- subset(SYDNEY.index2, cue!='Devalued')

TELAVIV.index2 <- ddply(TELAVIV.index, .(ID), transform, pressFreq = pressFreq-pressFreq[cue=="Devalued"])
TELAVIV.index2 <- subset(TELAVIV.index2, cue!='Devalued')

# get the mean and the std 
estimate.caltech = summaryBy(pressFreq ~ group, data = CALTECH.index2,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.hamburg = summaryBy(pressFreq ~ group, data = HAMBURG.index2,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.caltech2 = summaryBy(pressFreq ~ group, data = CALTECH2.index2,
                              FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.sydney  = summaryBy(pressFreq ~ group, data = SYDNEY.index2,
                             FUN = function(x) { c(m = mean(x), s = sd(x)) } )

estimate.telaviv  = summaryBy(pressFreq ~ group, data = TELAVIV.index2,
                              FUN = function(x) { c(m = mean(x), s = sd(x)) } )

# build database for meta-analysis
site           = c ("Pasadena1 (2017-Sept) "                      , "Hamburg (2018-Jan)"                           ,"Pasadena2 (2018-May)"                          ,"Sydeny (2018-May)"                             ,"Tel-Aviv (2018-June)")
year           = c ("2017-sept"                                   , "2018-jan"                                     ,"2018-may"                                      ,"2018-may"                                      , "2018-june")
food_cons      = c ("bysession"                                   , "byrun"                                        , "byrun"                                        ,"byrun"                                         , "byrun")

mean_moderate  = c (estimate.caltech$pressFreq.m[1]               , estimate.hamburg$pressFreq.m[1]                , estimate.caltech2$pressFreq.m[1]                , estimate.sydney$pressFreq.m[1]                , estimate.telaviv$pressFreq.m[1]) # mean difference prepost for moderate training
mean_extensive = c (estimate.caltech$pressFreq.m[2]               , estimate.hamburg$pressFreq.m[2]                , estimate.caltech2$pressFreq.m[2]                , estimate.sydney$pressFreq.m[2]                , estimate.telaviv$pressFreq.m[2]) # mean difference prepost for extinsive trainig
std_moderate   = c (estimate.caltech$pressFreq.s[1]               , estimate.hamburg$pressFreq.s[1]                , estimate.caltech2$pressFreq.s[1]                , estimate.sydney$pressFreq.s[1]                , estimate.telaviv$pressFreq.s[1])
std_extensive  = c (estimate.caltech$pressFreq.s[2]               , estimate.hamburg$pressFreq.s[2]                , estimate.caltech2$pressFreq.s[2]                , estimate.sydney$pressFreq.s[2]                , estimate.telaviv$pressFreq.s[2])
n_moderate     = c (length(which(CALTECH.index2$group == '1-day')), length(which(HAMBURG.index2$group == '1-day')) , length(which(CALTECH2.index2$group == '1-day')) , length(which(SYDNEY.index2$group == '1-day')) , length(which(TELAVIV.index2$group == '1-day')))
n_extensive    = c (length(which(CALTECH.index2$group == '3-day')), length(which(HAMBURG.index2$group == '3-day')) , length(which(CALTECH2.index2$group == '3-day')) , length(which(SYDNEY.index2$group == '3-day')) , length(which(TELAVIV.index2$group == '3-day')))

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
dev.print(pdf, file.path(figures_path,'S_Figure_forest.pdf'))
dev.off()


# ---------------------  LINEAR MIXED MODEL 
change.inter = lmer(normPressFreq ~ group*cue*prepost + site + itemxcondition + (cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
change.simple = lmer(normPressFreq ~ (group+cue+prepost)^2 + site + itemxcondition + (cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(change.inter, change.simple)

# check 1 there is no difference before and that there is a difference after devaluation
PRE  <- subset(CHANGE, prepost == 'pre')
POST <- subset(CHANGE, prepost == 'post')

pre.check = lmer (normPressFreq ~ group * cue + site + itemxcondition + (cue+itemxcondition|ID), data = PRE)
anova(pre.check)

post.check =  lmer (normPressFreq ~ group * cue + site + itemxcondition + (cue+itemxcondition|ID), data = POST)
anova(post.check)

# check 2 learning trajectory in pre session since there is a strong effect of itemxcondition in PRE
bg_b = ddply(PRE,.(itemxcondition,group,cue),summarise,normPressFreq=mean(normPressFreq))

ggplot(bg_b, aes(itemxcondition, fill = cue, color = cue)) +
  geom_point(aes(y=normPressFreq),alpha=0.9)+
  geom_line(aes(y=normPressFreq),alpha=0.9)+
  facet_grid(~group)+
  theme_bw()+
  ylim (0,1)+
  labs(
    title = '',
    x = 'Trial',
    y = "Normalised pressing"
  ) 

# ----- assumptions check
plot(fitted(change.inter),residuals(change.inter)) # note heteroscedastisity and the impact of the 0 values
qqnorm(residuals(change.inter))
hist(residuals(change.inter))

# ---------------------- FIGURE 2
# plot distribution of effect of interest to see how the 0 responses affected our tageted effect: 
bg_b = ddply(CHANGE,.(ID,group),summarise,normChangeBehav=mean(normChangeBehav))

behav =data.frame(bg_b$normChangeBehav)
behav$group = bg_b$group
behav$typeMeasure <- 'changeBehavior'
colnames(behav) [1] <- "Normscore"


pp <- ggplot(behav, aes(Normscore, fill = group)) +
  geom_histogram(aes(y=..density..),alpha=0.3,binwidth=0.2)+
  geom_density(alpha = 0.1)+
  facet_grid(~group)+
  theme_bw()+
  labs(
    title = '',
    x = 'Normalized Change in Behavior',
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

pdf(file.path(figures_path,'Figure_histograms.pdf'))
print(ppp)
dev.off()

# ----------------------------------------- THE ROLE OF INTERINDIVIDUAL DIFFERENCES -----------------------------------------

# -------------------------------------------------------- STAI 
change.stai  = lmer(normPressFreq~ group*cue*prepost*ANXIETY + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+prepost+ANXIETY)^3 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
anova(change.stai, change.basic)

# ----- assumptions check
plot(fitted(change.stai),residuals(change.stai)) # note heteroscedastisity and the impact of the 0 values
qqnorm(residuals(change.stai))
hist(residuals(change.stai))

# ----- follow-up simple slope approach

# Anxiety -1 SD
CHANGE$ANX_pSD <- scale(CHANGE$ANXIETY, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.pSD = lmer(normPressFreq~ group*cue*prepost*ANX_pSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(sslop.pSD)

# complete model
mdl.complex = lmer(normPressFreq~ group + cue + prepost + ANX_pSD  + 
                     group:cue+ group:prepost + group:ANX_pSD +
                     cue:prepost+cue:ANX_pSD+prepost:ANX_pSD+
                     group:cue:prepost + group:cue:ANX_pSD + group:prepost:ANX_pSD +
                     cue:prepost:ANX_pSD + group:cue:prepost:ANX_pSD +
                     itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

# model without the 3-way interaction of interest
mdl.simple = lmer(normPressFreq~ group + cue + prepost + ANX_pSD  + 
                    group:cue+ group:prepost + group:ANX_pSD +
                    cue:prepost+cue:ANX_pSD+prepost:ANX_pSD -
                    group:cue:prepost + group:cue:ANX_pSD + group:prepost:ANX_pSD +
                    cue:prepost:ANX_pSD + group:cue:prepost:ANX_pSD +
                    itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(mdl.complex, mdl.simple)

# Anxiety +1 SD
CHANGE$ANX_mSD <- scale(CHANGE$ANXIETY, scale = T) - 1 # here I'm going to test at + 1SD (so people that are high in anxiety)
sslop.mSD = lmer(normPressFreq ~ group*cue*prepost*ANX_mSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(sslop.mSD)

# complete model
mdl.complex = lmer(normPressFreq~ group + cue + prepost + ANX_mSD  + 
                     group:cue+ group:prepost + group:ANX_mSD +
                     cue:prepost+cue:ANX_mSD+prepost:ANX_mSD+
                     group:cue:prepost + group:cue:ANX_mSD + group:prepost:ANX_mSD +
                     cue:prepost:ANX_mSD + group:cue:prepost:ANX_mSD +
                     itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

# model without the 3-way interaction of interest
mdl.simple = lmer(normPressFreq~ group + cue + prepost + ANX_mSD  + 
                    group:cue+ group:prepost + group:ANX_mSD +
                    cue:prepost+cue:ANX_mSD+prepost:ANX_mSD -
                    group:cue:prepost + group:cue:ANX_mSD + group:prepost:ANX_mSD +
                    cue:prepost:ANX_mSD + group:cue:prepost:ANX_mSD +
                    itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(mdl.complex, mdl.simple)


# --------------------------------------------Figure 2

Anxiety.means <- aggregate(CHANGE$normChangeBehav, by = list(CHANGE$ID, CHANGE$group, CHANGE$site, CHANGE$ANX_pSD, CHANGE$ANX_mSD, CHANGE$ANXIETY), FUN='mean', na.rm = T) # extract means
colnames(Anxiety.means) <- c('ID','group','site', 'ANX_pSD', 'ANX_mSD','ANXIETY', 'normChangeBehav')

# ADJUSTED MEANS this is not exactly what we test with the lmer
acqC1.aov      <- aov_car(normChangeBehav  ~ group*ANXIETY +Error(ID), data = Anxiety.means, observed = c("ANXIETY"), factorize = F, anova_table = list(es = "pes"))
acqC1.adjmeans <- emmeans(acqC1.aov, specs = c("group"), by = "ANXIETY", at = list(ANXIETY= c(-1, 1)))
acqC1.adjmeans

# ACTUAL DATA
CHANGE$traitGroup <- ntile(CHANGE$ANXIETY, 2)
CHANGE$traitGroup <- factor(CHANGE$traitGroup)
CHANGE.means <- aggregate(CHANGE$normChangeBehav, by = list(CHANGE$ID, CHANGE$group, CHANGE$site, CHANGE$BIS_total, CHANGE$TICS_CSSS, CHANGE$ANXIETY, CHANGE$traitGroup), FUN='mean') # extract means
colnames(CHANGE.means) <- c('ID','group','site','BIS_total', 'TICS_CSSS', 'ANXIETY','traitGroup','normChangeBehav')

extreamGroupData.means = subset(CHANGE.means, traitGroup == min(as.numeric(CHANGE$traitGroup), na.rm = T) | traitGroup == max(as.numeric(CHANGE$traitGroup), na.rm = T) )
extreamGroupData.means$traitGroup = mapvalues(extreamGroupData.means$traitGroup, from = c("1", "2"), to = c("1: Lower-level", "2: Higher-level"))

makeSplitGroupPlot(extreamGroupData.means,'Anxiety Groups (STAI questionnaire)', 'Figure_STAI.pdf', figures_path)

# stats for figure
CHANGE.medianLow  <- subset (CHANGE, traitGroup == '1')
CHANGE.medianHigh <- subset (CHANGE, traitGroup == '2')

median.low  = lmer(normPressFreq~ group*cue*prepost + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE.medianLow, REML=FALSE)
median.low0 = lmer(normPressFreq~ (group+cue+prepost)^2 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE.medianLow, REML=FALSE)
anova(median.low, median.low0)

median.high = lmer(normPressFreq~ group*cue*prepost + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE.medianHigh, REML=FALSE)
median.high0 = lmer(normPressFreq~ (group+cue+prepost)^2 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE.medianHigh, REML=FALSE)
anova(median.high, median.high0)

# constrast analysis
# all habitual but low anxious low training
CHANGE$contrast[CHANGE$traitGroup == '1' & CHANGE$group =='1-day'] <-  3
CHANGE$contrast[CHANGE$traitGroup == '1' & CHANGE$group =='3-day'] <- -1
CHANGE$contrast[CHANGE$traitGroup == '2' & CHANGE$group =='1-day'] <- -1
CHANGE$contrast[CHANGE$traitGroup == '2' & CHANGE$group =='3-day'] <- -1

mdl.constrast = lmer(normPressFreq~ contrast*cue*prepost + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
mdl.null      = lmer(normPressFreq~ (contrast+cue+prepost)^2 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(mdl.constrast, mdl.null)

# -------------------------------------------------------- TICS
change.inter = lmer(normPressFreq~ group*cue*prepost*TICS_CSSS + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+prepost+TICS_CSSS)^3 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(change.inter, change.basic)

# ----- follow up testing
# Stress -1 SD
CHANGE$STRESS_pSD <- scale(CHANGE$TICS_CSSS, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low on stress)
sslop.pSD = lmer(normPressFreq~ group*cue*prepost*STRESS_pSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

# Stress +1 SD
CHANGE$STRESS_mSD <- scale(CHANGE$TICS_CSSS, scale = T) - 1 # here I'm going to test at - 1SD (so people that are high on stress)
sslop.mSD = lmer(normPressFreq~ group*cue*prepost*STRESS_mSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

# -------------------------------------------------------- BIS 
change.inter = lmer(normPressFreq~ group*cue*prepost*BIS_total + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)
change.basic = lmer(normPressFreq ~ (group+cue+prepost+BIS_total)^3 + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE, REML=FALSE)

anova(change.inter, change.basic)
# ----------------------------------------- PRINCIPLE COMPONENT ANALYSIS ---------------------------

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
# save the plot in the figures folder
dev.print(pdf, file.path(figures_path,'S_Figure_PCA.pdf'))
dev.off()

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
change.inter = lmer(normPressFreq~ group*cue*prepost*(RC1+RC2+RC3+RC4) + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = PCA_CHANGE, REML=FALSE)
anova(change.inter) 

# ----- assumptions check
plot(fitted(change.inter),residuals(change.inter)) # show this to ben
qqnorm(residuals(change.inter))
hist(residuals(change.inter))

# ******************************** linear approach follow up

# Social Isolation -1 SD
PCA_CHANGE$ANX_pSD <- scale(PCA_CHANGE$RC1, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
sslop.pSD = lmer(normPressFreq~ group*cue*prepost*ANX_pSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = PCA_CHANGE, REML=FALSE)

anova(sslop.pSD)

# Social Isolation +1 SD
PCA_CHANGE$ANX_mSD <- scale(PCA_CHANGE$RC1, scale = T) - 1 # here I'm going to test at + 1SD (so people that are high in anxiety)
sslop.mSD = lmer(normPressFreq ~ group*cue*prepost*ANX_mSD + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = PCA_CHANGE, REML=FALSE)

anova(sslop.mSD)

# ------------------------------- figure for social isolation and work load

# get means for the factors of interests
PCA_CHANGE.means <- aggregate(PCA_CHANGE$normChangeBehav, by = list(PCA_CHANGE$ID, PCA_CHANGE$group, PCA_CHANGE$site, PCA_CHANGE$RC1, PCA_CHANGE$RC4), FUN='mean') # extract means
colnames(PCA_CHANGE.means) <- c('ID','group','site', 'RC1', 'RC4','normChangeBehav')

# figure for RC1: Anxiety-Isolation
PCA_CHANGE.means$traitGroup <- ntile(PCA_CHANGE.means$RC1, 2)
PCA_CHANGE.means$traitGroup <- factor(PCA_CHANGE.means$traitGroup)

extreamGroupData.means = subset(PCA_CHANGE.means, traitGroup == min(as.numeric(PCA_CHANGE$traitGroup), na.rm = T) | traitGroup == max(as.numeric(PCA_CHANGE$traitGroup), na.rm = T) )
extreamGroupData.means$traitGroup = mapvalues(extreamGroupData.means$traitGroup, from = c("1", "3"), to = c("1: Low-level", "2: High-level"))

makeSplitGroupPlot(extreamGroupData.means,'Anxiety-Isolation Groups', 'S_Figure_PCAsocialIsolation.pdf', figures_path)

# figure for RC4: Work pressure
PCA_CHANGE.means$traitGroup <- ntile(PCA_CHANGE.means$RC4, 2)
PCA_CHANGE.means$traitGroup <- factor(PCA_CHANGE.means$traitGroup)

extreamGroupData.means = subset(PCA_CHANGE.means, traitGroup == min(as.numeric(PCA_CHANGE$traitGroup), na.rm = T) | traitGroup == max(as.numeric(PCA_CHANGE$traitGroup), na.rm = T) )
extreamGroupData.means$traitGroup = mapvalues(extreamGroupData.means$traitGroup, from = c("1", "2"), to = c("1: Low-level", "2: High-level"))

# makeIndividualDiffPlot(PCA_CHANGE.means, PCA_CHANGE.means$RC1,'Work Pressure', '','Figure_workPressure.pdf', figures_path)
makeSplitGroupPlot(extreamGroupData.means,'Work Pressure Groups', 'S_Figure_PCAworkLoad.pdf', figures_path)

