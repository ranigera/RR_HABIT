####################################################################################################
# R code for power analysis on the interaction and signle slope analyss effect of the PILOT DATA of: 
# "Does anxiety moderate training duration effects on habits in humans?  Determining the effects of 
#  trait anxiety on the experimental induction of habits in an instrumental outcome devaluation task"

## Last modified by Eva Pool on JANUARY 2019
## Verified by Rani Gera
####################################################################################################


# ----------------------------------------- PRELIMINARY STUFF --------------------------------------

# load libraries
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load(car, lme4, lmerTest, pbkrtest, ggplot2, dplyr, plyr, tidyr, multcomp, mvoutlier, HH, doBy, psych, pastecs, reshape, reshape2, 
               jtools, effects, simr, DescTools, MBESS, afex, ez, metafor, influence.ME,GPArotation)

require(lattice)


# Set path
home_path       <- '/Users/evapool/Documents/my_github/TASK_HABITS/RR/PILOT_DATA' # this will need to be made non-specific at the end (source the)
figures_path    <- file.path(home_path,'ANALYSIS','powerAnalysis', 'figures')
utilities_path  <- file.path(home_path,'ANALYSIS','interindividual','R')
setwd (home_path)

# source my utilites
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

#----------------------------------------------POWER ANALYSIS FOR INTERACTION------------------------------------------------

# remove variables containig nans for this analysis (for simr)
vars <- names(CHANGE) %in% c("ID","itemxcondition", "group","prepost","cue","ANXIETY","site","normPressFreq")
CHANGE.naomit <- CHANGE[vars] # DATABASE WITHOUT NANS
CHANGE.naomit <- na.omit(CHANGE.naomit)

# define our model
mdl.int =  lmer(normPressFreq~ group*cue*prepost*ANXIETY + itemxcondition + site + (1+cue*prepost+itemxcondition|ID), data = CHANGE.naomit, REML=FALSE)
#let's check that the simr we obtain the same main result
doTest(mdl.int, fcompare(~ (group+cue+prepost+ANXIETY)^3 + itemxcondition + site))

#let's calcualte the power we have in our sample
powerA <- powerSim(mdl.int, test = fcompare(~ (group+cue+prepost+ANXIETY)^3 + itemxcondition + site), nsim =50)
powerA # here we notice we only have 60% power, we need to simulate more participants

# let's add simulated participants to our sample to obtain more power
simulated.stai <- extend(mdl.int, along="ID", n=500) # let simulate double the amount of particpants
powerS <- powerSim(simulated.stai, test = fcompare(~ (group+cue+prepost+ANXIETY)^3 + itemxcondition + site), nsim = 500) 
powerS # Now we have enoungh power

#saveRDS(powerS, file.path(home_path,'ANALYSIS/powerAnalysis/powerInteraction.RDATA'))
save (powerS, file = file.path(home_path,'ANALYSIS/powerAnalysis/powerInteraction.RDATA'))

#let's now calcuate of how much we can reduce the sample to obtain a 90% power
powerD <- powerCurve(simulated.stai,along = "ID", test = fcompare(~ (group+cue+prepost+ANXIETY)^3 + itemxcondition + site), nsim =500)

save(powerD, file = file.path(home_path,'/ANALYSIS/powerAnalysis/powerInteractionCurve.RDATA'))

#----- figure
powerplot <- plot(powerD)
dev.print(pdf, file.path(figures_path,'Figure_powerInteraction.pdf'))

dev.off()

#------------ Figure for the manuscript
sink(file.path(home_path,'ANALYSIS/powerAnalysis/dataForPlotInteraction.txt'))
powerD # you will need to manually modify this file
sink()

db_plot <- read.delim(file.path(home_path,'ANALYSIS/powerAnalysis/dataForPlotInteraction.txt'), header = F, sep ='', skip = 2, nrows = 10) # read in dataset

db <- rename(db_plot, c("V1"="ID", "V2"="means", "V3"="CIlow","V4"="CIup"))


pp <- ggplot(db, aes(x = ID, y = means)) +
  geom_abline(slope=0, intercept=95, linetype = "dashed", color = "red", size = 1) +
  geom_abline(slope=0, intercept=80, linetype = "dashed", color = "gray", size = 0.5) +
  geom_line(data= db, stat = "identity",  color="#999999") + 
  geom_errorbar(data = db, stat = "identity", aes( ymin = CIlow , ymax = CIup), width= .3, color="black") + # add reduce witdh
  geom_point(data=db, stat = "identity", size = 4) +
  annotate('text', y = 97, x = 80, label = '95% power') +
  annotate('text', y = 83, x = 80, label = '80% power') +
  geom_vline(xintercept = 400, color= "gray") +
  theme_bw() +
  labs(
    title = '',
    x = "Number of participants",
    y = "Estimated power (%)"
  )

ppp <- pp + theme_linedraw(base_size = 12, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "#999999"),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "#999999"),
        legend.justification = c(1,1), legend.position = "right",
        legend.text  = element_text(size = 12),
        axis.title.y = element_text(size = 12)) 

pdf(file.path(figures_path,'Figure_power.pdf'))
print(ppp)
dev.off()


#----------------------------------------------POWER ANALYSIS FOR SINGLE SLOPE------------------------------------------------

# Anxiety -1 SD
CHANGE.naomit$ANX_pSD <- scale(CHANGE.naomit$ANXIETY, scale = T) + 1 # here I'm going to test at - 1SD (so people that are low in anxiety)
mdl.slp = lmer(normPressFreq~ group*cue*prepost*ANX_pSD + itemxcondition + site + (cue*prepost+itemxcondition|ID), data = CHANGE.naomit, REML=FALSE)

#let's check that the simr we obtain the same main result
doTest(mdl.slp, test=fixed("group:cue:prepost", method = "anova"), lmerTestType = 3)# this does not work it looks like it cannot read more than 

#let's calcualte the power we have in our sample
powerA1 <- powerSim(mdl.slp,test=fixed("group:cue:prepost", method = "anova"),lmerTestType=3, nsim =50)
powerA1 # here we notice we only have 62% power, we need to simulate more participants

# let's add simulated participants to our sample to obtain more power
simulated.complex <- extend(mdl.slp, along="ID", n=500) # let simulate double the amount of particpants
powerS1 <- powerSim(simulated.complex,test=fixed("group:cue:prepost", method = "anova"),lmerTestType=3, nsim =50)
powerS1 # Now we have enoungh power

save(powerS1, file = file.path(home_path,'ANALYSIS/powerAnalysis/powerSigleSlope.RDATA'))

#let's now calcuate of how much we can reduce the sample to obtain a 90% power
powerD1 <- powerCurve(simulated.complex,along = 'ID', test=fixed("group:cue:prepost", method = "anova"),lmerTestType=3, nsim =500)

save(powerD1, file = file.path(home_path,'ANALYSIS/powerAnalysis/powerSigleSlopeCurveß.RDATA'))

sink(file.path(home_path,'ANALYSIS/powerAnalysis/dataForPlotSingleSlope.txt'))
powerD1 # you will need to manually modify this file
sink()

#----- figure
powerplot <- plot(powerD1)
dev.print(pdf, file.path(figures_path,'Figure_powerSingleSlope.pdf'))
pdf (file.path(figures_path,'Figure_powerSingleSlope.pdf'))
print(powerplot)
dev.off()

ps <- lastResult()
ps$errors


