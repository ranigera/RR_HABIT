####################################################################################################
# R code to compile database for PILOT DATA 
# "Does anxiety moderate training duration effects on habits in humans?  Determining the effects of 
# trait anxiety on the experimental induction of habits in an instrumental outcome devaluation task"

## Last modified by Eva on NOVEMBER 2018
## Verified by Rani Gera
####################################################################################################


# load libraries
if(!require(pacman)) {
install.packages("pacman")
library(pacman)
}
pacman::p_load( ggplot2, dplyr, plyr, tidyr, pastecs, reshape, reshape2 
               )

# Set path
full_path       <- dirname(rstudioapi::getActiveDocumentContext()$path)
pos             <- regexpr("PILOT_DATA", full_path)
home_path       <- substr(full_path, 1, pos+9)
utilities_path  <- file.path(home_path,'ANALYSIS','interindividual','R')
setwd (home_path)

source (file.path(utilities_path, 'normalizeVariablesBehavior.R'))
source (file.path(utilities_path, 'normalizeVariablesQuestionnaires.R'))

#####################################################################################################
# ----------------------------------------- QUESTIONNAIRES ------------------------------------------

# get questionnaiores databases
Q.CALTECH1 <- read.delim(file.path(home_path,'DATABASES/CALTECH_V1_QUESTIONNARIES.txt'), header = T, sep ='') # read in dataset
Q.CALTECH1 <- normalizeVariablesQuestionnaires(Q.CALTECH1) # we need to normalize for each center individually
Q.HAMBURG  <- read.delim(file.path(home_path,'DATABASES/HAMBURG_QUESTIONNARIES.txt'), header = T, sep ='')    # read in dataset
Q.HAMBURG  <- normalizeVariablesQuestionnaires(Q.HAMBURG) # we need to normalize for each center individually
Q.TELAVIV  <- read.delim(file.path(home_path,'DATABASES/TELAVIV_QUESTIONNARIES.txt'), header = T, sep ='')    # read in dataset
Q.TELAVIV  <- normalizeVariablesQuestionnaires(Q.TELAVIV) # we need to normalize for each center individually
Q.CALTECH2 <- read.delim(file.path(home_path,'DATABASES/CALTECH_V2_QUESTIONNARIES.txt'), header = T, sep ='') # read in dataset
Q.CALTECH2 <- normalizeVariablesQuestionnaires(Q.CALTECH2) # we need to normalize for each center individually

tmp1 = join (Q.CALTECH1, Q.HAMBURG, type = "full")
tmp2 = join (tmp1, Q.TELAVIV, type = "full")
QUESTIONNAIRES = join (tmp2, Q.CALTECH2, type = "full")
# remove participant that have more than one missing data in the questionnaire of interest
QUESTIONNAIRES <- subset (QUESTIONNAIRES, !ID == '330' & !ID == '200' & !ID == '168' & !ID == '222')


#####################################################################################################
# ----------------------------------------- FREE OPERANT TASK ---------------------------------------

P.CALTECH1 <- read.delim(file.path(home_path,'DATABASES/CALTECH_V1.txt'), header = T, sep ='') # read in dataset
P.CALTECH1 <- normalizeVariablesBehavior(P.CALTECH1) # we need to normalize for each center individually
P.HAMBURG  <- read.delim(file.path(home_path,'DATABASES/HAMBURG.txt'), header = T, sep ='') # read in dataset
P.HAMBURG  <- normalizeVariablesBehavior(P.HAMBURG) # we need to normalize for each center individually
P.TELAVIV  <- read.delim(file.path(home_path,'DATABASES/TELAVIV.txt'), header = T, sep ='') # read in dataset
P.TELAVIV  <- normalizeVariablesBehavior(P.TELAVIV) # we need to normalize for each center individually
P.CALTECH2 <- read.delim(file.path(home_path,'DATABASES/CALTECH_V2.txt'), header = T, sep ='') # read in dataset
P.CALTECH2 <- normalizeVariablesBehavior(P.CALTECH2) # we need to normalize for each center individually
P.SYDNEY   <- read.delim(file.path(home_path,'DATABASES/SYDNEY.txt'), header = T, sep ='') # read in dataset
P.SYDNEY   <- normalizeVariablesBehavior(P.SYDNEY ) # we need to normalize for each center individually

tmp1 = join (P.CALTECH1, P.HAMBURG, type = "full")
tmp2 = join (tmp1, P.CALTECH2, type = "full")
tmp3 = join (tmp2, P.TELAVIV, type = 'full')
FREEOPERANT =  join (tmp3, P.SYDNEY, type = "full")
# remove participant based on pre-reg criteria
FREEOPERANT <- subset (FREEOPERANT,!ID == '234') # caltech 2 extream 
FREEOPERANT <- subset (FREEOPERANT,!ID == '299'  & !ID == '334' & !ID == '341' & !ID == '310' & !ID == '304' & !ID == '322' & !ID == '326' & !ID == '352' & !ID == '356' & !ID == '360' & !ID == '301') # automated exclusions in Telaviv


# merge
FULL   <- join (QUESTIONNAIRES, FREEOPERANT, type = "full")

# print database
write.table(FULL,file.path(home_path,'DATABASES/PILOT_DATABASE.txt'),sep="\t",row.names=FALSE)