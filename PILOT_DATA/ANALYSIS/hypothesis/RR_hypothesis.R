####################################################################################################
# R code for power analysis on the interaction and signle slope analyss effect of the PILOT DATA of: 
# "Does anxiety moderate training duration effects on habits in humans?  Determining the effects of 
#  trait anxiety on the experimental induction of habits in an instrumental outcome devaluation task"

## Last modified by Eva Pool on JANUARY 2019
####################################################################################################


# load libraries
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
pacman::p_load( ggplot2)

require(lattice)

# Set path
home_path       <- dirname(rstudioapi::getActiveDocumentContext()$path)
figures_path    <- file.path(home_path,'ANALYSIS','hypothesis', 'figures')


# create data to illustrate hypothesis
unit     <- c (3,      1,      1,      1)
training <- c('1-day','1-day','3-day','3-day')
anxiety  <- c('1: lower-level','2: higher-level','1: lower-level','2: higher-level')
database <- data.frame(unit, training, anxiety)


ann_text1 <- data.frame(training = "1-day", unit = 3.5, lab = "Text",
                        anxiety = factor('1: lower-level',levels = c("1: lower-level","2: higher-level")))

ann_text2 <- data.frame(training = "3-day", unit = 3.5, lab = "Text",
                        anxiety = factor('1: lower-level',levels = c("1: lower-level","2: higher-level")))

ann_text3 <- data.frame(training = "1-day", unit = 3.5, lab = "Text",
                       anxiety = factor('2: higher-level',levels = c("1: lower-level","2: higher-level")))

ann_text4 <- data.frame(training = "3-day", unit = 3.5, lab = "Text",
                        anxiety = factor('2: higher-level',levels = c("1: lower-level","2: higher-level")))


pp <- ggplot(database, aes (x = training, y = unit, fill = training, color = training)) +
  geom_bar(data =database, stat = "identity", alpha = .5) +
  facet_grid(~ anxiety) +
  ylim(0, 4)+ 
  geom_text(data = ann_text1,label = "+3", color = 'black', size = 10) +
  geom_text(data = ann_text2,label = "-1", color = 'black', size = 10) +
  geom_text(data = ann_text3,label = "-1", color = 'black', size = 10) +
  geom_text(data = ann_text4,label = "-1", color = 'black', size = 10) +
  
  theme_bw() +
  labs(
    title = '',
    x = 'General Anxiety', 
    y = "Behavioral change (Arbitrary Unit)"
  )
  
ppp <- pp + theme_linedraw(base_size = 14, base_family = "Helvetica")+ 
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.justification = c(1,1), legend.position = "right",
        legend.text = element_text(size = 14),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"))  


pdf (file.path(figures_path,'Figure_hypothesis.pdf'))
print(ppp)
dev.off()
