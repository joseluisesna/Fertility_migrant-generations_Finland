########################################################################################################################
# (5) FERTILITY AMONG MIGRANT DESCENDANTS (Matching)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Nov 25th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(ggplot2);library(dplyr);library(ggpubr);library(tidyverse);library(MatchIt);library(cobalt)
theme_set(theme_bw()) # bw theme

# DATA LOADING ----
load('data4.person.RData') # indiv
load('data4.person.year.RData') # dat

########################################################################################################################

# COARSENED EXACT MATCHING ----

# The "treatment group" are always migrants descendants (native-born people act as "control")
indiv[,treat := ifelse(background != 'Finland',1,0)]

# Let's define some variables as factors
indiv[,family_type := factor(ifelse(family_type %in% c('Married couple','Cohabitating couple'),'Couple','Single parent'))]
indiv[,family_status := ifelse(family_status %in% c('Student','Pensioner','Unemployed','Other (conscript)'),'Other or unknown',family_status)]
indiv[is.na(family_status)]$family_status <- 'Other or unknown' # quite some missing in this variable (this is convenient here)
indiv[,family_status := factor(family_status,levels = unique(indiv$family_status)[c(3,4,2,1,5)])] # Manual worker is ref category

# Remove individuals with NAs in variables used for matching
indiv[mig != 'Native-born',length(shnro)]
indiv[,table(background,mig)]
indiv <- indiv[complete.cases(indiv[,.(gender,syntyv,region,habitat,family_type,children,family_status)])]
indiv[mig != 'Native-born',length(shnro)] # we reduce from 21,203 migrant descendants to 20,020
indiv[,table(background,mig)]

# Let' create a subset of native individuals that resemble each ancestry group (Russia, Estonia, etc.)
indiv[,table(background)]
bckgs <- levels(indiv$background)
matched <- balance <- samples <- alist()

for(i in bckgs[-1]){
  samples[[i]] <- indiv[background %in% c('Finland',i)]
}

# let's check out the imbalance
for(i in names(samples)){
  balance[[i]] <- matchit(treat ~ gender + syntyv # gender and year of birth
                          + region + habitat 
                          + family_type + children + family_status,
                          data = samples[[i]],method=NULL,distance='glm',link='probit') # probit
}
lapply(balance,summary)

# CEM
for(i in names(samples)){
  matched[[i]] <- matchit(treat ~ gender + syntyv # gender and year of birth
                          + region + habitat 
                          + family_type + children + family_status,
                          data = samples[[i]],method='cem',k2k=TRUE)
}
lapply(matched,summary,un=FALSE)

# Save image
save(matched,file='data5.matched.RData')

#######################################################################################################################

# VISUALIZATOIN ----

# Function to visualize the matching results
visualize.match <- function(object){
  figs <- alist()

  # labels
  lbls <- c(paste("SES (",indiv[,rev(levels(family_status))],")",sep=''),
            "Number of children in parents' household", "Single parent",
            paste("Habitat (",indiv[,rev(levels(habitat))],')',sep=''),
            paste('Region (',indiv[,rev(levels(region))],')',sep=''),
            'Birth year', 'Gender (Woman)')
  
  figs[['balance']] <- love.plot(object,stats='m',abs=FALSE,
                                 drop.distance=FALSE,threshold=c(m=.01),
                                 shapes=c('circle filled','circle'),
                                 colors=c('grey30','grey10'),
                                 sample.names = c("All","Matched"),
                                 position = 'top') + scale_y_discrete(labels=lbls)
  
  # figure per variable
  figdat <- data.table(var = c('gender','syntyv','region','habitat','family_type','children','family_status'),
                       lbl = c('Gender','Birth year','Region in Finland','Habitat','Family type',
                               "Number of children in parents' household",'Socio-economic status'))
  
  for(i in 1:nrow(figdat)){
    figs[[paste(figdat$var[i])]] <- bal.plot(object,var.name=figdat$var[i],which='both',position='top') +
      scale_fill_manual(values=c('tomato','dodgerblue'),labels=c('Native-born','Migrant descendant')) +
      labs(title='',x=figdat$lbl[i],fill='') +
      theme(strip.background = element_rect(fill='black'),strip.text=element_text(color='white'),
            axis.text.x = element_text(angle=90))
    
  }
  # Output
  return(figs)
}

# Apply function to all 6 ancestry groups
vizz <- alist()
for(i in names(matched)){
  vizz[[i]] <- visualize.match(matched[[i]])
}

# Export results
for(i in names(matched)){
  tiff(filename=paste('fig.match.',i,'.tiff',sep=''),
       width=40,height=40,units='cm',compression='lzw',bg='white',res=1000)
  ggarrange(ggarrange(vizz[[i]][[1]],
                      ggarrange(vizz[[i]][[2]],vizz[[i]][[3]],vizz[[i]][[4]],
                                labels=LETTERS[2:4],ncol=1,common.legend = TRUE),
                      labels=c('A','','',''),ncol=2,widths = c(3,1)),
            ggarrange(vizz[[i]][[5]],vizz[[i]][[6]],vizz[[i]][[7]],vizz[[i]][[8]],
                      labels=LETTERS[5:8],nrow=1,common.legend = TRUE,legend='none'),
            nrow=2,heights = c(3,1))
  dev.off()
}

########################################################################################################################

# REDUCE THE PERSON-YEAR DATA TO MATCHED SAMPLE ----

# Add weights to individuals and remove unmatched individuals
for(i in names(samples)){
  samples[[i]][,weights := matched[[i]]$weights]
  samples[[i]] <- samples[[i]][weights == 1]
}

# Create different person-year datasets per ancestry group
person.year <- alist()
for(i in bckgs[-1]){
  person.year[[i]] <- dat[shnro %in% samples[[i]]$shnro]
}

# Add time-fixed dimensions to these person-year data
for(i in names(samples)){
  person.year[[i]] <- merge(person.year[[i]],
                            samples[[i]][,.(shnro,region,habitat,family_type,children,family_status,tenure_status)],
                            by='shnro',all.x=TRUE)
}

# Save image
save(person.year,file='data5.person.year.RData')

########################################################################################################################