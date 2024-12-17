########################################################################################################################
# (4) FERTILITY AMONG MIGRANT DESCENDANTS (Descriptive statistics)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Nov 15th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(ggplot2);library(dplyr);library(tidyverse);library(survival);library(RColorBrewer);
library(scales);library(stringr)
theme_set(theme_bw()) # bw theme

# DATA LOADING ----
load('data3.person.year.RData') # indivyear
load('data.person.year.RData') # person.year

########################################################################################################################

# BIRTH ONE ----
# Let's focus on entry into parenthood (first birth)

# Observations before parenthood
cum_zero <- indivyear %>%
  filter(children_cum == 0)
# Observation with first birth
first_child <- indivyear %>%
  filter(children_cum != 0) %>% # remove observations while childless
  group_by(shnro) %>%
  arrange(age) %>% # sort by age (the earliest observation only)
  filter(!duplicated(shnro))
# And put it together
dat <- rbind(cum_zero,first_child)
dat <- dat %>% arrange(shnro) # rearrange by person

# The event of interest (irrespective of number of children had)
dat[,event := as.numeric(children_born != 0)]
dat[,table(event,children_born)]

# For these analyses, the coverage is from age 15 to 34, those who became parent at age 15 or earlier are not tracked
dat <- dat[age >= 15][!(age == 15 & event == TRUE)]
# As a result, our sample has...
nrow(dat) # 8,890,442 observations (person-year)
dat[,length(unique(shnro))] # From 691,184 individuals
# If we divide the by ancestry, generation, and gender...
dat[age == 15,table(background,mig,gender)]
dat[,sum(event)] # And we observe 230,738 transitions to parenthood

########################################################################################################################

# KAPLAN-MEIER CURVES ----
# Let's retrieve the Kaplan-Meier estimates
categ <- dat[,unique(mig)]
km_sum <- km <- vector('list',length(categ))
names(km_sum) <- names(km) <- categ

# We need to set the start
dat[,start := age - 1]

for(i in names(km)){
  # for native, we only stratify by gender
  if(i == 'Native-born'){
    km[[i]] <- survfit(Surv(start,age,event)~gender,conf.type='log',
                       conf.int=.95,type='kaplan-meier',error='greenwood',
                       data=dat[mig == i])
  # for migrant and descendants with stratify by gender and background
  }else{
    km[[i]] <- survfit(Surv(start,age,event)~gender+background,conf.type='log',
                       conf.int=.95,type='kaplan-meier',error='greenwood',
                       data=dat[mig == i])
  }
}

# Extract results 
for(i in names(km)){
  # Function to get the strata correctly
  strata_lengths <- sapply(km[[i]]$strata,function(x) x)
  strata_labels <- rep(names(km[[i]]$strata),strata_lengths)
  
  km_sum[[i]] <- data.table(age = km[[i]]$time,
                            risk = km[[i]]$n.risk,
                            event = km[[i]]$n.event,
                            surv = km[[i]]$surv,
                            lower = km[[i]]$lower,
                            upper = km[[i]]$upper,
                            mig = paste(i),
                            strata = strata_labels)
}

# Now split the strata info
km_sum <- do.call(rbind,km_sum)
km_sum <- km_sum %>%
  separate(strata,into=c('gender','background'),sep=', ') %>%
  separate(gender,into=c(NA,'gender'),sep='=') %>%
  separate(background,into=c(NA,'background'),sep='=') %>%
  mutate(background = str_trim(background)) # remove blank spaces

km_sum <- as.data.table(km_sum)
km_sum[is.na(background)]$background <- 'Finland'

# We can add hazard probabilities
km_sum[,hazard := event/risk] # hazards: event / risk set
km_sum[,se.h := sqrt(hazard*(1-hazard)/risk)] # SE of the hazard

# Visualization
km_sum[,mig := factor(mig,levels = levels(indivyear$mig))]
km_sum[,background := factor(background,levels = levels(indivyear$background))]

tiff(filename='fig5.tiff',
     width=30,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=km_sum %>% filter(age %in% 15:31),
       aes(x=age,y=surv,ymin=lower,ymax=upper,col=mig,fill=mig)) +
  geom_ribbon(alpha=1/3,linetype='dashed') +
  geom_line() +
  facet_grid(gender~background) +
  scale_color_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_fill_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  xlab("Age") + ylab("Survival probabilities") + labs(color=NULL,fill=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# How about the hazards?
tiff(filename='fig6.tiff',
     width=30,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=km_sum %>% filter(age %in% 15:31),
       aes(x=age,y=hazard,col=mig,fill=mig)) +
  geom_smooth(method='loess',se=FALSE) +
  geom_line(linetype='dashed') +
  facet_grid(gender~background) +
  scale_color_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_fill_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  xlab("Age") + ylab("Hazard probabilities") + labs(color=NULL,fill=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# Let's see the evolution of the risk sets
tiff(filename='fig7.tiff',
     width=34,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=km_sum %>% filter(age %in% 15:31),
       aes(x=age,y=risk,fill=mig)) +
  geom_bar(stat='identity',position='stack',col='grey10') +
  scale_fill_manual(values=c('red1',brewer.pal(5,"Greys")[-1])) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = 'K')) +
  facet_wrap(gender~background,scales='free_y',nrow=2) +
  xlab("Age") + ylab("Risk set") + labs(color=NULL,fill=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# Remove unnecessary objects
rm(list=setdiff(ls(),c('dat','person.year')))

# EXCLUDE 1G MIGRANTS ----

# From this point on, let's focus on migrant descendants (for migrants we miss lots of info)
dat <- dat [mig != '1G migrant']
# This leaves with a sample of...
nrow(dat) # 8,657,627 observations
dat[,length(unique(shnro))] # From 637,801 individuals
dat[,sum(event)] # 215,074

########################################################################################################################

# COVARIATES ----

# Let's complement the data with covariates
person.year <- do.call(rbind,person.year)
person.year <- person.year[shnro %in% unique(dat$shnro)]

# TIME-FIXED DIMENTIONS: INFO BASED ON FAMILY OF ORIGIN ----
timefixvars <- merge(dat[,.(shnro,vuosi,syntyv,gender,mig,country,background,age,event)],
                     person.year[,.(shnro,vuosi,syntyv,
                                smkunta,skunta,mkunta,kunta,kuntaryhm, # region, municipality, habitat
                                peas,akoko_k,lkm_k,pekoko_k,pety, # family properties and size
                                hape,sose)], # SES
             by=c('shnro','vuosi','syntyv'),all=TRUE)
timefixvars <- timefixvars %>% rename(age_parent = age) # age at parenthood

# Let's create a person-level dataset with the age at parenthood (if happened, or the fist observation if it did not) 
indiv <- timefixvars %>%
  arrange(desc(event)) %>% filter(!duplicated(shnro)) %>% 
  select(shnro,syntyv,gender,mig,country,background,event,age_parent)
indiv[event == 0]$age_parent <- NA
# This is to look at migrant descendants only
id_mig <- indiv[mig != 'Native-born',shnro] # 21,203 individuals

# I will always look at age 14. If the info is missing, then at ages 15 and 16
timefixvars[,age := vuosi - syntyv] # age of the person in the observation
timefixvars <- timefixvars[age %in% 14:17]

# Function to draw future info (up to two years)
trackplustwo <- function(data,var,newvar,max_lead=3){
  var <- enquo(var)
  newvar <- enquo(newvar)
  
  data <- data %>%
    arrange(shnro,age) %>% # order by person and age
    group_by(shnro) %>% # group by person
    mutate(!!newvar := !!var) # copy original variable
  
  for(i in 1:max_lead){
    data <- data %>%
      mutate(!!newvar := ifelse(
       is.na(!!newvar),
       lead(!!newvar,n=i),
       !!newvar
      ))
  }
  return(ungroup(data))
}

# GEOGRAPHY ----
# Region of residence, municipality of residence, and habitat at ages 14-17
timefixvars <- trackplustwo(timefixvars,mkunta,region) # region of residence
timefixvars <- trackplustwo(timefixvars,kunta,municipality) # municipality of residence
timefixvars <- trackplustwo(timefixvars,kuntaryhm,habitat) # habitat (urban (1), semi-urban (2), rural (3))
timefixvars <- as.data.table(timefixvars)

# Let's bring this info to the indiv object
indiv <- merge(indiv,
               timefixvars %>% arrange(age) %>% filter(!duplicated(shnro)) %>% select(shnro,region,municipality,habitat),
               by='shnro')

# Define variables as factors (add labels)
indiv[,municipality := factor(municipality)]
indiv[,habitat := factor(habitat,labels=c('Urban','Semi-urban','Rural'))]
indiv[,region := factor(region,labels=c('Uusimaa','Southwest Finland','Satakunta','Kanta-Häme','Pirkanmaa','Päijät-Häme',
                                        'Kymenlaakso','South Karelia','South Savo','Northern Savonia','North Karelia',
                                        'Central Finland','South Ostrobothnia','Ostrobothnia','Central Ostrobothnia',
                                        'North Ostrobothnia','Kainuu','Lapland','Åland'))]

# Inspect
indiv[,prop.table(table(mig,habitat),margin=1)*100] # migrant descendants more often from urban areas
cor.test(indiv$age_parent,as.integer(indiv$habitat)) # some correlation: rural areas (3) have lower age at parenthood

# How many migrant descendants lost?
indiv[mig != 'Native-born' & (is.na(region) | is.na(habitat)),
      length(shnro)] # 138 cases (less than 1%)

# FAMILY TYPE AND SIZE ----
# Eliminate info if the person's status in the family is head or spouse (we want the info when child)
timefixvars <- timefixvars[!(peas %in% c(1,2,4,5))]
# Let's use the find whether the person come from a single-parent household, and the number of children in the household
timefixvars <- trackplustwo(timefixvars,pety,family_type) # Family type
timefixvars <- trackplustwo(timefixvars,lkm_k,children) # Number of children in the family
timefixvars <- as.data.table(timefixvars)

# Let's bring this info to the indiv object
indiv <- merge(indiv,
               timefixvars %>% arrange(age) %>% filter(!duplicated(shnro)) %>% select(shnro,family_type,children),
               by='shnro',all.x=TRUE)

# Define variables as factors (add labels)
indiv[,family_type := factor(ifelse(family_type == 2,'Married couple',
                                    ifelse(family_type %in% 3:4,'Single parent','Cohabitating couple')))]

# Inspect
indiv[,prop.table(table(mig,family_type),margin=1)*100] # migrant descendants more often from single-parent families
indiv %>% group_by(family_type) %>% summarize(mean(age_parent,na.rm=TRUE)) # no clear association here
cor.test(indiv$age_parent,indiv$children) # neg correlation: individuals from large families have earlier entry into parenthood

# Visualization
tiff(filename='fig8.tiff',
     width=15,height=12,units='cm',compression='lzw',bg='white',res=1000)
ggplot(indiv %>% filter(!is.na(children)),aes(x=children,y=age_parent)) +
  geom_boxplot(aes(x=as.factor(children)),fill='grey80') +
  stat_smooth(method='glm',color='blue',linetype='dashed') +
  facet_grid(~gender) +
  xlab("Number of children in the parents' household around age 14-16") + ylab("Age at parenthood") +
  scale_x_discrete(labels=c('1'='1','2'='2','3'='3','4'='4','5'='5','6'='6','7'='7+')) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# How many migrant descendants lost?
indiv[mig != 'Native-born' & (is.na(region) | is.na(habitat) | is.na(family_type) | is.na(children)),
      length(shnro)] # 1,183 cases (about 6% of all migrant descendants)

# SOCIO-ECONOMIC STATUS OF THE FAMILY ----
# Profession (persons aged 0 to 15 are classified based on the reference person in the household) (two year lag)
timefixvars <- trackplustwo(timefixvars,sose,family_status)
# Tenure status of dwelling (one year lag)
timefixvars <- trackplustwo(timefixvars,hape,tenure_status)
timefixvars <- as.data.table(timefixvars)

# Let's bring this info to the indiv object
indiv <- merge(indiv,
               timefixvars %>% arrange(age) %>% filter(!duplicated(shnro)) %>% select(shnro,family_status,tenure_status),
               by='shnro',all.x=TRUE)

# Define variables as factors (add labels)
indiv <- indiv %>%
  mutate(family_status = case_when(
    family_status %in% c(10,20) ~ 'Self-employed',
    family_status %in% 31:34 ~ 'Upper-level employee',
    family_status %in% 41:44 ~ 'Lower-level employee',
    family_status %in% 51:59 ~ 'Manual worker',
    family_status == 60 ~ 'Student',
    family_status == 70 ~ 'Pensioner',
    family_status == 81 ~ 'Unemployed',
    family_status == 82 ~ 'Other (conscript)',
    family_status == 99 ~ NA,
    TRUE ~ NA
  ))

indiv <- indiv %>%
  mutate(tenure_status = case_when(
    tenure_status %in% 1:2 ~ 'Owner',
    tenure_status %in% 3:5 ~ 'Rental',
    tenure_status %in% 6:7 ~ 'Other',
    tenure_status == 9 ~ NA,
    TRUE ~ NA
  ))

# Inspect
indiv[,prop.table(table(mig,family_status),margin=1)*100] # migrant descendants more often from unemployed backgrounds
indiv %>% group_by(family_status) %>% summarize(mean(age_parent,na.rm=TRUE)) # no clear association here
indiv[,prop.table(table(mig,tenure_status),margin=1)*100] # migrant descendants more likely to rent
indiv %>% group_by(tenure_status) %>% summarize(mean(age_parent,na.rm=TRUE)) # no clear association here either

# How many migrant descendants lost?
indiv[mig != 'Native-born' & (is.na(region) | is.na(habitat) | is.na(family_type) | is.na(children) | is.na(tenure_status)),
      length(shnro)] # 1,361 (about 6% of all migrant descendants)
indiv[mig != 'Native-born' & (is.na(region) | is.na(habitat) | is.na(family_type) | is.na(children) | is.na(tenure_status) | is.na(family_status)),
      length(shnro)] # 3,597 cases... this is almost 17% of the 21,023 migrant descendants

########################################################################################################################

# TIME-VARYING DIMENSIONS ----
timevarvars <- merge(dat[,.(shnro,vuosi,syntyv,gender,mig,country,background,start,age,event,observed)],
                     person.year[,.(shnro,vuosi,
                                    ututku_aste,ututku_ala,optuki, # education
                                    ptoim1,tyke,tyokk,kturaha_k, # employment and income
                                    sivs)], # partnership status
                     by=c('shnro','vuosi'),all.x=TRUE)

# EDUCATION ---- 
# One year lag!
timevarvars[shnro %in% id_mig,table(age,ututku_aste,useNA = 'always')] # level
timevarvars[shnro %in% id_mig,table(age,ututku_ala,useNA = 'always')] # field
# Two year lag!
timevarvars[shnro %in% id_mig,table(age,optuki,useNA = 'always')] # student financial help

timevarvars <- timevarvars %>%
  mutate(education = case_when(
    ututku_aste == 3 ~ 'Upper secondary',
    ututku_aste %in% 4:5 ~ 'Post-secondary or short-cycle tertiary',
    ututku_aste %in% 6:8 ~ "Bachelor's, master's, or doctoral",
    TRUE ~ 'Less than secondary or unknown'
  ))
timevarvars[shnro %in% id_mig,table(age,education)] 

# (UN)EMPLOYMENT ----
# One year lag!
timevarvars[shnro %in% id_mig,table(age,ptoim1,useNA = 'always')]

timevarvars <- timevarvars %>%
  mutate(activity = case_when(
    ptoim1 == 11 ~ 'Employed',
    ptoim1 == 12 ~ 'Unemployed',
    ptoim1 == 22 ~ 'Student',
    ptoim1 == 24 ~ 'Pensioner',
    ptoim1 == 25 ~ 'Conscript',
    ptoim1 %in% c(21,99) ~ 'Outside labor force',
    TRUE ~ NA
  ))
timevarvars[shnro %in% id_mig,table(age,activity,useNA = 'always')]

# Months unemployed (this variables changed in 2005!!)
timevarvars[shnro %in% id_mig,table(age,tyke,useNA = 'always')]
timevarvars[shnro %in% id_mig,table(age,tyokk,useNA = 'always')] # first look more reliable

# If not unemployed and months unemployed equal NA, set to zero
timevarvars[shnro %in% id_mig,table(activity,tyke,useNA = 'always')] 
timevarvars[,unemployment := tyke]
timevarvars[activity != "Unemployed" & is.na(unemployment)]$unemployment <- 0
timevarvars[shnro %in% id_mig,table(activity,unemployment,useNA = 'always')] 

# DISPOSABLE INCOME ----
# one year lag!
timevarvars[shnro %in% id_mig & age == 15,table(kturaha_k,useNA = 'always')]
timevarvars[shnro %in% id_mig & age == 20,table(kturaha_k,useNA = 'always')]
timevarvars[shnro %in% id_mig & age == 25,table(kturaha_k,useNA = 'always')]
timevarvars[,income := kturaha_k]

# MARITAL STATUS ----
timevarvars[shnro %in% id_mig,table(age,sivs,useNA = 'always')]

timevarvars <- timevarvars %>%
  mutate(marital_status = case_when(
    sivs == 1 ~ 'Unmarried',
    sivs == 2 ~ 'Married, in partnership, separated',
    sivs == 4 ~ 'Divorced',
    sivs == 5 ~ 'Widowed',
    TRUE ~ NA
  ))
timevarvars[shnro %in% id_mig,table(age,marital_status,useNA = 'always')]

########################################################################################################################

# Save image
dat <- timevarvars[,.(shnro,vuosi,syntyv,gender,mig,country,background,start,age,event,observed,
                      education,activity,unemployment,income,marital_status)]

save(indiv,file='data4.person.RData')
save(dat,file='data4.person.year.RData')

########################################################################################################################