########################################################################################################################
# (3) FERTILITY AMONG MIGRANT DESCENDANTS (Person-year format and sample selection)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Nov 9th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(ggplot2);library(dplyr);library(tidyverse);library(survival);library(RColorBrewer);library(scales)
theme_set(theme_bw()) # bw theme

# DATA LOADING ----
load('data.sample.RData') # indiv 
load('data.children.RData') # children
load('data.person.year.RData') # person.year
load('data2.sample.with.backgrounds.RData') # sample

########################################################################################################################

# MERGE INDIV AND BACKGROUNDS ----
indiv[,gender := factor(sukup,levels=1:2,labels=c('Man','Woman'))] # Let's call gender by the English terms
indiv <- merge(indiv[,.(shnro,gender)], # ID and gender
               sample[,.(shnro,syntyv,syntyp2,mig,ageatmigr,country,background,language)],
               by='shnro',all.x=TRUE)

########################################################################################################################

# PERSON-YEAR DATASET ----
# We need to know the entry and exit points of focal subjects based on mobility or demise

# Let's start by the person-year data altogether
person.year <- do.call(rbind,person.year)
# Let's find out if the person past away, and in which year
deceased <- person.year %>% 
  arrange(vuosi) %>%
  filter(!is.na(kuolv)) %>%
  filter(!duplicated(shnro)) %>%
  select(c(shnro,sukup,vuosi))
deceased[,table(sukup)] # 4,397 men and 1,630 women passed away
names(deceased) <- c('shnro','sukup','demise')

# Also, what is the last observation for every individual
lastobs <- person.year %>% 
  group_by(shnro) %>%
  summarize(lastobs = max(vuosi))

# Let's put all this together
indiv <- merge(indiv,deceased[,.(shnro,demise)],by='shnro',all.x=TRUE)
indiv <- merge(indiv,lastobs,by='shnro',all.x=TRUE)

# Using the dates of dead and last observation, we can determine the exit points for each person
indiv <- indiv %>% 
  mutate(exityear = pmin(demise,lastobs,na.rm=TRUE)) %>%
  mutate(end = exityear - syntyv)

# With the exit points, we can create person-year dataset
indivyear <- survSplit(data=indiv[,.(shnro,syntyv,gender,syntyp2,mig,country,background,language,ageatmigr,end)],
                       cut=1:indiv[,max(end)],
                       end='end',event='end',start='age')
indivyear <- as.data.table(indivyear) # data.table format
indivyear[,age := age + 1]
indivyear[,vuosi := syntyv + age] # we connect with other dataset through year (vuosi) not age

# CHILDREN ALLOCATION ----
# Now, let's connect the sample to their children

# First we turn the data into births per subject per year
parentyear <- children %>%
  filter(shnro %in% indiv$shnro) %>% # only children of focal subjects
  group_by(shnro,syntyv) %>%
  summarise(children_born = n())
nrow(parentyear) # 442,125 transitions 

# Merge data
indivyear <- merge(indivyear,parentyear,by.x=c('shnro','vuosi'),by.y=c('shnro','syntyv'),all.x=TRUE)
indivyear[is.na(children_born)]$children_born <- 0 # if no children that year, use NA
indivyear[,sum(children_born != 0)] # We keep 440,813 transitions

# Now, let's see the cumulative number of children (no matter their biological state)
indivyear <- indivyear %>% 
  group_by(shnro) %>%
  arrange(age) %>% # start is the age
  mutate(children_cum = cumsum(children_born)) %>%
  ungroup()

# Now, we can remove observations during childhood
indivyear <- as.data.table(indivyear)
indivyear <- indivyear[age >= 14] # keep 14 for lagging
indivyear[,sum(children_born != 0)] # Only a few transition lost: 440,798 transitions

# With all this, we can have age at parenthood and detect cases of individuals becoming parents before 15
earlyparent <- unique(indivyear[children_cum != 0 & age <= 15]$shnro) # 207 cases
indiv[shnro %in% earlyparent,table(mig,gender)] # Mostly Native born and migrants

# Let' remove observation before arrival to migrants
indivyear <- indivyear[is.na(ageatmigr) | (!is.na(ageatmigr) & age >= ageatmigr)]

# How many person year observations we have here?
indivyear[,table(mig,gender)]
indivyear[children_cum == 0,table(mig,gender)] # and before turning parents?

# IDENTIFY PERSON-YEAR OBSERVATIONS ABROAD ----

# Let' identify spells abroad
spellsobserved <- person.year[,.(shnro,vuosi)]
spellsobserved[,observed := TRUE]
indivyear <- merge(indivyear,spellsobserved,by=c('shnro','vuosi'),all.x=TRUE)
indivyear[is.na(observed)]$observed <- FALSE # unobserved spells

# How many and affecting whom?
indivyear[,table(observed),by=.(gender)] # 31051 spells for men (0.5%) and 35436 for women (0.7%)
indivyear[observed == FALSE,table(mig,gender)]

########################################################################################################################

# SAMPLE SELECTION ---

# How many individuals are lost in the process (not observed between ages 15 and 34)?
sum(!(indiv$shnro %in% unique(indivyear[age %in% 15:34 & observed == TRUE]$shnro))) # 6,290 cases

# Let's remove those not observed in the person-year data
indiv <- indiv[shnro %in% unique(indivyear[age %in% 15:34 & observed == TRUE]$shnro)] # 736,706
indivyear <- indivyear[shnro %in% unique(indivyear[age %in% 15:34 & observed == TRUE]$shnro)] # 11,056,521

# Let's remove those for whom we do not know their background
indiv <- indiv[!is.na(country)] # 735,501
indivyear <- indivyear[!is.na(country)] # 11,036,665 

# How many cases we have by generation and gender
indiv[,table(mig,gender)];indivyear[,table(mig,gender)]

# Albeit the category "2G migrant (mixed parents)" can be interesting, very few cases
indiv[mig == '2G migrant (mixed parents)',]$mig <- '2G migrant'
indivyear[mig == '2G migrant (mixed parents)',]$mig <- '2G migrant'
indiv[,table(mig,gender)];indivyear[,table(mig,gender)]

# There are very few 2G migrant compared to 1.5G. 
# Let's take those who arrived at very early ages (0 to 5), and reclassify them as 2G
indiv[mig == '1.5G migrant' & ageatmigr %in% 0:5,]$mig <- '2G migrant'
indivyear[mig == '1.5G migrant' & ageatmigr %in% 0:5,]$mig <- '2G migrant'
indiv[,table(mig,gender)];indivyear[,table(mig,gender)]

# Let's reorder the categories
indiv[,mig := factor(mig,levels=c('finnish','migrant','1.5G migrant','2G migrant','2.5G migrant'),
                     labels=c('Native-born','1G migrant','1.5G migrant','2G migrant','2.5G migrant'))]
indivyear[,mig := factor(mig,levels=c('finnish','migrant','1.5G migrant','2G migrant','2.5G migrant'),
                     labels=c('Native-born','1G migrant','1.5G migrant','2G migrant','2.5G migrant'))]

# How do we fare in terms of group sizes?
tiff(filename='fig1.tiff',
     width=25,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=indiv[!(mig %in% c('Native-born'))],
       aes(x=background)) +
  geom_bar(color='black',position='dodge',fill='grey') +
  geom_hline(yintercept = 100,linetype='dashed',color='red') +
  facet_grid(gender~mig,scales='free_y') +
  scale_y_log10() + 
  xlab("") + ylab("Count (logarithmic scale)") + labs(fill='') +
  theme(legend.position = 'none',axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# Let's find a selection of regions for which there is sufficient representation (about 100 men and 100 women per generation)
categories <- list()
for(i in levels(indiv$mig)[-1]){ 
  categories[[i]] <- indiv[gender == 'Man' & mig == i,           
                           as.data.table(table(background))][N >= 100]$background
}
Reduce(intersect,categories) # 6 groups!!

# Let's reduce our migrants to only those in the 6 groups (plus the natives)
indiv <- indiv[background %in% c('Finland',Reduce(intersect,categories))]
indivyear <- indivyear[background %in% c('Finland',Reduce(intersect,categories))]
indiv[,table(mig,gender)];indivyear[,table(mig,gender)]

# How do they do in terms of age
indiv[,table(syntyv,mig,background)] # Estonian and East African migrant-natives might be tricky

# Reorder categories
indiv[,background := factor(background,levels=c('Finland','Russia','Estonia','Central & Eastern Europe',
                                                'Middle East','East Africa','Southeast Asia'))]
indivyear[,background := factor(background,levels=c('Finland','Russia','Estonia','Central & Eastern Europe',
                                                'Middle East','East Africa','Southeast Asia'))]

# Let's showcase the most important nationalities
indiv <- indiv %>%
  mutate(country2 = case_when(
    # Finland
    country == 'Finland' ~ 'Finland',
    # Estonia
    country == 'Estonia' ~ 'Estonia',
    # Russia
    country == 'Russian Federation' ~ 'Russia',
    # Central and Eastern Europe
    country == 'Albania' ~ 'Other',
    country == 'Bosnia and Herzegovina' ~ 'Former Yugoslavia',
    country == 'Croatia' ~ 'Former Yugoslavia',
    country == 'North Macedonia' ~ 'Former Yugoslavia',
    country == 'Montenegro' ~ 'Former Yugoslavia',
    country == 'Serbia' ~ 'Former Yugoslavia',
    country == 'Slovenia' ~ 'Former Yugoslavia',
    country == 'Former Yugoslavia' ~ 'Former Yugoslavia',
    country == 'Former Serbia and Montenegro' ~ 'Former Yugoslavia',
    country == 'Bulgaria' ~ 'Other',
    country == 'Romania' ~ 'Romania',
    country == 'Czechia' ~ 'Other',
    country == 'Hungary' ~ 'Other',
    country == 'Poland' ~ 'Poland',
    country == 'Slovakia' ~ 'Other',
    country == 'Former Czechoslovakia' ~ 'Other',
    country == 'Latvia' ~ 'Other',
    country == 'Lithuania' ~ 'Other',
    country == 'Belarus' ~ 'Other',
    country == 'Moldova' ~ 'Other',
    country == 'Ukraine' ~ 'Ukraine',
    # Middle East
    country == 'Afghanistan' ~ 'Afghanistan',
    country == 'Bahrain' ~ 'Other',
    country == 'Iran' ~ 'Iran',
    country == 'Iraq' ~ 'Iraq',
    country == 'Jordan' ~ 'Other',
    country == 'Kuwait' ~ 'Other',
    country == 'Lebanon' ~ 'Other',
    country == 'Oman' ~ 'Other',
    country == 'Palestine' ~ 'Other',
    country == 'Qatar' ~ 'Other',
    country == 'Saudi Arabia' ~ 'Other',
    country == 'United Arab Emirates' ~ 'Other',
    country == 'Syria' ~ 'Syria',
    country == 'Turkey' ~ 'Turkey',
    country == 'Yemen' ~ 'Other',
    country == 'Former South Yemen' ~ 'Other',
    # East Africa
    country == 'Burundi' ~ 'Other',
    country == 'Comoros' ~ 'Other',
    country == 'Eritrea' ~ 'Eritrea',
    country == 'Ethiopia' ~ 'Ethiopia',
    country == 'Kenya' ~ 'Kenya',
    country == 'Madagascar' ~ 'Other',
    country == 'Mauritius' ~ 'Other',
    country == 'Mayotte' ~ 'Other',
    country == 'Mozambique' ~ 'Other',
    country == 'Reunion' ~ 'Other',
    country == 'Rwanda' ~ 'Other',
    country == 'Seychelles' ~ 'Other',
    country == 'Somalia' ~ 'Somalia',
    country == 'Tanzania' ~ 'Other',
    country == 'Uganda' ~ 'Other',
    # South East Asia
    country == 'Brunei Darussalam' ~ 'Other',
    country == 'Cambodia' ~ 'Other',
    country == 'East Timor' ~ 'Other',
    country == 'Indonesia' ~ 'Other',
    country == 'Malaysia' ~ 'Other',
    country == 'Laos' ~ 'Other',
    country == 'Philippines' ~ 'Philippines',
    country == 'Singapore' ~ 'Other',
    country == 'Thailand' ~ 'Thailand',
    country == 'Vietnam' ~ 'Vietnam'
  ))

indiv[,country2 := factor(country2,
                          levels=c('Finland','Russia','Estonia',
                                   'Poland','Romania','Ukraine','Former Yugoslavia',
                                   'Afghanistan','Iran','Iraq','Syria','Turkey',
                                   'Eritrea','Ethiopia','Kenya','Somalia',
                                   'Philippines','Thailand','Vietnam','Other'))]

# Choose color
clr <- c('bisque2','darkkhaki',
         brewer.pal(5,"Blues")[-1],brewer.pal(6,"Reds")[-1],
         brewer.pal(5,"Greens")[-1],brewer.pal(4,"Purples")[-1],'grey70')

# Visualization
tiff(filename='fig2.tiff',
     width=25,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=indiv[!(mig %in% c('Native-born'))],
       aes(x=background,fill=country2)) +
  geom_bar(color='grey10',position='stack') +
  facet_wrap(gender~mig,scales='free_y',nrow = 2, ncol=4) +
  scale_fill_manual(values = clr) +
  xlab("") + ylab("Count") + labs(fill='') +
  theme(legend.position = 'right',axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

########################################################################################################################

# FERTILITY PATTERNS ----

# Let's start by finding the average number of children per year and group
familysize <- indivyear %>% 
  group_by(background,mig,gender,age) %>% 
  summarize(number_children = mean(children_cum))

tiff(filename='fig3.tiff',
     width=30,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=filter(familysize,age<=31),aes(x=age,y=number_children,color=mig,shape=mig)) +
  geom_point() + geom_line() +
  scale_color_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_shape_manual(values=c(19,15,17,1,4)) +
  facet_grid(gender~background) +
  xlab("Age") + ylab("Number of children (average)") + labs(color=NULL,shape=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

# Now, proportion of individuals who became parents
areparents <- indivyear %>% 
  group_by(background,mig,gender,age) %>% 
  summarize(prop = sum(children_cum != 0)/n())

tiff(filename='fig4.tiff',
     width=30,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=filter(areparents,age<=31),aes(x=age,y=prop,color=mig,shape=mig)) +
  geom_point() + geom_line() +
  scale_color_manual(values=c('red',brewer.pal(5,"Greys")[-1])) +
  scale_shape_manual(values=c(19,15,17,1,4)) +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  facet_grid(gender~background) +
  xlab("Age") + ylab("Parent") + labs(color=NULL,shape=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()
 
########################################################################################################################

# Save image
save(indivyear,file='data3.person.year.RData')

########################################################################################################################