########################################################################################################################
# (1) FERTILITY AMONG MIGRANT DESCENDANTS (Data retrieval)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Nov 6th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(tidyverse)

########################################################################################################################

# FOLK DATA (1999-2019) ----
data2 <- data1 <- list()

for(i in 1999:2009){
  data1[[as.character(i)]] <- data.table(readRDS(paste("W:/Nordforsk/shared/data/Basic_yearwise_info/basic_",i,"_all_ids.RDS",sep='')))
}
# Sample selection
for(i in seq_along(data1)){
  # birth cohort 1980-94, only observations when they are 15 to 34 years of age
  data1[[i]] <- data1[[i]][(syntyv %in% 1985:1994) & ika %in% 14:34] # keep 14 for data lagging purposes
}

for(i in 2010:2019){
  data2[[as.character(i)]] <- data.table(readRDS(paste("W:/Nordforsk/shared/data/Basic_yearwise_info/basic_",i,"_all_ids.RDS",sep='')))
}
# Sample selection
for(i in seq_along(data2)){
  data2[[i]] <- data2[[i]][(syntyv %in% 1985:1994) & ika %in% 14:34]
}

# Put data altogether
person.year <- append(data1,data2)

# Count individuals rather than person-year observations
indiv <- do.call(rbind,person.year)
indiv <- indiv[order(vuosi,decreasing = TRUE)][!duplicated(shnro)] # pick only the last observation of each person

# Save object
save(person.year,file='data.person.year.RData')
save(indiv,file='data.sample.RData') # 742,966 individuals (our focal subjects)
rm(list=setdiff(ls(),'indiv'))

########################################################################################################################

# COHABITATION DATA (1987-2019) ----

# Data loading
cohab <- list()
cohab[['1']] <- fread("D:/ready-made/FOLK_aslii_8800a/folk_19872000_tua_aslii21tot_1.csv")
cohab[['2']] <- fread("D:/ready-made/FOLK_aslii_0110a/folk_20012010_tua_aslii21tot_1.csv")
cohab[['3']] <- fread("D:/ready-made/FOLK_aslii_11a/folk_20112020_tua_aslii21tot_1.csv")
cohab <- do.call(rbind,cohab)

# Removed data recorded before 2000 (a year earlier for lagging later) or not pertaining our focal subjects
cohab <- cohab[vuosi >= 1999]
cohab <- cohab[shnro %in% indiv$shnro | spuhnro %in% indiv$shnro]
# Many observations are redundant. Let's keep the last observation per unique couple
cohab <- cohab %>%
  mutate(IDcouple := paste(shnro,spuhnro,sep='-')) %>% # create a couple ID
  arrange(desc(vuosi)) %>%
  filter(!duplicated(IDcouple)) %>% # keep only the last observation per couple
  select(-IDcouple) # remove this couple label

# Find the start and end dates for every couple
cohab <- cohab %>%
  mutate(startdate = as.IDate(pmin(alku,ymuuttopv,vihkipvm,na.rm=TRUE))) %>% # marrying or moving together (whatever happened first)
  mutate(enddate = as.IDate(pmin(loppu,emuuttopv,asumero,kuolpv,pu_kuolpv,na.rm=TRUE)))

# All partnerships registered in Finland between 1999 and 2020
save(cohab,file='data.unions.RData')

########################################################################################################################

# PARENT-CHILD & MOBILITY DATA ----

# We need IDS of the partners of those in the sample to reconstruct their migratory background too
partnersid <- cohab[!(shnro %in% indiv[,shnro]),shnro]
connecttoparents <- c(indiv[,shnro],partnersid) # all those whose background we need
connecttoparents <- unique(connecttoparents) # 924,824 (either focal subjects or their (former) partners)

# We need to connect focal subjects and their partners to their parents to reconstruct their migratory background
parents <- fread("W:/Nordforsk/shared/data/u1708_a_folkllvv_19872021_f.csv")
muutt <- fread("W:/Nordforsk/shared/data/u1708_a_folkmuutto_1987_2021.csv")

# PARENTHOOD TRANSITIONS ----
children <- parents

# For parenthood transition, let's keep those children with a parent among the focal subjects (or their partners)
children <- children[shnro_m %in% connecttoparents | shnro_f %in% connecttoparents]
children[,length(unique(shnro))] # 355,622 children

# What we need are not children themselves but parenthood transitions
children <- children[,.(shnro,shnro_m,shnro_f,syntyv,vuosi)]
names(children) <- c('child','mom','dad','syntyv','vuosi')
children <- pivot_longer(children,cols=c('mom','dad'),
                              names_to='parent',values_to = 'shnro')
# Remove when Id of parent is missing
children <- children %>% filter(shnro != '')

# Now, for every parent-child, keep the earliest observation
children <- children %>%
  mutate(IDparchi = paste(shnro,child,sep='-')) %>%
  arrange(vuosi) %>%
  filter(!duplicated(IDparchi)) %>%
  select(c(shnro,child,syntyv))
children <- as.data.table(children)

# How many birth did our focal subjects had
children %>% 
  filter(shnro %in% indiv$shnro) %>% 
  summarize(n()) # 447,202 births
children %>% 
  filter(shnro %in% indiv$shnro) %>%
  filter(!duplicated(child)) %>% 
  summarize(n()) # 298,483
                
# Save image
save(children,file='data.children.RData')

# BACKGROUND RECONSTRUCTION ----

# Now, let's connect our sample (and their partners) to their parents
forbackg <- parents[shnro %in% connecttoparents]
forbackg[,length(unique(shnro))] # these are 809,153 individuals of the 924,824 we need parents' info
# How many unique parents are here?
forbackg[,length(unique(c(shnro_m,shnro_f)))] # 982,373

# Let's correct blank spaces "" with NAs
forbackg[shnro_m == "",]$shnro_m <- NA
forbackg[shnro_f == "",]$shnro_f <- NA

# Let's simplify these data
b1 <- forbackg[,.(shnro,vuosi,shnro_m)]
b2 <- forbackg[,.(shnro,vuosi,shnro_f)]
# Remove NAs
b1 <- na.omit(b1)
b2 <- na.omit(b2)
# Let's take the first observation
b1 <- b1[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,shnro_m)]
b2 <- b2[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,shnro_f)]
sample.parents <- merge(b1,b2,by='shnro',all=TRUE)

# Finally let's add all those for whom no info on parents is available
noinfoparents <- connecttoparents[!(connecttoparents %in% sample.parents[,shnro])]
sample.parents <- rbind(sample.parents,
                        data.table(shnro = noinfoparents,shnro_m=NA,shnro_f=NA))

# Save image
save(sample.parents,file='data.sample.partners.parents.RData')

########################################################################################################################

# MOBILITY DATA ----

# Now, let's identify in the mobility data traces of international migration
# sample of people whose mobility we need to identify
IDsall <- unique(c(sample.parents[,shnro],
                   sample.parents[,shnro_m],
                   sample.parents[,shnro_f])) # focal subjects, plus their partners and parents
muutt <- muutt[shnro %in% IDsall] # Remove data pertaining to other people different from those we need

# Remember that values 991, 997, 998 and 999 are forms of missing, so let's re-categorize them
muutt[kansa1_m %in% 991:999]$kansa1_m <- NA # nationality at migration
muutt[kansa2_m %in% 991:999]$kansa2_m <- NA # second nationality
muutt[svaltio_m %in% 991:999]$svaltio_m <- NA # country of birth
muutt[lahtomaakoodi %in% 991:999]$lahtomaakoodi <- NA # country of departure
muutt[kieli_m %in% 98:99]$kieli_m <- NA # language

# We can retrieve the migratory background: first and second nationality, country of birth, language (for former USSR)
bg1 <- muutt[,.(shnro,vuosi,kansa1_m)] # first nationality
bg2 <- muutt[,.(shnro,vuosi,kansa2_m)] # second nationality
bg3 <- muutt[,.(shnro,vuosi,svaltio_m)] # country of birth
bg4 <- muutt[,.(shnro,vuosi,lahtomaakoodi)] # country of departure 
bg5 <- muutt[,.(shnro,vuosi,kieli_m)] # language 
# Remove NAs
bg1 <- na.omit(bg1)
bg2 <- na.omit(bg2)
bg3 <- na.omit(bg3)
bg4 <- na.omit(bg4)
bg5 <- na.omit(bg5)

# Take the oldest observation per person for nationality (sometimes people adopt the Finnish nationality over time) and language
bg1 <- bg1[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,kansa1_m)]
bg2 <- bg2[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,kansa2_m)]
bg4 <- bg4[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,lahtomaakoodi)]
bg5 <- bg5[order(vuosi,decreasing = FALSE)][!duplicated(shnro)][,.(shnro,kieli_m)]
# But the newest observation for country of birth (sometimes registers correct this information)
bg3 <- bg3[order(vuosi,decreasing = TRUE)][!duplicated(shnro)][,.(shnro,svaltio_m)]

# Let's put all the information together
IDsall <- data.table(shnro = IDsall) # ids as datatable instead of vector
IDsall <- merge(IDsall,bg1,by='shnro',all.x = TRUE)
IDsall <- merge(IDsall,bg2,by='shnro',all.x = TRUE)
IDsall <- merge(IDsall,bg3,by='shnro',all.x = TRUE)
IDsall <- merge(IDsall,bg4,by='shnro',all.x = TRUE)
IDsall <- merge(IDsall,bg5,by='shnro',all.x = TRUE)

# Remove first row (NA)
IDsall <- IDsall[!is.na(shnro)]

# Let's add the additional info in muutt: date of migration, and marital status
add1 <- muutt[,.(shnro,muuttopv,sivs_m)]
add1 <- add1[order(muuttopv,decreasing = FALSE)][!duplicated(shnro)] # oldest observation only
IDsall <- merge(IDsall,add1,by='shnro',all.x=TRUE)

# rename variables
names(IDsall) <- c('shnro','nation1','nation2','birthcountry','departurecountry','language','migrationdate','maritalstatus') 
backgrounds <- IDsall

# Delete unnecessary objects
rm(list=setdiff(ls(),'backgrounds'))

# Finally, let's complement the information in MUTTo with the FOLK data
data3 <- data2 <- data1 <- list()

for(i in 1987:1997){
  data1[[as.character(i)]] <- data.table(readRDS(paste("W:/Nordforsk/shared/data/Basic_yearwise_info/basic_",i,"_all_ids.RDS",sep='')))
}
# Sample selection
for(i in seq_along(data1)){
  data1[[i]] <- data1[[i]][shnro %in% backgrounds[,shnro],
                           .(shnro,vuosi,syntyv,syntyp2)] # keep ID, year of birth, and migrant status
}

for(i in 1998:2008){
  data2[[as.character(i)]] <- data.table(readRDS(paste("W:/Nordforsk/shared/data/Basic_yearwise_info/basic_",i,"_all_ids.RDS",sep='')))
}
# Sample selection
for(i in seq_along(data2)){
  # only observations from parents and partners (those in backgrounds2)
  data2[[i]] <- data2[[i]][shnro %in% backgrounds[,shnro],
                           .(shnro,vuosi,syntyv,syntyp2)]
}

for(i in 2009:2019){
  data3[[as.character(i)]] <- data.table(readRDS(paste("W:/Nordforsk/shared/data/Basic_yearwise_info/basic_",i,"_all_ids.RDS",sep='')))
}
# Sample selection
for(i in seq_along(data3)){
  data3[[i]] <- data3[[i]][shnro %in% backgrounds[,shnro],
                           .(shnro,vuosi,syntyv,syntyp2)]
}

# Put data altogether
data <- append(data1,append(data2,data3))
rm(data1);rm(data2);rm(data3)
data <- do.call(rbind,data) # omit sub-folders
data <- na.omit(data) # remove in case of missing data
# Remove observations from the same individuals, keep only the first one
data <- data[order(vuosi,decreasing = FALSE)][!duplicated(shnro)]

# Now, let's merge the data
backgrounds <- merge(backgrounds,data[,.(shnro,syntyp2,syntyv)],by='shnro',all.x=TRUE)

# Save image
save(backgrounds,file='data.backgrounds.RData')

########################################################################################################################