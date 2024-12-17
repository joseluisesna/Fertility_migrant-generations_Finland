########################################################################################################################
# (7) FERTILITY AMONG MIGRANT DESCENDANTS (Analyses, multi-level)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Dec 10th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(ggplot2);library(dplyr);library(tidyverse);library(survival);library(RColorBrewer)
library(zoo);library(scales);library(modelsummary);library(sjPlot);library(performance)
theme_set(theme_bw()) # bw theme

# DATA LOADING ----
load('data5.person.year.RData') # person.year

# OPTIONS (do not use scientific notation)
options(scipen=999)
options(digits=4)

# Function for safe.lapply (in case it does not work for some ancestry groups)
safe.lapply <- function(X,FUN,...){
  lapply(X,function(x){
    result <- try(FUN(x,...),silent=TRUE)
    if(inherits(result,"try-error")){NULL} else{result}
  })
}

########################################################################################################################

# NEW KAPLAN-MEIER CURVES ----

# Extract the curver from matched datasets
km <- km_sum <- alist()
for(i in names(person.year)){
  km[[i]] <- survfit(Surv(start,age,event)~gender+mig,conf.type='log',
                     conf.int=.95,type='kaplan-meier',error='greenwood',
                     data=person.year[[i]])
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
                            background = paste(i),
                            strata = strata_labels)
}

# Now split the strata info
km_sum <- do.call(rbind,km_sum)
km_sum <- km_sum %>%
  separate(strata,into=c('gender','mig'),sep=', ') %>%
  separate(gender,into=c(NA,'gender'),sep='=') %>%
  separate(mig,into=c(NA,'mig'),sep='=') %>%
  mutate(mig = str_trim(mig)) # remove blank spaces

km_sum <- as.data.table(km_sum)

# Change order of appearance of values
km_sum[,background := factor(background,levels=names(person.year))]
km_sum[,mig := ifelse(mig == 'Native-born','Native-born (matched)',mig)]
km_sum[,mig := factor(mig,levels=c('Native-born (matched)','1.5G migrant','2G migrant','2.5G migrant'))]

# Visualization (survival prob)
tiff(filename='fig9.tiff',
     width=30,height=15,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=km_sum %>% filter(age %in% 15:31),
       aes(x=age,y=surv,ymin=lower,ymax=upper,col=mig,fill=mig)) +
  geom_ribbon(alpha=1/3,linetype='dashed',linewidth=.25) +
  geom_line() +
  facet_grid(gender~background) +
  scale_color_manual(values=c('red',brewer.pal(4,"Greys")[2:4])) +
  scale_fill_manual(values=c('red',brewer.pal(4,"Greys")[2:4])) +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  xlab("Age") + ylab("Survival probabilities") + labs(color=NULL,fill=NULL) +
  theme(legend.position = 'top',
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

########################################################################################################################

# VARIABLE PREPARATION ----

# Select variables
for(i in names(person.year)){
  person.year[[i]] <- person.year[[i]][,.(shnro,gender,syntyv,age,mig,event,observed,
                                          region,habitat,family_type,children,family_status,
                                          education,activity,unemployment,marital_status)]
  #person.year[[i]] <- na.omit(person.year[[i]])
}

# Entry and end points
for(i in names(person.year)){
  person.year[[i]][,tstart := age - 15]
  person.year[[i]][,tend := age - 14]
}

# Time-fixed dimensions
for(i in names(person.year)){
  # Migrant generation
  person.year[[i]][,onefive := factor(ifelse(mig == '1.5G migrant',1,0))]
  person.year[[i]][,second := factor(ifelse(mig == '2G migrant',1,0))]
  person.year[[i]][,twofive := factor(ifelse(mig == '2.5G migrant',1,0))]
  # And divided by gender (for additional analyses)
  person.year[[i]][,onefive_m := factor(ifelse(mig == '1.5G migrant' & gender == 'Man',1,0))]
  person.year[[i]][,second_m := factor(ifelse(mig == '2G migrant' & gender == 'Man',1,0))]
  person.year[[i]][,twofive_m := factor(ifelse(mig == '2.5G migrant' & gender == 'Man',1,0))]
  person.year[[i]][,onefive_f := factor(ifelse(mig == '1.5G migrant' & gender == 'Woman',1,0))]
  person.year[[i]][,second_f := factor(ifelse(mig == '2G migrant' & gender == 'Woman',1,0))]
  person.year[[i]][,twofive_f := factor(ifelse(mig == '2.5G migrant' & gender == 'Woman',1,0))]
  # Cohort: 85-89 vs 90-94
  person.year[[i]][,cohort := factor(ifelse(syntyv %in% 1990:94,1,0))]
  # Habitat
  person.year[[i]][,semiurban := factor(ifelse(habitat == 'Semi-urban',1,0))]
  person.year[[i]][,rural := factor(ifelse(habitat == 'Rural',1,0))]
  person.year[[i]][,capital := factor(ifelse(region == 'Uusimaa',1,0))] # a proxy for capital
  # Children
  person.year[[i]][,children.23 := factor(ifelse(children %in% c(2,3),1,0))]
  person.year[[i]][,children.4plus := factor(ifelse(children %in% c(4,5,6,7),1,0))]
  # Family type
  person.year[[i]][,single.parent := factor(ifelse(family_type == 'Single parent',1,0))]
  # Family status
  person.year[[i]][,manual.worker := factor(ifelse(family_status == 'Manual worker',1,0))]
  person.year[[i]][,employee := factor(ifelse(family_status %in% c('Lower-level employee','Upper-level employee'),1,0))]
  person.year[[i]][,lower.employee := factor(ifelse(family_status == 'Lower-level employee',1,0))]
  person.year[[i]][,upper.employee := factor(ifelse(family_status == 'Upper-level employee',1,0))]
  person.year[[i]][,self.employed := factor(ifelse(family_status == 'Self-employed',1,0))]
}

# Let's check we have enough cases
for(i in names(person.year)){
  print(person.year[[i]][,table(family_status,gender)])
}

# Time-varying dimensions 
# Education, although time-varying, will be turn into a time-fixed dimensions (max education reached)
for(i in names(person.year)){
  # first, let's make this variable numeric
  person.year[[i]][,education_num := ifelse(education == 'Less than secondary or unknown',1,
                                            ifelse(education == "Bachelor's, master's, or doctoral",3,2))]
  # max education reached
  person.year[[i]] <- person.year[[i]] %>%
    group_by(shnro) %>%
    mutate(education_max = max(education_num,na.rm=TRUE)) %>%
    ungroup()
  person.year[[i]] <- as.data.table(person.year[[i]])
  person.year[[i]][,education_max := as.factor(education_max)]
}

for(i in names(person.year)){
  # Proxy for high education and mid-level education
  person.year[[i]][,mid.edu := factor(ifelse(education_max == 2,1,0))]
  person.year[[i]][,high.edu := factor(ifelse(education_max == 3,1,0))]
  # Proxy for if employed
  person.year[[i]][,employed := factor(ifelse(activity %in% c('Employed'),1,0))]
}

# Let's check we have enough cases
for(i in names(person.year)){
  print(person.year[[i]][!duplicated(shnro),table(education_max,gender,useNA = 'always')])
}

# Months of unemployment in the prior 3 years
for(i in names(person.year)){
  person.year[[i]] <- person.year[[i]] %>%
    arrange(shnro,age) %>%
    group_by(shnro) %>%
    mutate(
      cum_unemployment = rollapply(
        unemployment,
        width = 3, # sum up the last 3 years
        FUN = sum,
        align = 'right',
        fill = 0 # fill with zero if the window is incomplete
      )    
    ) %>%
    ungroup()
  person.year[[i]] <- as.data.table(person.year[[i]])
  # Also its square-root version
  person.year[[i]][,cum_unemployment_sqrt := sqrt(cum_unemployment)]
}

# Marital status (also as with education, if ever formed a couple)
for(i in names(person.year)){
  # first, let's make this variable numeric
  person.year[[i]][,partner := ifelse(marital_status == 'Married, in partnership, separated',1,0)]
  # find if ever formed a couple...
  person.year[[i]] <- person.year[[i]] %>%
    group_by(shnro) %>%
    mutate(formed_couple = max(partner,na.rm = TRUE)) %>%
    ungroup()
  person.year[[i]] <- as.data.table(person.year[[i]])
  person.year[[i]][is.infinite(formed_couple)]$formed_couple <- NA # if not a number, set to missing
  person.year[[i]][,partner := as.factor(partner)]
  person.year[[i]][,formed_couple := as.factor(formed_couple)]
}

# Let's check we have enough cases
for(i in names(person.year)){
  print(person.year[[i]][!duplicated(shnro),table(formed_couple,gender,useNA = 'always')])
}

########################################################################################################################

# SUMMARY OF DATA (UNIVARIATE) ----

# Select variables 
for(i in names(person.year)){
  person.year[[i]] <- person.year[[i]][,.(shnro,tstart,tend,event,gender,observed,
                                          onefive,second,twofive,
                                          onefive_m,onefive_f,second_m,second_f,twofive_m,twofive_f,
                                          cohort,capital,semiurban,rural,children.23,children.4plus,single.parent,
                                          manual.worker,lower.employee,upper.employee,self.employed,
                                          education_max,mid.edu,high.edu,
                                          employed,cum_unemployment_sqrt,
                                          formed_couple,partner)]
}

# Print summary
for(i in names(person.year)){
  datasummary_skim(person.year[[i]][,-c('shnro','onefive_m','onefive_f','second_m','second_f','twofive_m','twofive_f',
                                        'education_max')],
                   output=paste('table.univar',i,'html',sep='.'))
}

# SUMMARY (BIVARIATE) ----
corr_data <- person.year

for(i in names(corr_data)){
  corr_data[[i]] <- corr_data[[i]][,-c('shnro','tstart','tend',
                                       'onefive_m','onefive_f','second_m','second_f','twofive_m','twofive_f',
                                       'education_max')]
  corr_data[[i]] <- corr_data[[i]][,gender := factor(ifelse(gender == 'Woman',1,0))] # make gender numeric
  # Make all variables numeric for this
  corr_data[[i]] <- as.data.table(lapply(corr_data[[i]],function(x){as.numeric(as.character(x))}))
  datasummary_correlation(corr_data[[i]],
                          output=paste('table.bivar',i,'html',sep='.'))
}

########################################################################################################################

# INTERCATION TERMS WITH TIME ----

for(i in names(person.year)){
  # I will use the logarithm of time
  person.year[[i]][,log.time := log(tend)]
  # Interaction terms
  person.year[[i]][,high.edu_t := (as.integer(high.edu)-1)*log(tend)]
  person.year[[i]][,mid.edu_t := (as.integer(mid.edu)-1)*log(tend)]
  person.year[[i]][,employed_t := (as.integer(employed)-1)*log(tend)]
  person.year[[i]][,cum_unemployment_sqrt_t := cum_unemployment_sqrt*log(tend)]
  person.year[[i]][,partner_t := (as.integer(partner)-1)*log(tend)]
  person.year[[i]][,formed_couple_t := (as.integer(formed_couple)-1)*log(tend)]
}

########################################################################################################################

# COX REGRESSION MODELS ----

# To allocate results
models1 <- models2 <- models3 <- models4 <- models5 <- models6 <- models7 <- models8 <- alist()

# Model 1
formula1 <- "Surv(tstart,tend,event) ~ strata(gender) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             frailty(shnro)"

for(i in names(person.year)){
  models1[[i]] <- coxph(formula = as.formula(formula1), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models1,cox.zph)
# check for collinearity
safe.lapply(models1,check_collinearity)
# Print models
#sjPlot::tab_model(models1,dv.labels=names(person.year),show.ci=FALSE,file='models1.docx')

# Model 2
formula2 <- "Surv(tstart,tend,event) ~ strata(gender) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             capital + semiurban + rural + 
             frailty(shnro)"

for(i in names(person.year)){
  models2[[i]] <- coxph(formula = as.formula(formula2), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models2,cox.zph)
# check for collinearity
safe.lapply(models2,check_collinearity)
# Print models
#sjPlot::tab_model(models2,dv.labels=names(person.year),show.ci=FALSE,file='models2.docx')

# Model 3
# Stratification by gender AND max education level reached (otherwise, problem of non-proportionality)
formula3 <- "Surv(tstart,tend,event) ~ strata(interaction(gender,education_max)) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             frailty(shnro)"

for(i in names(person.year)){
  models3[[i]] <- coxph(formula = as.formula(formula3), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models3,cox.zph)
# check for collinearity
safe.lapply(models3,check_collinearity)
# Print models
#sjPlot::tab_model(models3,dv.labels=names(person.year),show.ci=FALSE,file='models3.docx')

# Model 4
formula4 <- "Surv(tstart,tend,event) ~ strata(gender) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             employed + cum_unemployment_sqrt + employed_t + cum_unemployment_sqrt_t +
             frailty(shnro)"

for(i in names(person.year)){
  models4[[i]] <- coxph(formula = as.formula(formula4), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models4,cox.zph)
# check for collinearity
safe.lapply(models4,check_collinearity)
# Print models
#sjPlot::tab_model(models4,dv.labels=names(person.year),show.ci=FALSE,file='models4.docx')

# Model 5
formula5 <- "Surv(tstart,tend,event) ~ strata(gender) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             partner + partner_t +
             frailty(shnro)"

for(i in names(person.year)){
  models5[[i]] <- coxph(formula = as.formula(formula5), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models5,cox.zph)
# check for collinearity
safe.lapply(models5,check_collinearity)
# Print models
#sjPlot::tab_model(models5,dv.labels=names(person.year),show.ci=FALSE,file='models5.docx')

# Model 6 (alternative to partnership status)
formula6 <- "Surv(tstart,tend,event) ~ strata(interaction(gender,formed_couple)) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             frailty(shnro)"

for(i in names(person.year)){
  models6[[i]] <- coxph(formula = as.formula(formula6), 
                        data=person.year[[i]],
                        method='efron')
}# Check proportionality assumption
safe.lapply(models6,cox.zph)
# check for collinearity
safe.lapply(models6,check_collinearity)
# Print models
#sjPlot::tab_model(models6,dv.labels=names(person.year),show.ci=FALSE,file='models6.docx')

# All
formula7 <- "Surv(tstart,tend,event) ~ strata(interaction(gender,education_max)) +
             onefive + second + twofive + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             capital + semiurban + rural + 
             employed + cum_unemployment_sqrt + employed_t + cum_unemployment_sqrt_t + 
             partner + partner_t +
             frailty(shnro)"

for(i in names(person.year)){
  models7[[i]] <- coxph(formula = as.formula(formula7), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models7,cox.zph)
# check for collinearity
safe.lapply(models7,check_collinearity)
# Print models
#sjPlot::tab_model(models7,dv.labels=names(person.year),show.ci=FALSE,file='models7.docx')

# Model 8 (separate generation effects by gender)
formula8 <- "Surv(tstart,tend,event) ~ strata(interaction(gender)) +
             onefive_m + onefive_f + second_m + second_f + twofive_m + twofive_f + cohort +
             children.23 + children.4plus + single.parent + 
             manual.worker + lower.employee + upper.employee + self.employed +
             frailty(shnro)"

for(i in names(person.year)){
  models8[[i]] <- coxph(formula = as.formula(formula8), 
                        data=person.year[[i]],
                        method='efron')
}
# Check proportionality assumption
safe.lapply(models8,cox.zph)
# check for collinearity
safe.lapply(models8,check_collinearity)
# Print models
#sjPlot::tab_model(models8,dv.labels=names(person.year),show.ci=FALSE,file='models8.docx')

########################################################################################################################

# MODEL COMPARISON ----

for(i in names(models1)){
  print(abs(-2*logLik(models2[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models3[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models4[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models5[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models6[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models7[[i]]) - -2*logLik(models1[[i]])))
  print(abs(-2*logLik(models8[[i]]) - -2*logLik(models1[[i]])))
  #print(anova(models1[[i]],models2[[i]],models3[[i]],models4[[i]]))
}

# Goodness of fit (AIC)
data.table(background = names(models1),
           m1_AIC = unlist(safe.lapply(models1,AIC)),
           m2_AIC = unlist(safe.lapply(models2,AIC)),
           m3_AIC = unlist(safe.lapply(models3,AIC)),
           m4_AIC = unlist(safe.lapply(models4,AIC)),
           m5_AIC = unlist(safe.lapply(models5,AIC)),
           m6_AIC = unlist(safe.lapply(models6,AIC)),
           m7_AIC = unlist(safe.lapply(models7,AIC)),
           m8_AIC = unlist(safe.lapply(models8,AIC)))
# Goodness of fit (BIC)
data.table(background = names(models1),
           m1_BIC = unlist(safe.lapply(models1,BIC)),
           m2_BIC = unlist(safe.lapply(models2,BIC)),
           m3_BIC = unlist(safe.lapply(models3,BIC)),
           m4_BIC = unlist(safe.lapply(models4,BIC)),
           m5_BIC = unlist(safe.lapply(models5,BIC)),
           m6_BIC = unlist(safe.lapply(models6,BIC)),
           m7_BIC = unlist(safe.lapply(models7,BIC)),
           m8_BIC = unlist(safe.lapply(models8,BIC)))

########################################################################################################################

# SUMMARY OF RESULTS (VISUALIZATION) ----

# Let's draw the coefficients
# Let's write function to extra the coefficients
extract.info <- function(model){
  sum_coeff <- summary(model)
  rmrandomterm <- nrow(sum_coeff$coefficients)
  output <- data.table(mig = rownames(sum_coeff$coefficients)[1:rmrandomterm],
                       hr = sum_coeff$conf.int[,'exp(coef)'],
                       cilow = sum_coeff$conf.int[,'lower .95'],
                       ciup = sum_coeff$conf.int[,'upper .95'],
                       pval = sum_coeff$coefficients[,'p'][1:rmrandomterm])
  return(output)
}

# Extract results
results <- vector('list',length=7)

for(i in seq_along(results)){
  # This call
  obj <- get(paste('models',i,sep=''))
  for(j in names(obj)){
    results[[i]][[j]] <- extract.info(obj[[j]])
    results[[i]][[j]]$model <- i
    results[[i]][[j]]$background <- j
  }
  results[[i]] <- do.call(rbind,results[[i]])
}
results <- do.call(rbind,results)

# Let's keep only effects of interest
results <- results[mig %in% c('onefive1','second1','twofive1')]
# Change order of appearance of values
results[,mig := factor(mig,labels = c('1.5G migrant','2G migrant','2.5G migrant'))]
results[,background := factor(background,levels=names(person.year))]
results[,model := factor(model)]
results[,sig := factor(ifelse(pval < 0.05 & hr > 1,1,ifelse(pval < 0.05 & hr < 1,2,3)))]

# Visualization
tiff(filename='fig10ext.tiff',
     width=32,height=12,units='cm',compression='lzw',bg='white',res=1000)
ggplot(data=results[model != 6], # remember this is an alternative specification
       aes(x=model,y=hr,ymin=cilow,ymax=ciup,shape=model,color=sig)) +
  geom_hline(yintercept = 1,linetype='dashed') +
  geom_pointrange(position=position_dodge(width=0.75)) +
  facet_grid(mig~background) +
  scale_shape_manual(values=c(1,5,7,8,2,19),
                     labels=c(expression("Model 1 (Baseline)"),
                              expression("Model 2 (Habitat)"),
                              expression("Model 3 (Education (strata))"),
                              expression("Model 4 (Employment"[italic(j)] ~ ")"),
                              expression("Model 5 (Marital status"[italic(j)] ~ ")"),
                              expression("Model 6 (Full)"))) +
  scale_color_manual(values=c('blue2','red2','grey50'),
                     labels=c(expression(italic(p) < 0.05),
                              expression(italic(p) < 0.05),
                              expression(italic(p) >= 0.05))) +
  scale_y_log10() +
  labs(x=NULL,y='Hazard ratio',color=NULL,shape=NULL) +
  theme(legend.position = 'right',axis.text.x = element_blank(),
        strip.background = element_rect(fill='black'),strip.text=element_text(color='white'))
dev.off()

########################################################################################################################