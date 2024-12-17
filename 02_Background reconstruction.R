########################################################################################################################
# (2) FERTILITY AMONG MIGRANT DESCENDANTS (Background reconstruction)
# R script written by Jose Luis Estevez (University of Helsinki)
# Date: Nov 7th, 2024
########################################################################################################################

# CLEAN THE ENVIRONMENT ----
rm(list=ls())

# REQUIRED PACKAGES ----
library(data.table);library(dplyr)

# DATA LOADING ----
load('data.backgrounds.RData')
load("data.sample.partners.parents.RData")

########################################################################################################################

# FINNISH VS. MIGRANT (OR MIGRANT DESCENDANT) ----
# Divide Finnish, migrants or 2G migrants
backgrounds[,mig := ifelse(!is.na(syntyp2) & syntyp2 %in% c(11,12),'finnish',
                           ifelse(!is.na(syntyp2) & syntyp2 == 21,'2G migrant',
                                  ifelse(!is.na(syntyp2) & syntyp2 == 22,'migrant',NA)))]
backgrounds[,table(mig,useNA='always')] 

# WHICH COUNTRY?
backgrounds[mig == 'finnish',coe := 246] # Finland for the Finnish
# For migrants, look at first and second nationality, country of birth and departure
backgrounds[mig == 'migrant',
            coe := ifelse(!is.na(nation1) & nation1 != 246,nation1,
                          ifelse(!is.na(nation2) & nation2 != 246,nation2,
                                 ifelse(!is.na(birthcountry) & birthcountry != 246,birthcountry,
                                        ifelse(!is.na(departurecountry) & departurecountry != 246,departurecountry,NA))))]
backgrounds[,table(coe == 246,mig,useNA = 'always')]

# Let' add this information to the sample [remember many of those missing cases are probably parents]
sample <- sample.parents
sample <- merge(sample,backgrounds,by='shnro',all.x=TRUE)
sample[,table(mig,useNA='always')] 
sample[,table(mig,coe==246,useNA='always')] # we only miss a coe of 260 migrants

# PARENTS' DATA FOR MIGRANT DESCENDANTS (FOREIGN ASCENDANCY) ----

# Now let's add information about parents, if available, to find people with foreign ascendancy too
backgrounds <- backgrounds[,.(shnro,mig,coe)]
names(backgrounds) <- c('shnro_m','mig_m','coe_m') # mother's information
sample <- merge(sample,backgrounds,by='shnro_m',all.x=TRUE)
names(backgrounds) <- c('shnro_f','mig_f','coe_f') # father's information
sample <- merge(sample,backgrounds,by='shnro_f',all.x=TRUE)

# Using this information, we can detect individuals with mixed backgrounds
# There are 4 cases with both parents foreigners
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant'] 
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant']$coe 
# By looking at the parents, we can infer their backgrounds
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant']$coe_m
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant']$coe_f
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant']$coe <- c(804,643,NA,233)
# Change their status to 2G migrants
sample[mig == 'finnish' & mig_m == 'migrant' & mig_f == 'migrant']$mig <- '2G migrant'

# THE 2.5 GENERATION MIGRANTS ----
# Those with at least one migrant parent, re-categorized as 2.5G migrants
sample[mig == 'finnish' & (mig_m == 'migrant' | mig_f == 'migrant')]$mig <- '2.5G migrant' 
sample[,table(mig,syntyp2,useNA='always')]

# For 2.5G migrant we consider their country of ethnicity, the country of ethnicity of the migrant parent
sample[mig == '2.5G migrant',coe := ifelse(!is.na(coe_m) & coe_m != 246,coe_m,
                                           ifelse(!is.na(coe_f) & coe_f != 246,coe_f,NA))]
sample[,table(mig,coe==246,useNA='always')] # We retrieve 14,673, but miss still 1,110

# THE 2ND GENERATION MIGRANTS ----
# First, check whether both parents share the same origin or not
sample[mig == '2G migrant']
sample[mig == '2G migrant' & coe_m == coe_f] # In 2,590/3,436 cases (75.4%) both parents have the same origin
sample[mig == '2G migrant' & coe_m != coe_f] # In 506 cases (14.7%), they come from different countries
sample[mig == '2G migrant' & coe_m != coe_f]$mig <- '2G migrant (mixed parents)' # Let's single these out

# We will give priority to the mother's origin
sample[mig %in% c('2G migrant','2G migrant (mixed parents)'),
       coe := ifelse(!is.na(coe_m) & coe_m != 246,coe_m,
                     ifelse(!is.na(coe_f) & coe_f != 246,coe_f,NA))]

sample[,table(mig,coe==246,useNA='always')] # We miss the origin of 250 of these 2G migrants

# THE MIGRANTS ----
# For migrants with uncertain origin, we can rely on their parents too
sample[mig == 'migrant' & is.na(coe) & (!is.na(coe_m) | !is.na(coe_f))]

sample[mig == 'migrant' & is.na(coe),
       coe := ifelse(!is.na(coe_m) & coe_m != 246,coe_m,
                     ifelse(!is.na(coe_f) & coe_f != 246,coe_f,NA))]

sample[,table(mig,coe==246,useNA='always')] # we add the coe to an additional cases using parents' info

# IT REMAINS TO DIFFERENTIATE THE 1.5 GENERATION BY CHECKING DATE OF MIGRATION ----
sample[mig == 'migrant',migrationyear := substr(migrationdate,1,4)] # year of migration
sample[mig == 'migrant',ageatmigr := as.numeric(migrationyear) - as.numeric(syntyv)]
sample[mig == 'migrant' & ageatmigr <= 15,]$mig <- '1.5G migrant' # If arrived with 15 or less

# See the full sample
sample[,table(mig,coe==246,useNA='always')]

########################################################################################################################

# COUNTRY NAMES INSTEAD OF CODES ----

# Increase memory for this procedure
Sys.setenv(R_MAX_VSIZE=1000)

sample <- sample %>%
  mutate(country = case_when(
    coe == '4' ~ 'Afghanistan',
    coe == '8' ~ 'Albania',
    coe == '12' ~ 'Algeria',
    coe == '16' ~ 'American Samoa',
    coe == '20' ~ 'Andorra',
    coe == '24' ~ 'Angola',
    coe == '660' ~ 'Anguilla',
    coe == '10' ~ 'Antarctica',
    coe == '28' ~ 'Antigua and Barbuda',
    coe == '32' ~ 'Argentina',
    coe == '51' ~ 'Armenia',
    coe == '533' ~ 'Aruba',
    coe == '36' ~ 'Australia',
    coe == '40' ~ 'Austria',
    coe == '31' ~ 'Azerbaijan',
    coe == '44' ~ 'Bahamas',
    coe == '48' ~ 'Bahrain',
    coe == '50' ~ 'Bangladesh',
    coe == '52' ~ 'Barbados',
    coe == '112' ~ 'Belarus',
    coe == '56' ~ 'Belgium',
    coe == '84' ~ 'Belize',
    coe == '204' ~ 'Benin',
    coe == '60' ~ 'Bermuda',
    coe == '64' ~ 'Bhutan',
    coe == '68' ~ 'Bolivia',
    coe == '535' ~ 'Bonaire, Sint Eustatius and Saba',
    coe == '70' ~ 'Bosnia and Herzegovina',
    coe == '72' ~ 'Botswana',
    coe == '74' ~ 'Bouvet Island',
    coe == '76' ~ 'Brazil',
    coe == '86' ~ 'British Indian Ocean Territory',
    coe == '96' ~ 'Brunei Darussalam',
    coe == '100' ~ 'Bulgaria',
    coe == '854' ~ 'Burkina Faso',
    coe == '108' ~ 'Burundi',
    coe == '116' ~ 'Cambodia',
    coe == '120' ~ 'Cameroon',
    coe == '124' ~ 'Canada',
    coe == '132' ~ 'Cabo Verde',
    coe == '136' ~ 'Cayman Islands',
    coe == '140' ~ 'Central African Republic',
    coe == '148' ~ 'Chad',
    coe == '152' ~ 'Chile',
    coe == '156' ~ 'China',
    coe == '162' ~ 'Christmas Island',
    coe == '166' ~ 'Cocos (Keeling) Islands',
    coe == '170' ~ 'Colombia',
    coe == '174' ~ 'Comoros',
    coe == '178' ~ 'Congo (Congo-Brazzaville)',
    coe == '180' ~ 'Congo, The Democratic Republic of',
    coe == '184' ~ 'Cook Islands',
    coe == '188' ~ 'Costa Rica',
    coe == '384' ~ 'Ivory Coast',
    coe == '191' ~ 'Croatia',
    coe == '192' ~ 'Cuba',
    coe == '531' ~ 'Curacao',
    coe == '196' ~ 'Cyprus',
    coe == '203' ~ 'Czechia',
    coe == '208' ~ 'Denmark',
    coe == '262' ~ 'Djibouti',
    coe == '212' ~ 'Dominica',
    coe == '214' ~ 'Dominican Republic',
    coe == '626' ~ 'East Timor',
    coe == '218' ~ 'Ecuador',
    coe == '818' ~ 'Egypt',
    coe == '222' ~ 'El Salvador',
    coe == '226' ~ 'Equatorial Guinea',
    coe == '232' ~ 'Eritrea',
    coe == '233' ~ 'Estonia',
    coe == '231' ~ 'Ethiopia',
    coe == '238' ~ 'Falkland Islands (Malvinas)',
    coe == '234' ~ 'Faroe Islands',
    coe == '242' ~ 'Fiji',
    coe == '246' ~ 'Finland',
    coe == '250' ~ 'France',
    coe == '254' ~ 'French Guiana',
    coe == '258' ~ 'French Polynesia',
    coe == '260' ~ 'French Southern Territories',
    coe == '266' ~ 'Gabon',
    coe == '270' ~ 'Gambia',
    coe == '268' ~ 'Georgia',
    coe == '276' ~ 'Germany',
    coe == '288' ~ 'Ghana',
    coe == '292' ~ 'Gibraltar',
    coe == '300' ~ 'Greece',
    coe == '304' ~ 'Greenland',
    coe == '308' ~ 'Grenada',
    coe == '312' ~ 'Guadeloupe',
    coe == '316' ~ 'Guam',
    coe == '320' ~ 'Guatemala',
    coe == '831' ~ 'Guernsey',
    coe == '324' ~ 'Guinea',
    coe == '624' ~ 'Guinea-Bissau',
    coe == '328' ~ 'Guyana',
    coe == '332' ~ 'Haiti',
    coe == '334' ~ 'Heard Island and McDonald Island',
    coe == '336' ~ 'Holy See (Vatican City State)',
    coe == '340' ~ 'Honduras',
    coe == '344' ~ 'Hong Kong',
    coe == '348' ~ 'Hungary',
    coe == '352' ~ 'Iceland',
    coe == '356' ~ 'India',
    coe == '360' ~ 'Indonesia',
    coe == '364' ~ 'Iran',
    coe == '368' ~ 'Iraq',
    coe == '372' ~ 'Ireland',
    coe == '833' ~ 'Isle of Man',
    coe == '376' ~ 'Israel',
    coe == '380' ~ 'Italy',
    coe == '388' ~ 'Jamaica',
    coe == '392' ~ 'Japan',
    coe == '832' ~ 'Jersey',
    coe == '400' ~ 'Jordan',
    coe == '398' ~ 'Kazakhstan',
    coe == '404' ~ 'Kenya',
    coe == '296' ~ 'Kiribati',
    coe == '408' ~ 'North Korea',
    coe == '410' ~ 'South Korea',
    coe == '414' ~ 'Kuwait',
    coe == '417' ~ 'Kyrgyzstan',
    coe == '418' ~ 'Laos',
    coe == '428' ~ 'Latvia',
    coe == '422' ~ 'Lebanon',
    coe == '426' ~ 'Lesotho',
    coe == '430' ~ 'Liberia',
    coe == '434' ~ 'Libya',
    coe == '438' ~ 'Liechtenstein',
    coe == '440' ~ 'Lithuania',
    coe == '442' ~ 'Luxenbourg',
    coe == '446' ~ 'Macao',
    coe == '807' ~ 'North Macedonia',
    coe == '450' ~ 'Madagascar',
    coe == '454' ~ 'Malawi',
    coe == '458' ~ 'Malaysia',
    coe == '462' ~ 'Maldives',
    coe == '466' ~ 'Mali',
    coe == '470' ~ 'Malta',
    coe == '584' ~ 'Marshall Islands',
    coe == '474' ~ 'Martinique',
    coe == '478' ~ 'Mauritania',
    coe == '480' ~ 'Mauritius',
    coe == '175' ~ 'Mayotte',
    coe == '484' ~ 'Mexico',
    coe == '583' ~ 'Micronesia',
    coe == '498' ~ 'Moldova',
    coe == '492' ~ 'Monaco',
    coe == '496' ~ 'Mongolia',
    coe == '499' ~ 'Montenegro',
    coe == '500' ~ 'Montserrat',
    coe == '504' ~ 'Morocco',
    coe == '508' ~ 'Mozambique',
    coe == '104' ~ 'Myanmar',
    coe == '516' ~ 'Namibia',
    coe == '520' ~ 'Nauru',
    coe == '524' ~ 'Nepal',
    coe == '528' ~ 'Netherlands',
    coe == '540' ~ 'New Caledonia',
    coe == '554' ~ 'New Zealand',
    coe == '558' ~ 'Nicaragua',
    coe == '562' ~ 'Niger',
    coe == '566' ~ 'Nigeria',
    coe == '570' ~ 'Niue',
    coe == '574' ~ 'Norfolk Island',
    coe == '580' ~ 'Northern Mariana Islands',
    coe == '578' ~ 'Norway',
    coe == '512' ~ 'Oman',
    coe == '586' ~ 'Pakistan',
    coe == '585' ~ 'Palau',
    coe == '275' ~ 'Palestine',
    coe == '591' ~ 'Panama',
    coe == '598' ~ 'Papua New Guinea',
    coe == '600' ~ 'Paraguay',
    coe == '604' ~ 'Peru',
    coe == '608' ~ 'Philippines',
    coe == '612' ~ 'Pitcairn',
    coe == '616' ~ 'Poland',
    coe == '620' ~ 'Portugal',
    coe == '630' ~ 'Puerto Rico',
    coe == '634' ~ 'Qatar',
    coe == '638' ~ 'Reunion',
    coe == '642' ~ 'Romania',
    coe == '643' ~ 'Russian Federation',
    coe == '646' ~ 'Rwanda',
    coe == '652' ~ 'Saint Barthelemy',
    coe == '654' ~ 'Saint Helena',
    coe == '659' ~ 'Saint Kitts and Nevis',
    coe == '662' ~ 'Saint Lucia',
    coe == '663' ~ 'Saint Martin',
    coe == '666' ~ 'Saint Pierre and Miquelon',
    coe == '670' ~ 'Saint Vincent and the Grenadines',
    coe == '882' ~ 'Samoa',
    coe == '674' ~ 'San Marino',
    coe == '678' ~ 'Sao Tome and Principe',
    coe == '682' ~ 'Saudi Arabia',
    coe == '686' ~ 'Senegal',
    coe == '688' ~ 'Serbia',
    coe == '690' ~ 'Seychelles',
    coe == '694' ~ 'Sierra Leone',
    coe == '702' ~ 'Singapore',
    coe == '534' ~ 'Sint Maarten',
    coe == '703' ~ 'Slovakia',
    coe == '705' ~ 'Slovenia',
    coe == '90' ~ 'Solomon Islands',
    coe == '706' ~ 'Somalia',
    coe == '710' ~ 'South Africa',
    coe == '239' ~ 'South Georgia and the Sandwich Islands',
    coe == '728' ~ 'South Sudan',
    coe == '724' ~ 'Spain',
    coe == '144' ~ 'Sri Lanka',
    coe == '729' ~ 'Sudan',
    coe == '740' ~ 'Suriname',
    coe == '744' ~ 'Svalbard and Jan Mayen',
    coe == '748' ~ 'Eswatini',
    coe == '752' ~ 'Sweden',
    coe == '756' ~ 'Switzerland',
    coe == '760' ~ 'Syria',
    coe == '158' ~ 'Taiwan',
    coe == '762' ~ 'Tajikistan',
    coe == '834' ~ 'Tanzania',
    coe == '764' ~ 'Thailand',
    coe == '768' ~ 'Togo',
    coe == '772' ~ 'Tokelau',
    coe == '776' ~ 'Tonga',
    coe == '780' ~ 'Trinidad and Tobago',
    coe == '788' ~ 'Tunisia',
    coe == '792' ~ 'Turkey',
    coe == '795' ~ 'Turkmenistan',
    coe == '796' ~ 'Turks and Caicos Islands',
    coe == '798' ~ 'Tuvalu',
    coe == '800' ~ 'Uganda',
    coe == '804' ~ 'Ukraine',
    coe == '784' ~ 'United Arab Emirates',
    coe == '826' ~ 'United Kingdom',
    coe == '840' ~ 'United States',
    coe == '581' ~ 'United States Minor Outlying Islands',
    coe == '858' ~ 'Uruguay',
    coe == '860' ~ 'Uzbekistan',
    coe == '876' ~ 'Wallis and Futuna',
    coe == '548' ~ 'Vanuatu',
    coe == '862' ~ 'Venezuela',
    coe == '732' ~ 'Western Sahara',
    coe == '704' ~ 'Vietnam',
    coe == '92' ~ 'Virgin Islands, British',
    coe == '850' ~ 'Virgin Islands, United States',
    coe == '887' ~ 'Yemen',
    coe == '894' ~ 'Zambia',
    coe == '716' ~ 'Zimbabwe',
    coe == '248' ~ 'Aland Islands',
    # Former codes
    coe == '736' ~ 'Former Sudan',
    coe == '886' ~ 'Former South Yemen',
    coe == '200' ~ 'Former Czechoslovakia',
    coe == '810' ~ 'Former Soviet Union',
    coe == '890' ~ 'Former Yugoslavia',
    coe == '891' ~ 'Former Serbia and Montenegro',
    coe == '991' ~ 'Without citizenship',
    TRUE ~ 'Unknown'
  ))

# Let's turn Unknown to NA
sample[country == 'Unknown']$country <- NA

# Eventually, we managed to find the country of ...
sample[!is.na(country),length(shnro)] # 921,412
sample[!is.na(country),length(shnro)] / nrow(sample) # 99.6% 

# Amount of missing per category: Finnish, migrant, 2G migrant, etc.
sample[,table(mig,country=='Finland',useNA='always')]
sample[,round(prop.table(table(mig,country=='Finland',useNA='always'),1),3)] 

# Most common backgrounds
head(sample[,as.data.table(table(country))][order(N,decreasing = TRUE)],12)

# Note that we have the category Soviet, which I will divide based on language (not perfect solution)
sample[country == 'Former Soviet Union']
sample[country == 'Former Soviet Union',table(language)] # There are 434 Estonian speakers here, also Finnish speakers...
sample[country == 'Former Soviet Union',table(mig,language)] # but the Finnish speakers are mostly migrant descendants

# Use current republics
sample[country == 'Former Soviet Union' & language == 'be']$country <- 'Belarus'
sample[country == 'Former Soviet Union' & language == 'uk']$country <- 'Ukraine'
sample[country == 'Former Soviet Union' & language %in% c('mo','ro')]$country <- 'Moldova'
sample[country == 'Former Soviet Union' & language == 'lt']$country <- 'Lithuania'
sample[country == 'Former Soviet Union' & language == 'lv']$country <- 'Latvia'
sample[country == 'Former Soviet Union' & language == 'et']$country <- 'Estonia'
sample[country == 'Former Soviet Union' & language == 'ka']$country <- 'Georgia'
sample[country == 'Former Soviet Union' & language == 'hy']$country <- 'Armenia'
sample[country == 'Former Soviet Union' & language == 'az']$country <- 'Azerbaijan'
sample[country == 'Former Soviet Union' & language == 'kk']$country <- 'Kazakhstan'
sample[country == 'Former Soviet Union' & language == 'ky']$country <- 'Kyrgyzstan'
sample[country == 'Former Soviet Union' & language == 'tg']$country <- 'Tajikistan'
sample[country == 'Former Soviet Union' & language == 'tk']$country <- 'Turkmenistan'
sample[country == 'Former Soviet Union' & language == 'uz']$country <- 'Uzbekistan'
sample[country == 'Former Soviet Union']$country <- 'Russian Federation' # Otherwise, just Russia

# Most common backgrounds among migrant descendants
head(sample[mig %in% c('1.5G migrant','2G migrant','2G migrant (mixed parents)','2.5 migrant'),
            as.data.table(table(country))][order(N,decreasing = TRUE)],12)

# Final sample
sample[,table(mig,country=='Finland',useNA='always')]

########################################################################################################################

# Let's reduce the countries into broader categories

sample <- sample %>%
  mutate(background = case_when(
    # Finland
    country == 'Finland' ~ 'Finland',
    country == 'Aland Islands' ~ 'Finland',
    # Estonia
    country == 'Estonia' ~ 'Estonia',
    # Russia
    country == 'Russian Federation' ~ 'Russia',
    # Scandinavia
    country == 'Denmark' ~ 'Scandinavia',
    country == 'Faroe Islands' ~ 'Scandinavia',
    country == 'Iceland' ~ 'Scandinavia',
    country == 'Greenland' ~ 'Scandinavia',
    country == 'Norway' ~ 'Scandinavia',
    country == 'Sweden' ~ 'Scandinavia',
    # Western Europe
    country == 'Austria' ~ 'Western Europe',
    country == 'Germany' ~ 'Western Europe',
    country == 'Liechtenstein' ~ 'Western Europe',
    country == 'Switzerland' ~ 'Western Europe',
    country == 'Belgium' ~ 'Western Europe',
    country == 'Netherlands' ~ 'Western Europe',
    country == 'Luxenbourg' ~ 'Western Europe',
    country == 'Jersey' ~ 'Western Europe',
    country == 'Guernsey' ~ 'Western Europe',
    country == 'Ireland' ~ 'Western Europe',
    country == 'Isle of Man' ~ 'Western Europe',
    country == 'United Kingdom' ~ 'Western Europe',
    country == 'France' ~ 'Western Europe',
    # Southern Europe
    country == 'Andorra' ~ 'Southern Europe',
    country == 'Portugal' ~ 'Southern Europe',
    country == 'Spain' ~ 'Southern Europe',
    country == 'Malta' ~ 'Southern Europe',
    country == 'Monaco' ~ 'Southern Europe',
    country == 'Holy See (Vatican City State)' ~ 'Southern Europe',
    country == 'Italy' ~ 'Southern Europe',
    country == 'San Marino' ~ 'Southern Europe',
    country == 'Gibraltar' ~ 'Southern Europe',
    country == 'Cyprus' ~ 'Southern Europe',
    country == 'Greece' ~ 'Southern Europe',
    # Central and Eastern Europe
    country == 'Czechia' ~ 'Central & Eastern Europe',
    country == 'Hungary' ~ 'Central & Eastern Europe',
    country == 'Poland' ~ 'Central & Eastern Europe',
    country == 'Slovakia' ~ 'Central & Eastern Europe',
    country == 'Former Czechoslovakia' ~ 'Central & Eastern Europe',
    country == 'Latvia' ~ 'Central & Eastern Europe',
    country == 'Lithuania' ~ 'Central & Eastern Europe',
    country == 'Belarus' ~ 'Central & Eastern Europe',
    country == 'Moldova' ~ 'Central & Eastern Europe',
    country == 'Ukraine' ~ 'Central & Eastern Europe',
    country == 'Albania' ~ 'Central & Eastern Europe',
    country == 'Bosnia and Herzegovina' ~ 'Central & Eastern Europe',
    country == 'Croatia' ~ 'Central & Eastern Europe',
    country == 'North Macedonia' ~ 'Central & Eastern Europe',
    country == 'Montenegro' ~ 'Central & Eastern Europe',
    country == 'Serbia' ~ 'Central & Eastern Europe',
    country == 'Slovenia' ~ 'Central & Eastern Europe',
    country == 'Former Yugoslavia' ~ 'Central & Eastern Europe',
    country == 'Former Serbia and Montenegro' ~ 'Central & Eastern Europe',
    country == 'Bulgaria' ~ 'Central & Eastern Europe',
    country == 'Romania' ~ 'Central & Eastern Europe',
    # Caucasus
    country == 'Azerbaijan' ~ 'Caucasus',
    country == 'Armenia' ~ 'Caucasus',
    country == 'Georgia' ~ 'Caucasus',
    # Central Asia
    country == 'Kazakhstan' ~ 'Central Asia',
    country == 'Kyrgyzstan' ~ 'Central Asia',
    country == 'Tajikistan' ~ 'Central Asia',
    country == 'Turkmenistan' ~ 'Central Asia',
    country == 'Uzbekistan' ~ 'Central Asia',
    # Israel
    country == 'Israel' ~ 'Israel',
    # Middle East
    country == 'Afghanistan' ~ 'Middle East',
    country == 'Bahrain' ~ 'Middle East',
    country == 'Iran' ~ 'Middle East',
    country == 'Iraq' ~ 'Middle East',
    country == 'Jordan' ~ 'Middle East',
    country == 'Kuwait' ~ 'Middle East',
    country == 'Lebanon' ~ 'Middle East',
    country == 'Oman' ~ 'Middle East',
    country == 'Palestine' ~ 'Middle East',
    country == 'Qatar' ~ 'Middle East',
    country == 'Saudi Arabia' ~ 'Middle East',
    country == 'United Arab Emirates' ~ 'Middle East',
    country == 'Syria' ~ 'Middle East',
    country == 'Turkey' ~ 'Middle East',
    country == 'Yemen' ~ 'Middle East',
    country == 'Former South Yemen' ~ 'Middle East',
    # South Asia
    country == 'Bangladesh' ~ 'South Asia',
    country == 'Bhutan' ~ 'South Asia',
    country == 'India' ~ 'South Asia',
    country == 'Nepal' ~ 'South Asia',
    country == 'Myanmar' ~ 'South Asia',
    country == 'Pakistan' ~ 'South Asia',
    country == 'Sri Lanka' ~ 'South Asia',
    # East Asia
    country == 'China' ~ 'East Asia',
    country == 'Hong Kong' ~ 'East Asia',
    country == 'Japan' ~ 'East Asia',
    country == 'North Korea' ~ 'East Asia',
    country == 'South Korea' ~ 'East Asia',
    country == 'Macao' ~ 'East Asia',
    country == 'Mongolia' ~ 'East Asia',
    country == 'Taiwan' ~ 'East Asia',
    # South East Asia
    country == 'Brunei Darussalam' ~ 'Southeast Asia',
    country == 'Cambodia' ~ 'Southeast Asia',
    country == 'East Timor' ~ 'Southeast Asia',
    country == 'Indonesia' ~ 'Southeast Asia',
    country == 'Malaysia' ~ 'Southeast Asia',
    country == 'Laos' ~ 'Southeast Asia',
    country == 'Philippines' ~ 'Southeast Asia',
    country == 'Singapore' ~ 'Southeast Asia',
    country == 'Thailand' ~ 'Southeast Asia',
    country == 'Vietnam' ~ 'Southeast Asia',    
    # Pacific
    country == 'American Samoa' ~ 'Pacific',
    country == 'Christmas Island' ~ 'Pacific',
    country == 'Cook Islands' ~ 'Pacific',
    country == 'Fiji' ~ 'Pacific',
    country == 'French Polynesia' ~ 'Pacific',
    country == 'Guam' ~ 'Pacific',
    country == 'Kiribati' ~ 'Pacific',
    country == 'Marshall Islands' ~ 'Pacific',
    country == 'Micronesia' ~ 'Pacific',
    country == 'Nauru' ~ 'Pacific',
    country == 'New Caledonia' ~ 'Pacific',
    country == 'Niue' ~ 'Pacific',
    country == 'Norfolk Island' ~ 'Pacific',
    country == 'Northern Mariana Islands' ~ 'Pacific',
    country == 'Palau' ~ 'Pacific',
    country == 'Papua New Guinea' ~ 'Pacific',
    country == 'Pitcairn' ~ 'Pacific',
    country == 'Samoa' ~ 'Pacific',
    country == 'Solomon Islands' ~ 'Pacific',
    country == 'Tonga' ~ 'Pacific',
    country == 'Tuvalu' ~ 'Pacific',
    country == 'Tokelau' ~ 'Pacific',
    country == 'Vanuatu' ~ 'Pacific',
    country == 'Wallis and Futuna' ~ 'Pacific',
    # North America, Australia, New Zealand
    country == 'Canada' ~ 'US, Canada, Australia, New Zealand',
    country == 'United States' ~ 'US, Canada, Australia, New Zealand',
    country == 'United States Minor Outlying Islands' ~ 'US, Canada, Australia, New Zealand',
    country == 'Australia' ~ 'US, Canada, Australia, New Zealand',
    country == 'New Zealand' ~ 'US, Canada, Australia, New Zealand',
    # Latin America
    country == 'Anguilla' ~ 'Latin America',
    country == 'Antigua and Barbuda' ~ 'Latin America',
    country == 'Bahamas' ~ 'Latin America',
    country == 'Barbados' ~ 'Latin America',
    country == 'Argentina' ~ 'Latin America',
    country == 'Aruba' ~ 'Latin America',
    country == 'Belize' ~ 'Latin America',
    country == 'Bermuda' ~ 'Latin America',
    country == 'Bolivia' ~ 'Latin America',
    country == 'Bonaire, Sint Eustatius and Saba' ~ 'Latin America',
    country == 'Brazil' ~ 'Latin America',
    country == 'Cayman Islands' ~ 'Latin America',
    country == 'Chile' ~ 'Latin America',
    country == 'Cocos (Keeling) Islands' ~ 'Latin America',
    country == 'Colombia' ~ 'Latin America',
    country == 'Costa Rica' ~ 'Latin America',
    country == 'Cuba' ~ 'Latin America',
    country == 'Curacao' ~ 'Latin America',
    country == 'Dominica' ~ 'Latin America',
    country == 'Dominican Republic' ~ 'Latin America',
    country == 'Ecuador' ~ 'Latin America',
    country == 'El Salvador' ~ 'Latin America',
    country == 'Falkland Islands (Malvinas)' ~ 'Latin America',
    country == 'French Guiana' ~ 'Latin America',
    country == 'Grenada' ~ 'Latin America',
    country == 'Guadeloupe' ~ 'Latin America',
    country == 'Guatemala' ~ 'Latin America',
    country == 'Guyana' ~ 'Latin America',
    country == 'Haiti' ~ 'Latin America',
    country == 'Honduras' ~ 'Latin America',
    country == 'Jamaica' ~ 'Latin America',
    country == 'Maldives' ~ 'Latin America',
    country == 'Martinique' ~ 'Latin America',
    country == 'Mexico' ~ 'Latin America',
    country == 'Montserrat' ~ 'Latin America',
    country == 'Nicaragua' ~ 'Latin America',
    country == 'Panama' ~ 'Latin America',
    country == 'Paraguay' ~ 'Latin America',
    country == 'Peru' ~ 'Latin America',
    country == 'Puerto Rico' ~ 'Latin America',
    country == 'Saint Barthelemy' ~ 'Latin America',
    country == 'Saint Kitts and Nevis' ~ 'Latin America',
    country == 'Saint Lucia' ~ 'Latin America',
    country == 'Saint Martin' ~ 'Latin America',
    country == 'Saint Vincent and the Grenadines' ~ 'Latin America',
    country == 'Sint Maarten' ~ 'Latin America',
    country == 'Suriname' ~ 'Latin America',
    country == 'Trinidad and Tobago' ~ 'Latin America',
    country == 'Turks and Caicos Islands' ~ 'Latin America',
    country == 'Uruguay' ~ 'Latin America',
    country == 'Venezuela' ~ 'Latin America',
    country == 'Virgin Islands, British' ~ 'Latin America',
    country == 'Virgin Islands, United States' ~ 'Latin America',
    # North Africa
    country == 'Algeria' ~ 'North Africa',
    country == 'Egypt' ~ 'North Africa',
    country == 'Libya' ~ 'North Africa',
    country == 'Mauritania' ~ 'North Africa',
    country == 'Morocco' ~ 'North Africa',
    country == 'South Sudan' ~ 'North Africa',
    country == 'Sudan' ~ 'North Africa',
    country == 'Tunisia' ~ 'North Africa',
    country == 'Western Sahara' ~ 'North Africa',
    country == 'Former Sudan' ~ 'North Africa',
    # West Africa
    country == 'Benin' ~ 'West Africa',
    country == 'Burkina Faso' ~ 'West Africa',
    country == 'Cabo Verde' ~ 'West Africa',
    country == 'Chad' ~ 'West Africa',
    country == 'Djibouti' ~ 'West Africa',
    country == 'Ivory Coast' ~ 'West Africa',
    country == 'Gambia' ~ 'West Africa',
    country == 'Ghana' ~ 'West Africa',
    country == 'Guinea' ~ 'West Africa',
    country == 'Guinea-Bissau' ~ 'West Africa',
    country == 'Liberia' ~ 'West Africa',
    country == 'Mali' ~ 'West Africa',
    country == 'Niger' ~ 'West Africa',
    country == 'Nigeria' ~ 'West Africa',
    country == 'Senegal' ~ 'West Africa',
    country == 'Sierra Leone' ~ 'West Africa',
    country == 'Togo' ~ 'West Africa',
    # East Africa
    country == 'Burundi' ~ 'East Africa',
    country == 'Comoros' ~ 'East Africa',
    country == 'Eritrea' ~ 'East Africa',
    country == 'Ethiopia' ~ 'East Africa',
    country == 'Kenya' ~ 'East Africa',
    country == 'Madagascar' ~ 'East Africa',
    country == 'Mauritius' ~ 'East Africa',
    country == 'Mayotte' ~ 'East Africa',
    country == 'Mozambique' ~ 'East Africa',
    country == 'Reunion' ~ 'East Africa',
    country == 'Rwanda' ~ 'East Africa',
    country == 'Seychelles' ~ 'East Africa',
    country == 'Somalia' ~ 'East Africa',
    country == 'Tanzania' ~ 'East Africa',
    country == 'Uganda' ~ 'East Africa',
    # Central & Southern Africa
    country == 'Cameroon' ~ 'Central and Southern Africa',
    country == 'Central African Republic' ~ 'Central and Southern Africa',
    country == 'Congo (Congo-Brazzaville)' ~ 'Central and Southern Africa',
    country == 'Congo, The Democratic Republic of' ~ 'Central and Southern Africa',
    country == 'Equatorial Guinea' ~ 'Central and Southern Africa',
    country == 'Gabon' ~ 'Central and Southern Africa',
    country == 'Sao Tome and Principe' ~ 'Central and Southern Africa',
    country == 'Angola' ~ 'Central and Southern Africa',
    country == 'Botswana' ~ 'Central and Southern Africa',
    country == 'Lesotho' ~ 'Central and Southern Africa',
    country == 'Malawi' ~ 'Central and Southern Africa',
    country == 'Namibia' ~ 'Central and Southern Africa',
    country == 'South Africa' ~ 'Central and Southern Africa',
    country == 'Eswatini' ~ 'Central and Southern Africa',
    country == 'Zambia' ~ 'Central and Southern Africa',
    country == 'Zimbabwe' ~ 'Central and Southern Africa',
    # Other
    country == 'Antarctica' ~ 'Other',
    country == 'Bouvet Island' ~ 'Other',
    country == 'British Indian Ocean Territory' ~ 'Other',
    country == 'French Southern Territories' ~ 'Other',
    country == 'Heard Island and McDonald Island' ~ 'Other',
    country == 'Saint Helena' ~ 'Other',
    country == 'Saint Pierre and Miquelon' ~ 'Other',
    country == 'South Georgia and the Sandwich Islands' ~ 'Other',
    country == 'Svalbard and Jan Mayen' ~ 'Other',
    # No citizenship
    country == 'Without citizenship' ~ '',
    TRUE ~ 'Unknown'
  ))

# Redefine unknown as NA
sample[background == 'Unknown']$background <- NA

# Most common backgrounds
head(sample[,as.data.table(table(background))][order(N,decreasing = TRUE)],12)
# Seems we have now more representation among certain groups, but do we have enough for differentiating across generations...
head(sample[mig %in% c('1.5G migrant','2G migrant','2G migrant (mixed parents)','2.5 migrant'),
            as.data.table(table(background))][order(N,decreasing = TRUE)],12)

########################################################################################################################

# Save image
save(sample,file='data2.sample.with.backgrounds.RData')

########################################################################################################################