# This is the R code to accompany the report 'Examining the intersection of NCDs and COVID-19: Lessons and opportunities from emerging data' which is available at: https://eiuperspectives.economist.com/healthcare/examining-intersection-ncds-and-covid-19-lessons-and-opportunities-emerging-data
# The analyses uses data from the Global Burden of Disease, which is not freely available to all - the accompanying datasets (which had been partly preprocessed) have not therefore been shared
# The code was prepared by Peter Tennant (Causal Insights Ltd) for the Economist Intelligence Unit and run using Microsoft R Open 4.0.2  

#### PREPARATION ####

#Load required packages
require(MASS)
library(tidyverse)
library(dplyr)
require(readr) 
require(ggplot2)
require(forecast)
require(zoo)
require(bestglm)
require(rJava)
require(glmulti)

#Define list of LMICs for examination
List_countryid  <- list("AFG", "ALB", "DZA", "AGO", "ARG", "ARM", "AZE", "BGD", "BLR", "BLZ", "BEN", "BTN", "BOL", "BIH", "BWA", "BRA", "BGR", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "CHN", "COL", "COM", "COD", "COG", "CRI", "CIV", "CUB", "DJI", "DOM", "ECU", "EGY", "SLV", "GNQ", "ERI", "SWZ", "ETH", "FJI", "GAB", "GMB", "GEO", "GHA", "GTM", "GIN", "GNB", "GUY", "HTI", "HND", "IND", "IDN", "IRN", "IRQ", "JAM", "JOR", "KAZ", "KEN", "KGZ", "LBN", "LSO", "LBR", "LBY", "MDG", "MWI", "MYS", "MDV", "MLI", "MRT", "MEX", "MNG", "MNE", "MAR", "MOZ", "MMR", "NAM", "NPL", "NIC", "NER", "NGA", "MKD", "PAK", "PNG", "PRY", "PER", "PHL", "RUS", "RWA", "STP", "SEN", "SRB", "SLE", "SOM", "ZAF", "LKA", "SDN", "SUR", "SYR", "TJK", "TZA", "THA", "TGO", "TUN", "TUR", "UGA", "UKR", "UZB", "VEN", "VNM", "YEM", "ZMB", "ZWE")

#Import 'outcome' (COVID) data and reduce to relevant list of countries
COVID_data      <- read_csv("...")
COVID_data      <- subset(COVID_data, countryid %in% List_countryid)

#Reduce to relevant variables:
COVID_data <- COVID_data %>% 
  select(-c(total_cases_per_million, new_cases_per_million, total_deaths_per_million, new_deaths_per_million, population_density, median_age, aged_65_older, aged_70_older, life_expectancy)) %>% 
  group_by(countryid) %>%
  mutate(day_in_data      = row_number()) %>% 
  ungroup(countryid)

# Interpolate missing values:
# Create a blank dataframe to save smoothed values
COVID_smoothed <- tibble(data.frame(countryid = character(0), day_in_data=numeric(0),  s_stringency_index=numeric(0), s_reproduction_rate=numeric(0), stringsAsFactors=F))

#Create list of countries requiring interpolation:
List_imputing  <- list("AFG", "ALB", "DZA", "AGO", "ARG", "AZE", "BGD", "BLR", "BLZ", "BEN", "BTN", "BOL", "BIH", "BWA", "BRA", "BGR", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "CHN", "COL", "COD", "COG", "CRI", "CIV", "CUB", "DJI", "DOM", "ECU", "EGY", "SLV", "ERI", "SWZ", "ETH", "GAB", "GMB", "GEO", "GHA", "GTM", "GIN", "GUY", "HTI", "HND", "IND", "IDN", "IRN", "IRQ", "JAM", "JOR", "KAZ", "KEN", "KGZ", "LBN", "LSO", "LBR", "LBY", "MDG", "MWI", "MYS", "MLI", "MRT", "MEX", "MNG", "MAR", "MOZ", "MMR", "NAM", "NPL", "NIC", "NER", "NGA", "PAK", "PNG", "PRY", "PER", "PHL", "RUS", "RWA", "SEN", "SRB", "SLE", "SOM", "ZAF", "LKA", "SDN", "SUR", "SYR", "TJK", "TZA", "THA", "TGO", "TUN", "TUR", "UGA", "UKR", "UZB", "VEN", "VNM", "YEM", "ZMB", "ZWE")

lapply(List_imputing, function(COUNTRY) {

Data_Time_Series <- ts(data = subset(COVID_data, countryid == COUNTRY))

#Create smoothed values:
imputed_stringency    <- na.interp(Data_Time_Series[,"stringency_index"], lambda = "auto")
imputed_rr            <- na.interp(Data_Time_Series[,"reproduction_rate"], lambda = "auto")
imputed_data          <- do.call(rbind, Map(data.frame, countryid=COUNTRY, day_in_data=Data_Time_Series[,"day_in_data"], s_stringency_index=imputed_stringency, s_reproduction_rate=imputed_rr))

COVID_smoothed        <<- tibble(rbind(COVID_smoothed, imputed_data))

})

# Merge interpolated values into main dataset
COVID_smoothed        <- tibble(remove_rownames(COVID_smoothed))
COVID_data            <- left_join(COVID_data, COVID_smoothed, by=c("countryid", "day_in_data"))

#Remove imputation material:
rm(COVID_smoothed, List_imputing)

#Replace remaining missing values with zeros:
COVID_data$new_cases                <- ifelse(is.na(COVID_data$new_cases), 0, COVID_data$new_cases)
COVID_data$new_cases_per_million    <- (COVID_data$new_cases/COVID_data$population)*1000000


COVID_data        <- COVID_data %>% 
  group_by(countryid) %>% 
  mutate(date = as.Date(date, format="%d/%m/%Y", origin="1970-01-01")) %>% 
  mutate(total_cases_prior     = lag(total_cases, n=18),
         mean_stringency       = round(mean(s_stringency_index, na.rm = TRUE), digits=1),
         mean_reproduction     = round(mean(s_reproduction_rate, na.rm = TRUE), digits=1),
         change_stringency     = s_stringency_index - lag(s_stringency_index),
         seven_day_average     = rollmean(new_cases_per_million, 7, fill=NA),
         first_lockdown_flag   = if_else(change_stringency>0 & s_stringency_index>=50 & date<as.Date("2020-09-01"), 1, 0),
         sum_first_flag        = cumsum(first_lockdown_flag),
         first_lockdown_enter  = if_else(change_stringency>0 & s_stringency_index>=50 & date<as.Date("2020-09-01") & sum_first_flag<=1, 1, 0),
         rate_at_first         = round(sum(if_else(first_lockdown_enter==1, seven_day_average, 0)), digits=3),
         second_lockdown_flag  = if_else(change_stringency>0 & s_stringency_index>=50 & date>=as.Date("2020-09-01") & date<as.Date("2021-04-01"), 1, 0),
         sum_second_flag       = cumsum(second_lockdown_flag),
         second_lockdown_enter = if_else(change_stringency>0 & s_stringency_index>=50 & date>=as.Date("2020-09-01") & date<as.Date("2021-04-01") & sum_second_flag<=1, 1, 0),
         rate_at_second        = round(sum(if_else(second_lockdown_enter==1, seven_day_average, 0))), digits=3) %>%
  filter(date == as.Date("2021-04-01")) %>% 
  select(-c(new_cases, new_deaths, reproduction_rate, stringency_index, day_in_data, s_stringency_index, s_reproduction_rate, new_cases_per_million, change_stringency, seven_day_average, first_lockdown_flag, sum_first_flag, first_lockdown_enter, second_lockdown_flag, sum_second_flag, second_lockdown_enter))

#Import 'exposure' data from the Global Burden of Disease 2019
NCD_data <- read_csv("...")

#Convert into multiple columns per country 
NCD_data <- NCD_data %>% 
  unite(ncds, ncd_group, measure, sep="_") %>% 
  spread(ncds, n)

#Reduce to relevant countries only
NCD_data        <- subset(NCD_data, countryid %in% List_countryid)

#Import other country-level data and reduce to relevant
Country_Info_data <- read_csv("...")
Country_Info_data <- subset(Country_Info_data, countryid %in% List_countryid)

#Import age and sex data
Age_Sex_data <- read_csv("...")
Age_Sex_data <- subset(Age_Sex_data, countryid %in% List_countryid)

#Drop with missing age/sex data:
Age_Sex_data <- Age_Sex_data %>% 
  filter(is.na(m80)==FALSE)

#Import total testing data
Test_data <- read_csv("...")

#Drop rows with missing data
Test_data        <- Test_data[complete.cases(Test_data$cumulativetests),] 
  
#Retain highest value
Test_data        <- Test_data %>% 
  group_by(countryid) %>% 
  mutate(maxtests = max(cumulativetests)) %>% 
  summarise(maxtests=first(maxtests))

Test_data$maxtests  <- ifelse(Test_data$maxtests == -Inf, NA, Test_data$maxtests)
Test_data <- subset(Test_data, countryid %in% List_countryid)

#Join datasets
COVID_NCD_data  <- inner_join(COVID_data, NCD_data, by=c("countryid"))
COVID_NCD_data  <- inner_join(COVID_NCD_data, Country_Info_data, by=c("countryid"))
COVID_NCD_data  <- inner_join(COVID_NCD_data, Age_Sex_data, by=c("countryid"))

#Make new sex and age variables:
COVID_NCD_data    <- COVID_NCD_data %>% 
  mutate(m0_49    = (m0_14 + m15_19 + m20_24 + m25_29 + m30_34 + m35_39 + m40_44 + m45_49),
         m50_69   = (m50_54  + m55_59 + m60_64 + m65_69),
         m70      = (m70_74 + m75_79 + m80),
         f0_49    = (f0_14 + f15_19 + f20_24 + f25_29 + f30_34 + f35_39 + f40_44 + f45_49),
         f50_69   = (f50_54  + f55_59 + f60_64 + f65_69),
         f70      = (f70_74 + f75_79 + f80)) 

#### DESCRIPTIVE ANALYSIS ####

#Save some important summary numbers:  

mean_population                <- mean(COVID_NCD_data$pop_total)
mean_total_cases               <- mean(COVID_NCD_data$total_cases)
mean_total_cases_prior         <- mean(COVID_NCD_data$total_cases_prior)
mean_total_deaths              <- mean(COVID_NCD_data$total_deaths)
var_total_deaths               <- var(COVID_NCD_data$total_deaths)
mean_ncd_deaths                <- mean(COVID_NCD_data$all_ncds_deaths)
mean_ncd_deaths_pop            <- mean_ncd_deaths/mean_population
mean_ncd_dalys                 <- mean(COVID_NCD_data$all_ncds_dalys)
mean_ncd_dalys_pop             <- mean_ncd_dalys/mean_population
mean_cfr                       <- mean(COVID_NCD_data$total_deaths/COVID_NCD_data$total_cases)
variance_cfr                   <- var(COVID_NCD_data$total_deaths/COVID_NCD_data$total_cases)

#Describe the sample

#Create variables for summarising
COVID_NCD_data <- COVID_NCD_data %>% 
  mutate(ncd_deaths_1k       = all_ncds_deaths/1000,
         ncd_mortality       = all_ncds_deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         ncd_dalys_1k        = all_ncds_dalys/1000,
         ncd_daly_ratio      = all_ncds_dalys/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_cases_1k      = total_cases/1000,
         covid_incidence     = total_cases/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*1000,
         covid_deaths_1k     = total_deaths/1000,
         covid_mortality     = total_deaths/(m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)*100000,
         covid_cfr           = (total_deaths/total_cases)*100,
         total_pop_1m        = (m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70)/1000,
         total_males         = (m0_49 + m50_69 + m70)/1000/total_pop_1m*100,
         total_0_49          = (m0_49 + f0_49)/1000/total_pop_1m*100,
         total_50_69         = (m50_69 + f50_69)/1000/total_pop_1m*100,
         total_70            = (m70 + f70)/1000/total_pop_1m*100,
         gdp_1k              = gdp_per_capita/1000,
         doctors_1k          = doctors_10k/10,
         country_area_100k   = country_area/1000,
         id                  = row_number())

#List of variables to be summarised
List_Cont_Vars  <- list("ncd_deaths_1k", "ncd_mortality", "ncd_dalys_1k", "ncd_daly_ratio", "total_deaths", "covid_mortality", "covid_cases_1k", "covid_incidence", "covid_cfr", "total_pop_1m", "total_males", "total_0_49", "total_50_69", "total_70", "gdp_1k", "percent_urban", "smoking2018", "obesity2016", "doctors_1k", "hospital_beds_per_thousand", "healthcare_exp", "ph_exp", "mean_stringency", "mean_reproduction", "country_area_100k", "pollution")

#Create black dataframe to store summary information
Cont_Sum    <- data.frame(variable=character(0), min=numeric(0), max=numeric(0), median=numeric(0), q1=numeric(0), q3=numeric(0), n=numeric(0), stringsAsFactors=F)

for (i in 1:length(List_Cont_Vars)) {
  
  matrix        <-  cbind(NULL, List_Cont_Vars[i], as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0, na.rm=TRUE)), as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=1, na.rm=TRUE)), as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.5, na.rm=TRUE)), as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.25, na.rm=TRUE)), as.numeric(quantile(COVID_NCD_data[,unlist(List_Cont_Vars[i])], probs=0.75, na.rm=TRUE)), as.numeric(length(COVID_NCD_data$id[which(is.na(COVID_NCD_data[,unlist(List_Cont_Vars[i])])==FALSE)])))
  Cont_Sum      <-  rbind(Cont_Sum, as.numeric(matrix[1,]))
}

#### PRIMARY ANALYSIS ####

#Now to determine the best adjustment set:
COVID_NCD_data$all_ncds_deaths <- round(COVID_NCD_data$all_ncds_deaths, digits=0)
COVID_NCD_data$all_ncds_dalys  <- round(COVID_NCD_data$all_ncds_dalys, digits=0)  

#First adjust for age, sex, and population size:
best_model <- glmulti(total_deaths ~ m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + mean_stringency + mean_reproduction + gdp_per_capita + hospital_beds_per_thousand + country_area + smoking2018 + obesity2016 + doctors_10k + healthcare_exp + ph_exp + percent_urban + pollution, data = COVID_NCD_data, fitfunction = "glm", offset=log(COVID_NCD_data$total_cases_prior), family=gaussian(link="log"), crit=bic, method="l", level=1)
summary(best_model)

# The best model - according to the BIC - includes the age/sex variables plus smoking2018 and healthcare_exp

#Run estimation models

#MODEL 1
#Age and sex adjustment only
#Deaths model 1
model1_deaths                  <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_deaths + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70, data=COVID_NCD_data, family=quasipoisson())
summary(model1_deaths)
results1_deaths                <- summary(model1_deaths)$coefficients

#Call coefficient and CIs
summary(model1_deaths)$coefficients[2,]
confint(model1_deaths, "all_ncds_deaths")
#Transformed for interpretation:
exp(summary(model1_deaths)$coefficients[2]*0.1*mean_ncd_deaths)
exp(confint(model1_deaths, "all_ncds_deaths")*0.1*mean_ncd_deaths)

#DALYs model 1
model1_dalys                  <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_dalys + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70, data=COVID_NCD_data, family=quasipoisson())
summary(model1_dalys)
results1_dalys                <-  summary(model1_dalys)$coefficients

#Call coefficient and CIs
summary(model1_dalys)$coefficients[2,]
confint(model1_dalys, "all_ncds_dalys")
#Transformed for interpretation:
exp(summary(model1_dalys)$coefficients[2]*0.1*mean_ncd_dalys)
exp(confint(model1_dalys, "all_ncds_dalys")*0.1*mean_ncd_dalys)

#MODEL 2
#Age, sex, and bic variables
#Deaths model 2
model2_deaths              <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_deaths + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_data, family=quasipoisson(), na.action=na.exclude)
summary(model2_deaths)
results2_deaths            <- summary(model2_deaths)$coefficients

#Call coefficient and CIs
summary(model2_deaths)$coefficients[2,]
confint(model2_deaths, "all_ncds_deaths")
#Transformed for interpretation:
exp(summary(model2_deaths)$coefficients[2]*0.1*mean_ncd_deaths)
exp(confint(model2_deaths, "all_ncds_deaths")*0.1*mean_ncd_deaths)

#DALYs model 2
model2_dalys              <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_dalys + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_data, family=quasipoisson())
summary(model2_dalys)
results2_dalys            <- summary(model2_dalys)$coefficients

#Call coefficient and CIs
summary(model2_dalys)$coefficients[2,]
confint(model2_dalys, "all_ncds_dalys")
#Transformed for interpretation:
exp(summary(model2_dalys)$coefficients[2]*0.1*mean_ncd_dalys)
exp(confint(model2_dalys, "all_ncds_dalys")*0.1*mean_ncd_dalys)

#MODEL 3
#Adjusting for all confounders and competing exposures (apart from tests performed)
#Deaths model 3
model3_deaths              <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_deaths + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + mean_stringency + mean_reproduction + gdp_per_capita + hospital_beds_per_thousand + country_area + smoking2018 + obesity2016 + doctors_10k + healthcare_exp + ph_exp + percent_urban + pollution, data=COVID_NCD_data, family=quasipoisson())
summary(model3_deaths)
results3_deaths            <- summary(model3_deaths)$coefficients

#Call coefficient and CIs
summary(model3_deaths)$coefficients[2,]
confint(model3_deaths, "all_ncds_deaths")
#Transformed for interpretation:
exp(summary(model3_deaths)$coefficients[2]*0.1*mean_ncd_deaths)
exp(confint(model3_deaths, "all_ncds_deaths")*0.1*mean_ncd_deaths)

#DALYs model 3
model3_dalys              <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_dalys + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + mean_stringency + mean_reproduction + gdp_per_capita + hospital_beds_per_thousand + country_area + smoking2018 + obesity2016 + doctors_10k + healthcare_exp + ph_exp + percent_urban  + pollution, data=COVID_NCD_data, family=quasipoisson())
summary(model3_dalys)
results3_dalys            <- summary(model3_dalys)$coefficients

#Call coefficient and CIs
summary(model3_dalys)$coefficients[2,]
confint(model3_dalys, "all_ncds_dalys")
#Transformed for interpretation:
exp(summary(model3_dalys)$coefficients[2]*0.1*mean_ncd_dalys)
exp(confint(model3_dalys, "all_ncds_dalys")*0.1*mean_ncd_dalys)

#### SENSITIVITY ANALYSIS OF TESTS DATA ####

COVID_NCD_TEST_data  <- inner_join(COVID_NCD_data, Test_data, by=c("countryid"))

#MODEL 2A
# Identical to Model 2, but in reduced sample with available data
#Deaths model 2a
model2a_deaths             <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_deaths + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_TEST_data, family=quasipoisson())
summary(model2a_deaths)
results2a_deaths           <- summary(model2a_deaths)$coefficients

#Call coefficient and CIs
summary(model2a_deaths)$coefficients[2,]
confint(model2a_deaths, "all_ncds_deaths")
#Transformed for interpretation:
exp(summary(model2a_deaths)$coefficients[2]*0.1*mean_ncd_deaths)
exp(confint(model2a_deaths, "all_ncds_deaths")*0.1*mean_ncd_deaths)

#DALYs model 2a
model2a_dalys           <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_dalys + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_TEST_data, family=quasipoisson())
summary(model2a_dalys)
results_dalys_no_test   <- summary(model2a_dalys)$coefficients

#Call coefficient and CIs
summary(model2a_dalys)$coefficients[2,]
confint(model2a_dalys, "all_ncds_dalys")
#Transformed for interpretation:
exp(summary(model2a_dalys)$coefficients[2]*0.1*mean_ncd_dalys)
exp(confint(model2a_dalys, "all_ncds_dalys")*0.1*mean_ncd_dalys)

#MODEL 2B
# Model 2A with adjustment for tests performed 
#Deaths model 2b
model2b_deaths             <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_deaths + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp + maxtests, data=COVID_NCD_TEST_data, family=quasipoisson())
summary(model2b_deaths)
results2b_deaths           <- summary(model2b_deaths)$coefficients

#Call coefficient and CIs
summary(model2b_deaths)$coefficients[2,]
confint(model2b_deaths, "all_ncds_deaths")
#Transformed for interpretation:
exp(summary(model2b_deaths)$coefficients[2]*0.1*mean_ncd_deaths)
exp(confint(model2b_deaths, "all_ncds_deaths")*0.1*mean_ncd_deaths)

#DALYs model 2b
model2b_dalys           <- glm(total_deaths ~ offset(log(total_cases)) + all_ncds_dalys + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp + maxtests, data=COVID_NCD_TEST_data, family=quasipoisson())
summary(model2b_dalys)
results2b_dalys         <- summary(model2b_dalys)$coefficients

#Call coefficient and CIs
summary(model2b_dalys)$coefficients[2,]
confint(model2b_dalys, "all_ncds_dalys")
#Transformed for interpretation:
exp(summary(model2b_dalys)$coefficients[2]*0.1*mean_ncd_dalys)
exp(confint(model2b_dalys, "all_ncds_dalys")*0.1*mean_ncd_dalys)

#### GRAPH ####

#Determine deaths/CFR predicted by age, sex, smoking, and healthcare_exp
model_predicted_deaths             <- glm(total_deaths ~ offset(log(total_cases)) + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_data, family=quasipoisson(), na.action=na.exclude)
COVID_NCD_data$pred_deaths         <- predict.glm(model_predicted_deaths, type="response")

# Create 'residual CFR' - difference from CFR expected from  age, sex, smoking, and healthcare_exp
COVID_NCD_data$res_cfr             <- ((COVID_NCD_data$total_deaths-COVID_NCD_data$pred_deaths)/COVID_NCD_data$total_cases)*100

#For interpretation, model relationship between NCD mortality ratio and COVID CFR
model_deaths_mortality             <- glm(total_deaths ~ offset(log(total_cases)) + ncd_mortality + m0_49 + m50_69 + m70 + f0_49 + f50_69 + f70 + smoking2018 + healthcare_exp, data=COVID_NCD_data, family=quasipoisson(), na.action=na.exclude)

#COVID_NCD_data$countryid_few <- ifelse(COVID_NCD_data$countryid == "IND" | COVID_NCD_data$countryid == "IDN" | COVID_NCD_data$countryid == "BRA" | COVID_NCD_data$countryid == "PRY" | COVID_NCD_data$countryid == "PAK" | COVID_NCD_data$countryid == "RWA" | COVID_NCD_data$countryid == "GMB", COVID_NCD_data$countryid, NA)   

ggplot(data = COVID_NCD_data) +
  theme_classic(base_size = 12,
                base_family="sans") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_text(colour="black", size=12),
        axis.text.y=element_text(colour="black", size=12),
        axis.title.x =element_text(size=14, face="bold"),
        axis.title.y =element_text(size=14, face="bold")) +
  geom_point(aes(x = ncd_mortality, y = res_cfr), colour="#960032", size=1) +
  geom_smooth(aes(x = ncd_mortality, y = res_cfr), colour="#960032", method="glm", fill="#960032", alpha=0.1, formula = y ~ exp(coefficients(model_deaths_mortality)[2]*x), size=1.2) +
  scale_x_continuous(name="NCD Annual Mortality Ratio (per 1,000)", expand=c(0,0), n.breaks=9, limits=c(2,18)) +
  scale_y_continuous(name="COVID-19 CFR (Observed - Predicted)", expand=c(0,0), n.breaks=10, limits=c(-5, 5))
