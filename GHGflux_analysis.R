###################################################################################
############################Greenhouse gas fluxes##################################
#################analysis of temperature and moisture effects######################
##########################Ränskälänkorpi study site################################
###################################################################################

# ---- The required R libraries are listed in separate document and loaded here ---- 
source("packages.R")

# ---- Data files are here ----
GHG <- read.csv('GHG flux.csv')
str(GHG)
##We are interested in 2022 and 2023
GHG <- subset(GHG, Year > 2021 & Year < 2024)
summary(GHG)

# ---- FTIR data ----
## the data will be used to first compare treatment levels and organization with GHG data frame
## when everything works well, normalized GHG fluxes will be attached to this file.
FTIR <- read.xlsx('Top_Spectra_selected_and_variables130_JR.xlsx')
str(FTIR)

# ---- GHG and FTIR data comparison ----
##Number of unique plots
unique(GHG$Plot_code)
unique(GHG$Location)
unique(FTIR$plot)
### none of plot identifiers match

##Treatments - canopy (3 levels) + disturbance (3 levels) - 4 replicates
unique(GHG$Canopy_treatment); unique(GHG$Disturbance_treatment)
unique(FTIR$CT1); unique(FTIR$ST1)

###Number of values for each combination of treatments divided by 4 (replicates)
GHG |> group_by(Year, Canopy_treatment, Disturbance_treatment) |> filter(!is.na(CO2flux_mg_m2_h)) |> 
  summarise(nobs = length(unique(CO2flux_mg_m2_h))/4)
####approximately 12 measurements for each site
FTIR |> group_by(year_sampling, CT1, ST1) |> filter(year_sampling > 2021 & !is.na(CO2_g_m2_yearly)) |> 
  summarise(nobs = length(unique(CO2_g_m2_yearly))/4)
####it should be 1 

####################################ALL LOOK GOOD - DATA CHECKING END HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

################################################################################
################################DATA ANALYSIS###################################
################################################################################
#========================CO2 flux
# ---- Data visualization (presence/absence of the relationship to temperature/moisture) ----
##1. histogram of temperatures in 5 cm depth + table with means, median and quartiles
ggplot(GHG, aes(x = T5)) + geom_histogram(binwidth = 0.9, aes(fill = factor(Year)), alpha = 0.8) + theme_test() +
  facet_grid(Canopy_treatment ~ Disturbance_treatment) +
  labs(x = 'Temperature (˚C)', fill = 'Year')

GHG |> group_by(Year, Canopy_treatment, Disturbance_treatment) |> 
  summarise(Center = mean(T5, na.rm = T), Median = median(T5, na.rm = T), 
            Q75 = quantile(T5, probs = c(0.75), na.rm = T), 
            Q25 = quantile(T5, probs = c(0.25), na.rm = T),
            nmissing = length(is.na(T5)))
####2023 is warmer than 2022 by almost 2 degrees and there are differences between treatments (life forest is unsurprisingly colder)

##2. histogram of soil moistures - I do not know what units WT has  
ggplot(GHG, aes(x = WT)) + geom_histogram(binwidth = 0.9, aes(fill = factor(Year)), alpha = 0.8) + theme_test() +
  facet_grid(Canopy_treatment ~ Disturbance_treatment) +
  labs(x = 'Soil moisture (uknown units)', fill = 'Year')

GHG |> group_by(Year, Canopy_treatment, Disturbance_treatment) |> 
  summarise(Center = mean(WT, na.rm = T), Median = median(WT, na.rm = T), 
            Q75 = quantile(WT, probs = c(0.75), na.rm = T), 
            Q25 = quantile(WT, probs = c(0.25), na.rm = T),
            nmissing = length(is.na(WT)))
####Differences in soil moisture are striking - especially between canopy treatments.

##3. Temperature dependence of CO2 flux
ggplot(GHG, aes(T5, CO2flux_mg_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Temperature (˚C)', y = expression(paste(CO[2]~flux~(mg~C-CO[2]~m^{-2}~h^{-2}))), fill = 'Year')
####There is a relationship - seems to differ between among years and treatments

##4. Moisture dependence of CO2 flux
ggplot(GHG, aes(WT, CO2flux_mg_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Moisture (unknown units)', y = expression(paste(CO[2]~flux~(mg~C-CO[2]~m^{-2}~h^{-2}))), fill = 'Year')
####The relationship is not straightforward but certain tendencies are subjectively present

# ---- Fitting temperature and moisture response functions ----
#Temperature function is straightforward - Arrhenius equation:
##                                          flux = A*exp(-Ea/R/T)
##                                          in which T is temperature in K, Ea is activation energy in J/mol, R is gas constant (8.314 J/K/mol), and A is normalization constant in mg of CO2 per m^2 and hour - this defines the temperature independent flux
##                                         the function can be linearized - ln(flux) = ln(A) - Ea*(1/R/T) and A estimated

#Moisture function is less obvious - it has unimodal character in theory - I will use 2nd order polynomial, but test its significance

#Results will be stored here (one estimate per year and unique treatment combination):
CO2_res <- as.data.frame(GHG |> group_by(Year, Canopy_treatment, Disturbance_treatment) |> 
                           summarise(A = NA, #normalization constant
                                     Ea = NA, #activation energy
                                     p_temp = NA, #p for temperature response function
                                     s1 = NA, #first order slope to soil moisture
                                     s2 = NA, #second order slope to soil moisture
                                     p_moisture = NA, #p for moisture response function
                                     R2 = NA #coefficient of determination of the final model
                                     ))

#Temperature is converted to 1/R/T
GHG$invRT <- 1/8.314/(GHG$T5 + 273.15)

#The for loop is used to do all the work
for(i in unique(GHG$Year)){
  for(j in unique(GHG$Canopy_treatment)){
    for(k in unique(GHG$Disturbance_treatment)){
      tryCatch({
        
      }, error = function(e){print(paste("Something is wrong with:", c(i, j, k)))})
    }
  }
}