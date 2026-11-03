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
  labs(x = 'Temperature (˚C)', y = expression(paste(CO[2]~flux~(mg~C-CO[2]~m^{-2}~h^{-1}))), fill = 'Year')
####There is a relationship - seems to differ between among years and treatments

##4. Moisture dependence of CO2 flux
ggplot(GHG, aes(WT, CO2flux_mg_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Moisture (unknown units)', y = expression(paste(CO[2]~flux~(mg~C-CO[2]~m^{-2}~h^{-1}))), fill = 'Year')
####The relationship is not straightforward but certain tendencies are subjectively present

# ---- Fitting temperature and moisture response functions ----
#Temperature function is straightforward - exponential function
##                                          flux = A*exp(a*T)
##                                          in which T is temperature in 5cm depth, a is temperature sensitivity (exp(a*10) defines Q10 parameter), and A is normalization constant in mg of CO2 per m^2 and hour - this defines the temperature independent flux
##                                         the function can be linearized - ln(flux) = ln(A) - a*T and A estimated
##                                         !!!ln(flux) does not work with negative numbers, which does not make sense anyway. The negative fluxes will be ignored

#Moisture function is less obvious - it has unimodal character in theory - I will use 2nd order polynomial, but test its significance

#Results will be stored here (one estimate per year and unique treatment combination):
CO2_res <- as.data.frame(GHG |> group_by(Year, Canopy_treatment, Disturbance_treatment) |> 
                           summarise(A = NA, #normalization constant
                                     Ase = NA, #standard error of estimate - I will need it to generate 4 replicates by sampling from random distribution
                                     a = NA, #temperature sensitivity 
                                     ase = NA, #standard error
                                     p_temp = NA, #p for temperature response function
                                     s1 = NA, #first order slope to soil moisture
                                     s2 = NA, #second order slope to soil moisture
                                     p_s1 = NA, #p for moisture response function - first polynomial term
                                     p_s2 = NA, #p for moisture response function - second polynomial term
                                     R2 = NA #coefficient of determination of the final model
                                     ))

#The for loop is used to do all the work
for(i in unique(GHG$Year)){
  for(j in unique(GHG$Canopy_treatment)){
    for(k in unique(GHG$Disturbance_treatment)){
      tryCatch({#tryCatch is specified here to skip potential errors in the loop
        
        #Temperature effect
        mT <- lm(log(CO2flux_mg_m2_h) ~ T5, data = GHG, subset = (Year == i & Canopy_treatment == j & Disturbance_treatment == k))
        #Moisture effect
        mWT <- update(mT, . ~ . + poly(WT, 2, raw = T))
        #if else loop is used to extract coefficients from mT or mWT depending on the statistical significance of the polynomial term
        if(summary(mWT)$coefficients[15] < 0.05 | summary(mWT)$coefficients[16] < 0.05){
          #this is a vector with coefficients of mWT model
          res_out <- c(exp(summary(mWT)$coefficients[1]), #A
                       exp(summary(mWT)$coefficients[5]), #Ase
                       summary(mWT)$coefficients[2], #a
                       summary(mWT)$coefficients[6], #ase
                       summary(mWT)$coefficients[14], #p_temp
                       summary(mWT)$coefficients[c(3,4)], #s1 and s2
                       summary(mWT)$coefficients[c(15, 16)], #p_s1 and p_s2
                       summary(mWT)$r.squared) #R2
          CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & CO2_res$Disturbance_treatment == k),
                  c(4:ncol(CO2_res))] <- res_out
          
        } else {
          res_out <- c(exp(summary(mT)$coefficients[1]), #A
                       exp(summary(mT)$coefficients[3]), #Ase
                       summary(mT)$coefficients[2], #a
                       summary(mT)$coefficients[4], #ase
                       summary(mT)$coefficients[7], #p_temp
                       NA, NA, #s1 and s2
                       NA, NA, #p_s1 and p_s2
                       summary(mT)$r.squared) #R2
          CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & CO2_res$Disturbance_treatment == k),
                  c(4:ncol(CO2_res))] <- res_out
        }
        
      }, error = function(e){print(paste("Something is wrong with:", c(i, j, k)))})
    }
  }
}

####Notice that all temperature effects are highly significant.
####The numbers make sense - CO2 flux increases with residue addition and decreases with topsoil removal
####That brings me to one important consideration. The CO2 flux scales with organic C stocks so I suggest to normalize flux to a unit of organic C so additional variability goes to organic C quality and thus, FTIR
####Moisture effect is significant in 3 out of 18 cases so I will refit all data with temperature function only
for(i in unique(GHG$Year)){
  for(j in unique(GHG$Canopy_treatment)){
    for(k in unique(GHG$Disturbance_treatment)){
      tryCatch({#tryCatch is specified here to skip potential errors in the loop
        
        #Temperature effect
        mT <- lm(log(CO2flux_mg_m2_h) ~ T5, data = GHG, subset = (Year == i & Canopy_treatment == j & Disturbance_treatment == k))
       
        res_out <- c(exp(summary(mT)$coefficients[1]), #A
                     exp(summary(mT)$coefficients[3]), #Ase
                     summary(mT)$coefficients[2], #a
                     summary(mT)$coefficients[4], #ase
                     summary(mT)$coefficients[7], #p_temp
                     NA, NA, #s1 and s2
                     NA, NA, #p_s1 and p_s2
                     summary(mT)$r.squared) #R2
        CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & CO2_res$Disturbance_treatment == k),
                c(4:ncol(CO2_res))] <- res_out
       
      }, error = function(e){print(paste("Something is wrong with:", c(i, j, k)))})
    }
  }
}

#========================CH4 flux
##1. Temperature dependence of CH4 flux
ggplot(GHG, aes(T5, CH4flux_ug_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Temperature (˚C)', y = expression(paste(CH[4]~flux~(mu*g~C-CH[4]~m^{-2}~h^{-1}))), fill = 'Year')
####There is no apparent relationship

##2. Moisture dependence of CH4 flux
ggplot(GHG, aes(WT, CH4flux_ug_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Moisture (unknown units)', y = expression(paste(CH[4]~flux~(mu*g~C-CH[4]~m^{-2}~h^{-1}))), fill = 'Year')
####There is no apparent relationship
####It does not seem necessary to normalize fluxes. It is, however, important to state in manuscript that all fluxes are negative - methane is consumed by methanotrophs - methane is source of energy and C - it should be FTIR independent

#========================N2O flux
##1. Temperature dependence of N2O flux
ggplot(GHG, aes(T5, N2Oflux_ug_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Temperature (˚C)', y = expression(paste(N[2]*O~flux~(mu*g~N-N[2]*O~m^{-2}~h^{-1}))), fill = 'Year')
####There is no apparent relationship

##2. Moisture dependence of N2O flux
ggplot(GHG, aes(WT, N2Oflux_ug_m2_h)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Moisture (unknown units)', y = expression(paste(N[2]*O~flux~(mu*g~N-N[2]*O~m^{-2}~h^{-1}))), fill = 'Year')
####There is no apparent relationship
####It does not seem necessary to normalize N2O fluxes too.

################################################################################
###########################Merging with FTIR data###############################
################################################################################
is.factor(FTIR$CT1)
FTIR$CT1 <- as.factor(FTIR$CT1)
#I need two new columns with factor names and levels corresponding to GHG dataframe. I hope I did not mess up
FTIR_full <- FTIR |> mutate(Canopy_treatment = factor(CT1, levels = c('CC', 'CO', 'MT'),
                                                      labels = unique(GHG$Canopy_treatment)),
                            Disturbance_treatment = factor(ST1, levels = c('CO', 'TLR', 'M'),
                                                           labels = unique(GHG$Disturbance_treatment)))
#I will add normalized CO2 flux in mg C-CO2 per m^2 and hour (i.e. temperature effect filtered out) and Q10 (estimate ot temperature sensitivity)
FTIR_full$CO2_normalized <- numeric(length = nrow(FTIR))
FTIR_full$Q10_normalized <- numeric(length = nrow(FTIR))

#for loop will do the job again
for(i in unique(GHG$Year)){
  for(j in unique(GHG$Canopy_treatment)){
    for(k in unique(GHG$Disturbance_treatment)){
      Amean = CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & 
                         CO2_res$Disturbance_treatment == k), "A"]
      Ase = CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & 
                       CO2_res$Disturbance_treatment == k), "Ase"]
      amean = CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & 
                         CO2_res$Disturbance_treatment == k), "a"]
      ase = CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & 
                       CO2_res$Disturbance_treatment == k), "ase"]
      #I will add the numbers to corresponding depths - 20 - 30 in case of TLR or 0 - 10 in case of CO and M
      if(k == "Topsoil_removal"){
        #The respective subset should have 4 replicates but I do not want to risk anything so I define nrep
        nrep = nrow(FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                                 FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '20-30'), ])
        
        #Generate normalized CO2 flux by sampling from random normal distribution
        CO2draws <- rnorm(nrep, mean = Amean, sd = Ase)
        #Generate Q10 by sampling from random normal distribution
        Q10draws <- rnorm(nrep, mean = amean, sd = ase)
        
        #Paste these to respective FTIR subset
        FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                     FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '20-30'), 'CO2_normalized'] <- CO2draws
        FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                     FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '20-30'), 'Q10_normalized'] <- exp(Q10draws*10)
      } else {
        #The respective subset should have 4 replicates but I do not want to risk anything so I define nrep
        nrep = nrow(FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                                 FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '0-10'), ])
        
        #Generate normalized CO2 flux by sampling from random normal distribution
        CO2draws <- rnorm(nrep, mean = Amean, sd = Ase)
        #Generate Q10 by sampling from random normal distribution
        Q10draws <- rnorm(nrep, mean = amean, sd = ase)
        
        #Paste these to respective FTIR subset
        FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                     FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '0-10'), 'CO2_normalized'] <- CO2draws
        FTIR_full[(FTIR_full$year_sampling == i & FTIR_full$Canopy_treatment == j & 
                     FTIR_full$Disturbance_treatment == k & FTIR_full$Depth == '0-10'), 'Q10_normalized'] <- exp(Q10draws*10)
      }
    }
  }
}

FTIR_full$CO2_normalized[FTIR_full$CO2_normalized == 0] <- NA
FTIR_full$Q10_normalized[FTIR_full$Q10_normalized == 0] <- NA

#======================exporting to csv
write.csv(FTIR_full, row.names = F, 'DATA_normalizedFluxes.csv')
