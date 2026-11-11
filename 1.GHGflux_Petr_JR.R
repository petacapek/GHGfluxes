###################################################################################
############################Greenhouse gas fluxes##################################
#################analysis of temperature and moisture effects######################
##########################Ränskälänkorpi study site################################
###################################################################################
#
# Petr's script, unchanged, with ADDED parts (marked "ADDED"):
#   - CO2 flux converted to mg CO2-C (x 12/44) for the new fits
#   - linear temperature function next to the exponential one (better R2 kept)
#   - flux normalised at 10, 15 and 20 deg C (R10, R15, R20) per group
#   - plots of both fits for each year x canopy x disturbance group
#   - R10/R15/R20 per gram of soil C (C% and bulk density of each sample)
###################################################################################

# ---- The required R libraries are listed in separate document and loaded here ----
# source("packages.R")   # Petr's file with the libraries; the ones needed are loaded here:

# ---- Data files are here ----
GHG <- read.csv('1.GHG flux.csv')
str(GHG)
##We are interested in 2022 and 2023
GHG <- subset(GHG, Year > 2021 & Year < 2024)
summary(GHG)

# ---- FTIR data ----
## the data will be used to first compare treatment levels and organization with GHG data frame
## when everything works well, normalized GHG fluxes will be attached to this file.
FTIR <- read.csv('3.DATA_spectra.csv', check.names = FALSE)  #CHANGED FROM XLSX TO CSV
str(FTIR[, 1:24])

# ADDED: settings for the new parts
first_spectral_column <- 25            # spectra start here in DATA_spectra.csv
c_column              <- "C"           # total carbon (%)
bd_column             <- "bd_column"   # bulk density (g cm-3)
layer_cm              <- 10            # thickness of the soil layer (cm)
p_limit               <- 0.05          # significance limit of the temperature fit
output_dir            <- "output_R10_canopy"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

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
                       summary(mT)$coefficients[7], #p_temp   (NOTE: [7] is the p of the intercept; [8] is the p of T5)
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
                     summary(mT)$coefficients[7], #p_temp   (NOTE: [7] is the p of the intercept; [8] is the p of T5)
                     NA, NA, #s1 and s2
                     NA, NA, #p_s1 and p_s2
                     summary(mT)$r.squared) #R2
        CO2_res[(CO2_res$Year == i & CO2_res$Canopy_treatment == j & CO2_res$Disturbance_treatment == k),
                c(4:ncol(CO2_res))] <- res_out
        
      }, error = function(e){print(paste("Something is wrong with:", c(i, j, k)))})
    }
  }
}

# =============================================================================
# ADDED: CONVERSION, EXPONENTIAL AND LINEAR FIT, R10 / R15 / R20 PER GROUP
# =============================================================================
# convert from CO2flux_mg_m2_h to Flux in mg CO2-C m-2 h-1 (x 12/44)
GHG$CO2flux_mgC_m2_h <- GHG$CO2flux_mg_m2_h * 12 / 44

# For each Year x Canopy x Disturbance group (4 plots pooled, as above):
#   exponential: log(flux) = log(A) + a*T   ->  flux = exp(log(A) + a*T) = A*exp(a*T)
#   linear     : flux = b0 + b1*T
# Both fits use the same rows (flux > 0 and T5 present).
# R2 of both on the original flux scale; the higher R2 is kept.
# If the kept fit is not significant (p of T5 >= p_limit), the mean flux is R10.

CO2_res$n      <- NA
CO2_res$R2_exp <- NA
CO2_res$p_exp  <- NA
CO2_res$R2_lin <- NA
CO2_res$p_lin  <- NA
CO2_res$chosen <- NA
CO2_res$R10    <- NA
CO2_res$R15    <- NA
CO2_res$R20    <- NA

curves <- data.frame()   # fitted curves for the plots

for (g in 1:nrow(CO2_res)) {
  
  # Measurements of this group
  d <- GHG[GHG$Year == CO2_res$Year[g] &
             GHG$Canopy_treatment == CO2_res$Canopy_treatment[g] &
             GHG$Disturbance_treatment == CO2_res$Disturbance_treatment[g] &
             !is.na(GHG$CO2flux_mgC_m2_h) & !is.na(GHG$T5) & GHG$CO2flux_mgC_m2_h > 0, ]
  
  if (nrow(d) < 5) {
    next
  }
  
  m_exp <- lm(log(CO2flux_mgC_m2_h) ~ T5, data = d)
  m_lin <- lm(CO2flux_mgC_m2_h ~ T5, data = d)
  
  # R2 on the original flux scale
  ss_tot <- sum((d$CO2flux_mgC_m2_h - mean(d$CO2flux_mgC_m2_h))^2)
  r2_exp <- 1 - sum((d$CO2flux_mgC_m2_h - exp(fitted(m_exp)))^2) / ss_tot
  r2_lin <- 1 - sum((d$CO2flux_mgC_m2_h - fitted(m_lin))^2) / ss_tot
  
  p_exp <- summary(m_exp)$coefficients["T5", "Pr(>|t|)"]
  p_lin <- summary(m_lin)$coefficients["T5", "Pr(>|t|)"]
  
  # Better fit; mean flux if it is not significant
  if (r2_exp >= r2_lin) {
    chosen <- "exponential"
    p_chosen <- p_exp
  } else {
    chosen <- "linear"
    p_chosen <- p_lin
  }
  if (p_chosen >= p_limit) {
    chosen <- "none (mean flux)"
  }
  
  # Flux at 10, 15 and 20 deg C (mg CO2-C m-2 h-1)
  # exponential: exp(intercept + slope * T) = A * exp(a * T)
  # linear     : intercept + slope * T
  if (chosen == "exponential") {
    R <- exp(coef(m_exp)["(Intercept)"] + coef(m_exp)["T5"] * c(10, 15, 20))
  } else if (chosen == "linear") {
    R <- coef(m_lin)["(Intercept)"] + coef(m_lin)["T5"] * c(10, 15, 20)
  } else {
    R <- c(mean(d$CO2flux_mgC_m2_h), NA, NA)
  }
  
  CO2_res$n[g]      <- nrow(d)
  CO2_res$R2_exp[g] <- round(r2_exp, 3)
  CO2_res$p_exp[g]  <- signif(p_exp, 3)
  CO2_res$R2_lin[g] <- round(r2_lin, 3)
  CO2_res$p_lin[g]  <- signif(p_lin, 3)
  CO2_res$chosen[g] <- chosen
  CO2_res$R10[g]    <- R[1]
  CO2_res$R15[g]    <- R[2]
  CO2_res$R20[g]    <- R[3]
  
  # Curves for the plots
  t_grid <- seq(min(d$T5), max(d$T5), length.out = 50)
  curves <- rbind(curves,
                  data.frame(Year = CO2_res$Year[g],
                             Canopy_treatment = CO2_res$Canopy_treatment[g],
                             Disturbance_treatment = CO2_res$Disturbance_treatment[g],
                             T5 = t_grid,
                             exponential = exp(coef(m_exp)["(Intercept)"] + coef(m_exp)["T5"] * t_grid),
                             linear = coef(m_lin)["(Intercept)"] + coef(m_lin)["T5"] * t_grid))
}

print(CO2_res)
write.csv(CO2_res, file.path(output_dir, "CO2_fits_per_group.csv"), row.names = FALSE)

# =============================================================================
# ADDED: PLOTS OF THE EXPONENTIAL (red) AND LINEAR (blue) FITS, one figure per year
# =============================================================================
# Points coloured by plot, to see the 4 plots pooled in each group
if (!("plot" %in% names(GHG))) {
  GHG$plot <- GHG$Location
}

for (y in c(2022, 2023)) {
  
  points_y <- GHG[GHG$Year == y & !is.na(GHG$CO2flux_mgC_m2_h) & !is.na(GHG$T5) &
                    GHG$CO2flux_mgC_m2_h > 0, ]
  curves_y <- curves[curves$Year == y, ]
  res_y    <- CO2_res[CO2_res$Year == y & !is.na(CO2_res$chosen), ]
  
  points_y$panel <- paste(points_y$Canopy_treatment, points_y$Disturbance_treatment)
  curves_y$panel <- paste(curves_y$Canopy_treatment, curves_y$Disturbance_treatment)
  res_y$panel    <- paste(res_y$Canopy_treatment, res_y$Disturbance_treatment)
  
  panel_titles <- paste0(res_y$panel, ": ", res_y$chosen,
                         "\nR2 exp ", res_y$R2_exp, " | R2 lin ", res_y$R2_lin)
  names(panel_titles) <- res_y$panel
  
  p_fit <- ggplot() +
    geom_point(data = points_y, aes(x = T5, y = CO2flux_mgC_m2_h, colour = as.character(plot)),
               size = 1.5, alpha = 0.7) +
    geom_line(data = curves_y, aes(x = T5, y = exponential), colour = "red", linewidth = 1) +
    geom_line(data = curves_y, aes(x = T5, y = linear), colour = "blue", linewidth = 1) +
    facet_wrap(~ panel, scales = "free_y", ncol = 3, labeller = labeller(panel = panel_titles)) +
    labs(x = "Soil temperature T5 (deg C)", y = "CO2 flux (mg CO2-C m-2 h-1)", colour = "Plot",
         title = paste("CO2 flux vs T5 per group,", y, "(red = exponential, blue = linear)")) +
    theme_bw() +
    theme(strip.text = element_text(size = 8))
  
  ggsave(file.path(output_dir, paste0("fits_per_group_", y, ".png")),
         p_fit, width = 13, height = 10, dpi = 300)
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

# =============================================================================
# ADDED: R10 / R15 / R20 OF THE GROUP AND PER GRAM OF SOIL CARBON FOR EACH SAMPLE
# =============================================================================
# Same depth rule as above: 20-30 cm for Topsoil_removal, 0-10 cm for the others
use_row <- (FTIR_full$Disturbance_treatment == "Topsoil_removal" & FTIR_full$Depth == "20-30") |
  (FTIR_full$Disturbance_treatment != "Topsoil_removal" & FTIR_full$Depth == "0-10")

# Row of CO2_res (same year, canopy and disturbance) for each sample
row_in_res <- match(paste(FTIR_full$year_sampling, FTIR_full$Canopy_treatment, FTIR_full$Disturbance_treatment),
                    paste(CO2_res$Year, CO2_res$Canopy_treatment, CO2_res$Disturbance_treatment))
row_in_res[!use_row] <- NA

FTIR_full$CO2_fit <- CO2_res$chosen[row_in_res]
FTIR_full$R10     <- CO2_res$R10[row_in_res]     # mg CO2-C m-2 h-1
FTIR_full$R15     <- CO2_res$R15[row_in_res]
FTIR_full$R20     <- CO2_res$R20[row_in_res]

# Carbon stock of each sample (g C m-2) = C% / 100 x BD (g cm-3) x 10 cm x 10000 cm2 m-2
FTIR_full$C_stock_g_m2 <- as.numeric(FTIR_full[[c_column]]) / 100 *
  as.numeric(FTIR_full[[bd_column]]) * layer_cm * 10000

# Per gram of soil carbon (mg CO2-C g-1 C h-1)
FTIR_full$R10_per_gC <- FTIR_full$R10 / FTIR_full$C_stock_g_m2
FTIR_full$R15_per_gC <- FTIR_full$R15 / FTIR_full$C_stock_g_m2
FTIR_full$R20_per_gC <- FTIR_full$R20 / FTIR_full$C_stock_g_m2

#======================exporting to csv
write.csv(FTIR_full, row.names = F, 'DATA_normalizedFluxes.csv')

# =============================================================================
# ADDED: SAME DATA WITH THE NEW COLUMNS BEFORE THE SPECTRA 
# =============================================================================
# In DATA_normalizedFluxes.csv the new columns are AFTER the spectra.
# Here: metadata | new columns | spectra
n_original  <- ncol(FTIR)
new_names   <- names(FTIR_full)[(n_original + 1):ncol(FTIR_full)]

FTIR_ordered <- cbind(FTIR_full[, 1:(first_spectral_column - 1)],
                      FTIR_full[, new_names],
                      FTIR_full[, first_spectral_column:n_original])

write.csv(FTIR_ordered, row.names = F, file.path(output_dir, 'DATA_spectra_R10_canopy.csv'))
