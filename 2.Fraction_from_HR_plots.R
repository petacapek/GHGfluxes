# =============================================================================
# HETEROTROPHIC FRACTION PER CANOPY (TRENCHING DATA)
# APPLIED TO THE TOTAL RESPIRATION OF THE DISTURBANCE CONTROL PLOTS
# =============================================================================
#    holisoils_soil_respiration_dataset.csv (canopy-level measurements)
#     trenched = FALSE : total respiration (TR), e.g. point "2"
#     trenched = TRUE  : heterotrophic respiration (HR), e.g. points "2A", "2B", "2C"
#
#   Only the plots in trenched_plots_to_use are used
#
#   1. HR per campaign = mean of the trenched replicates ("2A", "2B", "2C") of the plot
#   2. HR paired with the TR of the same plot and campaign 
#   3. Fraction per plot and year:   f = sum(HR) / sum(TR)
#   4. Fraction per canopy and year: mean of the plots of that canopy

#R10: Respiration normalized at 10 degrees
#R10_per_gC: Respiration normalized at 10 degrees multiplied by bulk density and carbon

#   5. Control samples (CO, 0-10 cm) of the spectra file:
#        HR10_est        = f_canopy x R10          (mg CO2-C m-2 h-1)
#        HR10_per_gC_est = f_canopy x R10_per_gC   (mg CO2-C g-1 C h-1)
#      R10 and R10_per_gC are read from the spectra file made by script GHGflux_Petr_JR.R

## =============================================================================

library(ggplot2)

# -----------------------------------------------------------------------------
# 0. SETTINGS
# -----------------------------------------------------------------------------
holisoils_file        <- "2.holisoils_soil_respiration_dataset.csv"
site_to_use           <- "Ransk"    # value of the column "siteid" for Ränskälänkorpi
years_to_use          <- c(2022, 2023)
trenched_plots_to_use <- c("2", "3", "31", "42", "52")   # CC: 2, 3 | CO: 31, 42 | CCF: 52
max_days_apart        <- 2         # HR and TR within 2 days = same campaign

# Spectra file WITH the columns R10 and R10_per_gC:
#   script 24 (Petr, group fits): "output_R10_canopy/DATA_spectra_R10_canopy.csv", spectra from column 37
#   script 15 (fits per plot)   : "output_R10/DATA_spectra_with_R10.csv", spectra from column 34
spectra_file          <- "output_R10_canopy/DATA_spectra_R10_canopy.csv"
first_spectral_column <- 38

output_dir            <- "output_HR10"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# -----------------------------------------------------------------------------
# STEP 1: READ THE HOLISOILS DATA
# -----------------------------------------------------------------------------
holi <- read.csv(holisoils_file, stringsAsFactors = FALSE, check.names = FALSE)

print(table(holi$siteid))

if (!(site_to_use %in% holi$siteid)) {
  stop("site_to_use is not a value of siteid.")
}

# !!! I added moisture and temperature data. We should check if the tr and hr have the same thermal characteristics since we know the temperature effect is significant. 
# !!! Maybe, f_canopy is temperature sensitive itself (could be according to theory). If it is, some weird fluxes we saw could be partially explained
d <- data.frame(date     = as.Date(holi$date, format = "%m/%d/%Y"),
                siteid   = holi$siteid,
                point    = as.character(holi$point),
                trenched = toupper(as.character(holi$trenched)),
                flux     = as.numeric(holi$flux),
                moisture  = as.numeric(holi$tsmoisture),
                temp      = as.numeric(holi$t05))

if (sum(is.na(d$date)) > 0) {
  stop(" Check the date format (month/day/year)")
}

d$Year <- as.numeric(format(d$date, "%Y"))
d <- d[d$siteid == site_to_use & d$Year %in% years_to_use & !is.na(d$flux), ]

# TR or HR from the column "trenched"
d$type <- NA
d$type[d$trenched == "FALSE"] <- "TR"
d$type[d$trenched == "TRUE"]  <- "HR"


# Plot name: TR point "2" -> "2";  HR point "2A" -> "2" (last letter removed)
d$plot <- d$point
d$plot[d$type == "HR"] <- sub("[A-Za-z]$", "", d$point[d$type == "HR"])


d <- d[d$plot %in% trenched_plots_to_use, ]


for (p in trenched_plots_to_use) {
  if (sum(d$plot == p & d$type == "HR") == 0) {
    stop(paste("Plot", p, "has no trenched (HR) measurements."))
  }
}

tr <- d[d$type == "TR", ]
hr <- d[d$type == "HR", ]

# -----------------------------------------------------------------------------
# STEP 2: PAIR EACH HR MEASUREMENT WITH THE TR MEASUREMENT OF ITS CAMPAIGN
# -----------------------------------------------------------------------------
hr$campaign_date <- as.Date(NA)

for (i in 1:nrow(hr)) {
  tr_dates  <- tr$date[tr$plot == hr$plot[i]]
  days_diff <- abs(as.numeric(tr_dates - hr$date[i]))
  
  if (length(days_diff) > 0) {
    if (min(days_diff) <= max_days_apart) {
      hr$campaign_date[i] <- tr_dates[which.min(days_diff)]
    }
  }
}

hr <- hr[!is.na(hr$campaign_date), ]

# Mean HR of the trenched replicated per plot and campaign
hr_mean <- aggregate(flux ~ plot + campaign_date, data = hr, FUN = mean)
names(hr_mean) <- c("plot", "date", "HR")

# Mean TR per plot and date 
tr_pair <- aggregate(flux ~ plot + date, data = tr, FUN = mean)
names(tr_pair) <- c("plot", "date", "TR")

# Only plots and dates with both TR and HR are kept
paired <- merge(tr_pair, hr_mean, by = c("plot", "date"))
paired$Year  <- as.numeric(format(paired$date, "%Y"))
paired$ratio <- paired$HR / paired$TR
paired <- paired[order(paired$plot, paired$date), ]

write.csv(paired, file.path(output_dir, "HR_TR_paired_per_date.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# STEP 3: FRACTION PER PLOT AND YEAR
# -----------------------------------------------------------------------------
fraction_plot <- data.frame()

for (p in unique(paired$plot)) {
  for (y in years_to_use) {
    rows <- paired$plot == p & paired$Year == y
    if (sum(rows) > 0) {
      fraction_plot <- rbind(fraction_plot,
                             data.frame(plot = p, Year = y,
                                        n_dates = sum(rows),
                                        sum_HR = sum(paired$HR[rows]),
                                        sum_TR = sum(paired$TR[rows]),
                                        f = sum(paired$HR[rows]) / sum(paired$TR[rows])))
    }
  }
}

# -----------------------------------------------------------------------------
# STEP 4: FRACTION PER CANOPY AND YEAR (canopy CT1 taken from the spectra file)
# -----------------------------------------------------------------------------
spec <- read.csv(spectra_file, stringsAsFactors = FALSE, check.names = FALSE)
spec$plot <- as.character(spec$plot)


fraction_plot$CT1 <- spec$CT1[match(fraction_plot$plot, spec$plot)]


fraction_canopy <- aggregate(f ~ CT1 + Year, data = fraction_plot, FUN = mean)
names(fraction_canopy) <- c("CT1", "Year", "f_canopy")


write.csv(fraction_plot, file.path(output_dir, "HR_fraction_per_plot_year.csv"), row.names = FALSE)
write.csv(fraction_canopy, file.path(output_dir, "HR_fraction_per_canopy_year.csv"), row.names = FALSE)

# Plot: HR / TR on each date, to see whether the fraction changes in the season
paired$day_of_year <- as.numeric(format(paired$date, "%j"))
paired$Year_label  <- as.character(paired$Year)

p_ratio <- ggplot(paired, aes(x = day_of_year, y = ratio, colour = Year_label)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~ plot) +
  labs(x = "Day of year", y = "HR / TR", colour = "Year",
       title = "Heterotrophic fraction on each campaign date") +
  theme_bw()

ggsave(file.path(output_dir, "HR_TR_ratio_per_date.png"), p_ratio,
       width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# STEP 5: HR10 OF THE CONTROL SAMPLES (CO, 0-10 cm) = f_canopy x R10
# -----------------------------------------------------------------------------
if (!("R10" %in% names(spec)) | !("R10_per_gC" %in% names(spec))) {
  stop("The spectra file has no R10 / R10_per_gC columns")
}

# Fraction of the canopy and year of each sample (only control samples at 0-10 cm)
row_in_fraction <- match(paste(spec$CT1, spec$year_sampling),
                         paste(fraction_canopy$CT1, fraction_canopy$Year))

is_control_0_10 <- spec$ST1 == "CO" & spec$Depth == "0-10"
row_in_fraction[!is_control_0_10] <- NA

f_used <- fraction_canopy$f_canopy[row_in_fraction]

new_columns <- data.frame(f_canopy        = f_used,
                          HR10_est        = f_used * spec$R10,          # mg CO2-C m-2 h-1
                          HR10_per_gC_est = f_used * spec$R10_per_gC)   # mg CO2-C g-1 C h-1

# Table of the control samples, to check
check <- data.frame(spec[, c("year_sampling", "plot", "CT1", "ST1", "Depth", "R10", "R10_per_gC")],
                    new_columns)
check <- check[!is.na(check$f_canopy), ]
check <- check[order(check$year_sampling, check$CT1, check$plot), ]

write.csv(check, file.path(output_dir, "HR10_control_plots.csv"), row.names = FALSE)

# Metadata | new columns | spectra
spectra_out <- cbind(spec[, 1:(first_spectral_column - 1)],
                     new_columns,
                     spec[, first_spectral_column:ncol(spec)])

output_spectra_file <- sub(".csv", "_HR10.csv", basename(spectra_file), fixed = TRUE)
write.csv(spectra_out, file.path(output_dir, output_spectra_file), row.names = FALSE)


# -----------------------------------------------------------------------------
# Temperature sensitivities of heterotrophic and total respiration
# -----------------------------------------------------------------------------
# Here I paire your dataframes across all observations to see the entire variability
paired_all <- merge(tr[, c('flux', 'temp', 'plot', 'moisture', 'Year','date')],
                    hr[, c('flux', 'temp', 'plot', 'moisture',  'campaign_date')],
                    by.x = c('plot','date'),
                    by.y = c('plot', 'campaign_date'))
cnames <- gsub('.x', '.total', names(paired_all), fixed = T)
cnames <- gsub('.y', '.hetero', cnames, fixed = T)

colnames(paired_all) <- cnames

# I will add a new column defining canopy treatments (CC: 2, 3 | CO: 31, 42 | CCF: 52)
paired_all <- paired_all |> mutate(Canopy_treatment = case_when(plot == 2 | plot == 3 ~ 'CC',
                                                                plot == 31 | plot == 42 ~ 'CO',
                                                                plot == 52 ~ 'CCF'))
# Now, I will plot the two fluxes against each other and I need to see some relationship below 1:1 line
paired_all |> filter(flux.total > 0 & flux.hetero > 0) |> 
  ggplot(aes(flux.total, flux.hetero)) + geom_point(cex = 6, pch = 21, fill = 'grey') +
  theme_test() + facet_grid(Year ~ plot, scales = 'free') +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(method = 'lm', se = F) 

# Print the slopes of the relationship between heterotrophic and total respiration for each treatment and year
for(i in unique(paired_all$Canopy_treatment)){
  for(j in c(2022, 2023)){
    print(
      c(j, i,
        coef(
          lm(flux.hetero ~ flux.total, data = paired_all, subset = Canopy_treatment == i & Year == j)
        )[[2]]
        )
    )
  }
}

# I will add another column identifying positive, none or negative difference in background temperatures of measurement
paired_all <- paired_all |> mutate(temp_diff = case_when(temp.hetero - temp.total > 0 ~ 'positive',
                                                         temp.hetero - temp.total < 0 ~ 'negative',
                                                         temp.hetero - temp.total == 0 ~ 'equal'))

paired_all |> filter(flux.total > 0 & flux.hetero > 0) |> 
  ggplot(aes(flux.total, flux.hetero)) + geom_point(cex = 6, pch = 21, aes(fill = temp_diff)) +
  theme_test() + facet_grid(Year ~ Canopy_treatment, scales = 'free') +
  geom_abline(intercept = 0, slope = 1) +
  stat_smooth(method = 'lm', se = F) 

paired_all |> filter(flux.total > 0 & flux.hetero > 0) |> 
  ggplot(aes(temp.hetero, flux.hetero)) + geom_point(cex = 6, pch = 21) +
  theme_test() + facet_grid( ~ Canopy_treatment, scales = 'free') 

paired_all |> filter(flux.total > 0 & flux.hetero > 0) |> 
  ggplot(aes(temp.total, flux.total)) + geom_point(cex = 6, pch = 21) +
  theme_test() + facet_grid(Year ~ Canopy_treatment, scales = 'free') 
