# =============================================================================
# 7. HETEROTROPHIC RESPIRATION FROM A TOTAL RESPIRATION AND TEMPERATURE MODEL
# =============================================================================
# The fraction of heterotrophic respiration is a logistic function of a temperature
# f(T) = Rh / (Rh + Ra) = 1 / (1 + exp(−(α + β·(T − 10)/10)))
# where 
# 
# α = ln(Rh10 / Ra10) is the log ratio of the two fluxes at 10 °C. At 10 °C, f = plogis(α).
# β = ln(Q10_h / Q10_a) describes how f changes with temperature. If β = 0, f doesn’t depend on temperature.
# 
# logistic curve is convenient because it is always between 0 and 1, exactly what you expect from fraction of heterotrophic respiration
#
# Knowing parameters α, β and temperature, f(T) is fully defined and heterotrophic respiration (HR) can be calculated as: 
# HR = TR × plogis(α + β·(T − 10)/10)
#
# =============================================================================
source("packages.R")

# Settings from the main script; defaults are used if the section is run alone
if (!exists("holisoils_file"))        holisoils_file        <- "2.holisoils_soil_respiration_dataset.csv"
if (!exists("site_to_use"))           site_to_use           <- "Ransk"
if (!exists("years_to_use"))          years_to_use          <- c(2022, 2023)
if (!exists("trenched_plots_to_use")) trenched_plots_to_use <- c("2", "3", "31", "42", "52")
if (!exists("max_days_apart"))        max_days_apart        <- 2
if (!exists("output_dir"))            output_dir            <- "output_HR10"
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!file.exists(holisoils_file)) {
  stop("Data file not found: ", holisoils_file, "  (working directory: ", getwd(), ")")
}

canopy_of_plot  <- c("2" = "CC", "3" = "CC", "31" = "CO", "42" = "CO", "52" = "CCF")

max_ci          <- 0.01          # drop chamber fits with ci above this (NULL = keep all)
max_temp_offset <- 5             # drop t05 deviating > 5 °C from the site median of that date
q10_short       <- 2.5           # short-term Q10 used for the main result
q10_short_range <- c(1, 1.5, 2, 2.5, 3, 4)   # sensitivity (1 = no correction)
q10a_prior      <- c(mean = log(2), sd = 0.5) # weak prior on log Q10 of Ra (6.4); NULL = none
n_boot          <- 500
set.seed(1)

q95 <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)

# -----------------------------------------------------------------------------
# 7.1 DATA AND CLEANING
# -----------------------------------------------------------------------------
raw <- read.csv(holisoils_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Exact duplicate rows removed:", sum(duplicated(raw)), "\n")
raw <- unique(raw)

m <- data.frame(date     = as.Date(raw$date, format = "%m/%d/%Y"),
                siteid   = raw$siteid,
                point    = as.character(raw$point),
                trenched = toupper(as.character(raw$trenched)),
                flux     = as.numeric(raw$flux),
                ci       = as.numeric(raw$ci),
                temp     = as.numeric(raw$t05))
m$Year <- as.numeric(format(m$date, "%Y"))
m <- m[m$siteid == site_to_use & m$Year %in% years_to_use, ]

m$type <- ifelse(m$trenched == "TRUE", "HR", "TR")
m$plot <- m$point
m$plot[m$type == "HR"] <- sub("[A-Za-z]$", "", m$point[m$type == "HR"])
m <- m[m$plot %in% trenched_plots_to_use, ]
m$canopy <- unname(canopy_of_plot[m$plot])

# quality flags
m$flag <- ""
m$flag[is.na(m$flux) | is.na(m$temp)] <- "missing"
if (!is.null(max_ci)) m$flag[m$flag == "" & m$ci > max_ci] <- "high_ci"
site_med <- ave(m$temp, m$date, FUN = function(x) median(x, na.rm = TRUE))
m$flag[m$flag == "" & abs(m$temp - site_med) > max_temp_offset] <- "temp_outlier"

cat("\nQuality flags (rows):\n")
print(table(m$flag, m$type))
write.csv(m[m$flag != "", ], file.path(output_dir, "removed_rows.csv"), row.names = FALSE)
m <- m[m$flag == "", ]

# -----------------------------------------------------------------------------
# 7.2 PAIRING: HR mean of replicates (with its mean temperature) + TR, per plot and campaign
# -----------------------------------------------------------------------------
tr <- m[m$type == "TR", ]
hr <- m[m$type == "HR", ]

hr$campaign_date <- as.Date(NA)
for (i in seq_len(nrow(hr))) {
  td <- tr$date[tr$plot == hr$plot[i]]
  dd <- abs(as.numeric(td - hr$date[i]))
  if (length(dd) > 0 && min(dd) <= max_days_apart) hr$campaign_date[i] <- td[which.min(dd)]
}
hr <- hr[!is.na(hr$campaign_date), ]

hr_a <- aggregate(cbind(HR = flux, T_HR = temp) ~ plot + canopy + campaign_date,
                  data = hr, FUN = mean)
names(hr_a)[names(hr_a) == "campaign_date"] <- "date"
tr_a <- aggregate(cbind(TR = flux, T_TR = temp) ~ plot + date, data = tr, FUN = mean)

pr <- merge(tr_a, hr_a, by = c("plot", "date"))
pr$Year <- as.numeric(format(pr$date, "%Y"))
pr$dT   <- pr$T_TR - pr$T_HR
pr <- pr[order(pr$canopy, pr$plot, pr$date), ]
if (nrow(pr) == 0) stop("No paired HR-TR campaigns found - check site_to_use, years and plots.")

cat("\nPaired plot-campaigns:", nrow(pr), "\n")
cat("dT = T_TR - T_HR: median", round(median(pr$dT), 2), "°C, range",
    round(min(pr$dT), 1), "to", round(max(pr$dT), 1), "\n")

# -----------------------------------------------------------------------------
# 7.3 FRACTION WITH SHORT-TERM TEMPERATURE CORRECTION (Step 3/4 logic)
#     f_plot = sum(HR_corr) / sum(TR) over the campaigns of a year
#     f_canopy = mean of the plots
# -----------------------------------------------------------------------------
frac_canopy <- function(d, q10) {
  d$HRc <- d$HR * q10^(d$dT / 10)
  fp <- aggregate(cbind(HRc, TR) ~ plot + canopy + Year, data = d, FUN = sum)
  fp$f <- fp$HRc / fp$TR
  aggregate(f ~ canopy + Year, data = fp, FUN = mean)
}

# sensitivity to Q10_short
sens <- do.call(rbind, lapply(q10_short_range, function(q) {
  x <- frac_canopy(pr, q); x$Q10_short <- q; x
}))
sens_wide <- reshape(sens, idvar = c("canopy", "Year"), timevar = "Q10_short", direction = "wide")
names(sens_wide) <- sub("^f\\.", "f_Q10_", names(sens_wide))
cat("\n=== f_canopy for different short-term Q10 (Q10 = 1: no correction) ===\n")
print(sens_wide, row.names = FALSE, digits = 3)

# main result with bootstrap CI (resampling campaign dates within canopy)
f_main <- frac_canopy(pr, q10_short)
boot <- replicate(n_boot, {
  b <- do.call(rbind, lapply(split(pr, pr$canopy), function(d) {
    pick <- sample(unique(d$date), replace = TRUE)
    do.call(rbind, lapply(pick, function(x) d[d$date == x, ]))
  }))
  r <- frac_canopy(b, q10_short)
  r$f[match(paste(f_main$canopy, f_main$Year), paste(r$canopy, r$Year))]
})
f_main$lower95 <- apply(boot, 1, q95)[1, ]
f_main$upper95 <- apply(boot, 1, q95)[2, ]
names(f_main)[names(f_main) == "f"] <- "f_canopy"
cat("\n=== f_canopy at Q10_short =", q10_short, "(bootstrap 95% CI) ===\n")
print(f_main, row.names = FALSE, digits = 3)

write.csv(sens_wide, file.path(output_dir, "f_canopy_Q10_sensitivity.csv"), row.names = FALSE)
write.csv(f_main,    file.path(output_dir, "f_canopy_temperature_corrected.csv"), row.names = FALSE)

pr$HR_corr <- pr$HR * q10_short^(pr$dT / 10)
write.csv(pr, file.path(output_dir, "HR_TR_paired_corrected.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 7.4. ESTIMATING PARAMETERS ALPHA AND BETA
# -----------------------------------------------------------------------------
exclude_calibration <- c("CC_2022")   # canopy_year combinations not used (TR collar ~ HR)

if (!exists("pr") || !"HR_corr" %in% names(pr)) {
  stop("Object 'pr' with column HR_corr not found: run sections 6.1-6.3 first ",
       "(or source() the whole file).")
}

cal <- pr[pr$TR > 0 & pr$HR_corr > 0, ]
cal$Tc <- (cal$T_TR - 10) / 10
cal <- cal[!paste(cal$canopy, cal$Year, sep = "_") %in% exclude_calibration, ]

fit_f <- function(d, by_year) {
  yrs <- sort(unique(d$Year))
  na  <- if (by_year) length(yrs) else 1
  ia  <- if (by_year) match(d$Year, yrs) else rep(1, nrow(d))
  sse <- function(p) sum((log(d$HR_corr) - log(d$TR) -
                            plogis(p[ia] + p[na + 1] * d$Tc, log.p = TRUE))^2)
  o <- optim(c(rep(0.5, na), 0.5), sse, method = "BFGS")
  names(o$par) <- c(if (by_year) paste0("alpha_", yrs) else "alpha", "beta")
  list(par = o$par, sigma = sqrt(o$value / (nrow(d) - length(o$par))),
       n = nrow(d), converged = o$convergence == 0)
}

f_models <- list()
coef_tab <- data.frame()

for (cn in sort(unique(cal$canopy))) {
  d <- cal[cal$canopy == cn, ]
  for (version in c("pooled", "year")) {
    by_year <- version == "year"
    fit <- fit_f(d, by_year)
    
    draws <- t(replicate(n_boot, {
      pick <- sample(unique(d$date), replace = TRUE)
      db   <- do.call(rbind, lapply(pick, function(x) d[d$date == x, ]))
      if (by_year && length(unique(db$Year)) < length(unique(d$Year))) {
        rep(NA, length(fit$par))
      } else {
        fb <- try(fit_f(db, by_year), silent = TRUE)
        if (inherits(fb, "try-error") || !fb$converged) rep(NA, length(fit$par)) else fb$par
      }
    }))
    draws <- draws[complete.cases(draws), , drop = FALSE]
    colnames(draws) <- names(fit$par)
    
    f_models[[paste(cn, version, sep = "_")]] <- list(par = fit$par, draws = draws,
                                                      sigma = fit$sigma, canopy = cn,
                                                      version = version)
    ci <- apply(draws, 2, q95)
    coef_tab <- rbind(coef_tab, data.frame(
      canopy = cn, version = version, parameter = names(fit$par),
      estimate = fit$par, lower95 = ci[1, ], upper95 = ci[2, ],
      n = fit$n, sigma_log = fit$sigma, row.names = NULL))
  }
}

# f at 10 °C and Q10 ratio, for reading the table
coef_tab$f_at_10 <- ifelse(grepl("^alpha", coef_tab$parameter), plogis(coef_tab$estimate), NA)
coef_tab$Q10h_over_Q10a <- ifelse(coef_tab$parameter == "beta", exp(coef_tab$estimate), NA)

cat("\n=== Calibrated f(T) = plogis(alpha + beta * (T5 - 10) / 10) ===\n")
cat("    excluded from calibration:", paste(exclude_calibration, collapse = ", "), "\n")
print(coef_tab, row.names = FALSE, digits = 3)
write.csv(coef_tab, file.path(output_dir, "f_function_coefficients.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.5. ESTIMATING HETEROTROPHIC RESPIRATION FROM THE ORIGINAL FILE
# -----------------------------------------------------------------------------
# ---- Data files are here ----
GHG <- read.csv('1.GHG flux.csv')
str(GHG)
## We are interested in 2022 and 2023
GHG <- subset(GHG, Year > 2021 & Year < 2024)
summary(GHG)

# Parameters alpha and beta are added to GHG data - the pooled estimates
GHG$alpha <- numeric(length = nrow(GHG))
GHG$beta <- numeric(length = nrow(GHG))

for(i in 1:nrow(GHG)){
  if(GHG$Canopy_treatment[i] == "Non_harvest"){
    GHG$alpha[i] <- coef_tab[coef_tab$canopy == "CO" & coef_tab$version == 'pooled' & coef_tab$parameter == 'alpha', 'estimate']
    GHG$beta[i] <- coef_tab[coef_tab$canopy == "CO" & coef_tab$version == 'pooled' & coef_tab$parameter == 'beta', 'estimate']
  } else {
    if(GHG$Canopy_treatment[i] == 'CCF'){
      GHG$alpha[i] <- coef_tab[coef_tab$canopy == "CCF" & coef_tab$version == 'pooled' & coef_tab$parameter == 'alpha', 'estimate']
      GHG$beta[i] <- coef_tab[coef_tab$canopy == "CCF" & coef_tab$version == 'pooled' & coef_tab$parameter == 'beta', 'estimate']
    } else {
      GHG$alpha[i] <- coef_tab[coef_tab$canopy == "CC" & coef_tab$version == 'pooled' & coef_tab$parameter == 'alpha', 'estimate']
      GHG$beta[i] <- coef_tab[coef_tab$canopy == "CC" & coef_tab$version == 'pooled' & coef_tab$parameter == 'beta', 'estimate']
    }
  }
}

# Heterotrophic respiration is calculated - trenched plots in 2022 are not included in calculation because the flux represent heterotrophic respiration
GHG$HR <- GHG$CO2flux_mg_m2_h*plogis(GHG$alpha + GHG$beta*(GHG$T5 - 10)/10)

# plot of temperature dependence of heterotrophic respiration
ggplot(GHG, aes(T5, HR)) + geom_point(cex = 4, pch = 21, aes(fill = factor(Year))) + theme_test() +
  facet_wrap(Canopy_treatment ~ Disturbance_treatment, scales = 'free') +
  labs(x = 'Temperature (˚C)', y = expression(paste(CO[2]~flux~(mg~C-CO[2]~m^{-2}~h^{-1}))), fill = 'Year')

# =============================================================================
# ADDED: CONVERSION, EXPONENTIAL AND LINEAR FIT, R10 / R15 / R20 PER GROUP
# this part reuses your script
# =============================================================================
# convert from HR to  mg CO2-C m-2 h-1 (x 12/44)
GHG$HR <- GHG$HR * 12 / 44

# For each Year x Canopy x Disturbance group (4 plots pooled, as above):
#   exponential: log(flux) = log(A) + a*T   ->  flux = exp(log(A) + a*T) = A*exp(a*T)
#   linear     : flux = b0 + b1*T
# Both fits use the same rows (flux > 0 and T5 present).
# R2 of both on the original flux scale; the higher R2 is kept.
# If the kept fit is not significant (p of T5 >= p_limit), the mean flux is R10.

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
             !is.na(GHG$HR) & !is.na(GHG$T5) & GHG$HR > 0, ]
  
  if (nrow(d) < 5) {
    next
  }
  
  m_exp <- lm(log(HR) ~ T5, data = d)
  m_lin <- lm(HR ~ T5, data = d)
  
  # R2 on the original flux scale
  ss_tot <- sum((d$HR - mean(d$HR))^2)
  r2_exp <- 1 - sum((d$HR - exp(fitted(m_exp)))^2) / ss_tot
  r2_lin <- 1 - sum((d$HR - fitted(m_lin))^2) / ss_tot
  
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
    R <- c(mean(d$HR), NA, NA)
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
  
  points_y <- GHG[GHG$Year == y & !is.na(GHG$HR) & !is.na(GHG$T5) &
                    GHG$HR > 0, ]
  curves_y <- curves[curves$Year == y, ]
  res_y    <- CO2_res[CO2_res$Year == y & !is.na(CO2_res$chosen), ]
  
  points_y$panel <- paste(points_y$Canopy_treatment, points_y$Disturbance_treatment)
  curves_y$panel <- paste(curves_y$Canopy_treatment, curves_y$Disturbance_treatment)
  res_y$panel    <- paste(res_y$Canopy_treatment, res_y$Disturbance_treatment)
  
  panel_titles <- paste0(res_y$panel, ": ", res_y$chosen,
                         "\nR2 exp ", res_y$R2_exp, " | R2 lin ", res_y$R2_lin)
  names(panel_titles) <- res_y$panel
  
  p_fit <- ggplot() +
    geom_point(data = points_y, aes(x = T5, y = HR, colour = as.character(plot)),
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
