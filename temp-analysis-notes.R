
mod1 <- lm(CTmax.C ~ Species + Sex, data = rawTempData)
mod1b <- lm(CTmax.C ~ Species, data = rawTempData)

mod2 <- lm(CTmin.C~ Species + Sex, data = rawTempData)
mod2b <- lm(CTmin.C~Species, data = rawTempData)

mod3 <- lm(T.tot~ Species + Sex, data = rawTempData)
mod3b <- lm(T.tot~Species, data = rawTempData)

library(lmtest)
library(sandwich)

coeftest(mod1, vcov. = vcovHC(mod1, type= "HC3"))
confint(mod1)
coeftest(mod2, vcov. = vcovHC(mod2, type="HC3"))
confint(mod2)
coeftest(mod3, vcov. = vcovHC(mod2, type="HC3"))
confint(mod3)

coeftest(mod1b, vcov. = vcovHC(mod1b, type= "HC3"))
confint(mod1b)


simulationOutput <- simulateResiduals(mod1)

plot(simulationOutput)

plotResiduals(simulationOutput)
plotResiduals(simulationOutput, form = rawTempData$Sex)
plotResiduals(simulationOutput, form = rawTempData$Species)

simulationOutput2 <- simulateResiduals(mod2)

plot(simulationOutput2)

plotResiduals(simulationOutput2)
plotResiduals(simulationOutput2, form = rawTempData$Sex)
plotResiduals(simulationOutput2, form = rawTempData$Species)

simulationOutput3 <- simulateResiduals(mod3)

plot(simulationOutput3)

plotResiduals(simulationOutput3)
plotResiduals(simulationOutput3, form = rawTempData$Sex)
plotResiduals(simulationOutput3, form = rawTempData$Species)





m1r <- glm(CTmax.C ~ Sex + Species, family = gaussian, data=rawTempData)
m1br <- glm(CTmax.C ~ Species, family = gaussian, data = rawTempData)

model.sel(m1r, m1br)

plot(fitted(m1r), resid(m1r),
     xlab = "fitted values", 
     ylab = "residuals")
abline(h=0, lty=2)


plot(fitted(m1r), abs(resid(m1r)),
     xlab = "fitted values",
     ylab = "|residuals|")
lines(lowess(fitted(m1r), abs(resid(m1r))), col = "red")

dat <- rawTempData

df <- dat %>%
  group_by(Species, Sex) %>%
  summarise(
    mean_tt = mean(CTmax.C),
    var_tt = var(CTmax.C),
    .groups = "drop")

plot(df$mean_tt, df$var_tt) #slight negative trend, but highest variance for middle means
mod.var <- lm(var_tt~mean_tt, data = df)
abline(abline(a=coef(mod.var)[1], b=coef(mod.var)[2]))

df2 <- dat %>%
  group_by(Species, Sex) %>%
  summarise(
    mean_tt = mean(CTmin.C),
    var_tt = var(CTmin.C),
    .groups = "drop")

plot(df2$mean_tt, df2$var_tt) #very negative trend, lower variance for higher means
mod.var2 <- lm(var_tt~mean_tt, data = df2)
abline(abline(a=coef(mod.var2)[1], b=coef(mod.var2)[2]))

df3 <- dat %>%
  group_by(Species, Sex) %>%
  summarise(
    mean_tt = mean(T.tot),
    var_tt = var(T.tot),
    .groups = "drop")

plot(df3$mean_tt, df3$var_tt) #positive trend, with higher variance for higher means
mod.var3 <- lm(var_tt~mean_tt, data = df3)
abline(abline(a=coef(mod.var3)[1], b=coef(mod.var3)[2]))

library(e1071)

skewness(rawTempData$CTmax.C) #-0.51, fairly symmetric, slightly left skewed
plot(density(rawTempData$CTmax.C, na.rm = TRUE))
skewness(rawTempData$CTmin.C) #-0.60, fairly symmetric, slightly left skewed
plot(density(rawTempData$CTmin.C, na.rm = TRUE))
skewness(rawTempData$T.tot)   #0.094, fairly symmetric, ever so slightly right skewed
plot(density(rawTempData$T.tot, na.rm = TRUE))



m2r <- glm(CTmin.C ~ Sex + Species, family = gaussian, data=rawTempData)

plot(fitted(m2r), resid(m2r),
     xlab = "fitted values", 
     ylab = "residuals")
abline(h=0, lty=2)


plot(fitted(m2r), abs(resid(m2r)),
     xlab = "fitted values",
     ylab = "|residuals|")
lines(lowess(fitted(m2r), abs(resid(m2r))), col = "red")


m3r <- glm(T.tot ~ Sex + Species, family = gaussian, data=rawTempData)

plot(fitted(m3r), resid(m3r),
     xlab = "fitted values", 
     ylab = "residuals")
abline(h=0, lty=2)


plot(fitted(m3r), abs(resid(m3r)),
     xlab = "fitted values",
     ylab = "|residuals|")
lines(lowess(fitted(m3r), abs(resid(m3r))), col = "red")

dat %>%
  group_by(Species, Sex) %>%
  summarise(
    mean_tt = mean(CTmin.C),
    var_tt = var(CTmin.C),
    .groups = "drop")



library(dplyr)

# Function to simulate power and find minimum detectable sex difference
# Inputs:
#   nsim: number of simulations per effect size
#   species: vector of species names
#   nM, nF: vectors of males/females per species
#   sigma: residual SD from real LM
#   alpha: significance threshold
#   power_target: desired power (e.g., 0.8)
#   delta_range: vector of sex differences to test
#   trait_name: name of the trait (for labeling)
simulate_MDD <- function(nsim=1000, species, nM, nF, sigma, alpha=0.05,
                         power_target=0.8, delta_range=seq(0, 3, by=0.1),
                         trait_name="Trait") {
  
  power_results <- data.frame(delta = delta_range, power = NA)
  
  for(i in seq_along(delta_range)){
    sig_count <- 0
    for(j in 1:nsim){
      temp <- c()
      sex <- c()
      sp <- c()
      
      for(k in seq_along(species)){
        # simulate male and female data around delta/2 difference
        temp <- c(temp,
                  rnorm(nM[k], mean = delta_range[i]/2, sd = sigma),
                  rnorm(nF[k], mean = -delta_range[i]/2, sd = sigma))
        sex <- c(sex, rep(c("M","F"), times=c(nM[k], nF[k])))
        sp <- c(sp, rep(species[k], times=nM[k]+nF[k]))
      }
      
      df <- data.frame(temp=temp, sex=sex, species=sp)
      
      # Fit LM
      m <- lm(temp ~ sex + species, data=df)
      
      # extract p-value for sex
      pval <- summary(m)$coefficients["sexM","Pr(>|t|)"]
      if(!is.na(pval) & pval < alpha) sig_count <- sig_count + 1
    }
    
    power_results$power[i] <- sig_count / nsim
  }
  
  # Find minimum detectable delta (closest to target power)
  closest <- power_results$delta[which.min(abs(power_results$power - power_target))]
  
  cat("\n---", trait_name, "---\n")
  cat("Minimum detectable sex difference for", power_target*100, "% power: ", closest, "\n\n")
  
  return(list(power_curve = power_results, MDD = closest))
}

# ==========================
# Example usage for your data
species <- c("sailfin","mosquitofish","swordtail")
nM <- c(16, 9, 15)
nF <- c(16, 8, 15)

# Replace these sigma values with the residual SD from your real LMs
sigma_CTmax <- 1.025751   # summary(mod1)$sigma
sigma_CTmin <- 1.414415   # summary(mod2)$sigma
sigma_Ttot  <- 2.030924   # summary(mod3)$sigma

# Run simulations
MDD_CTmax <- simulate_MDD(species=species, nM=nM, nF=nF, sigma=sigma_CTmax, trait_name="CTmax")
MDD_CTmin <- simulate_MDD(species=species, nM=nM, nF=nF, sigma=sigma_CTmin, trait_name="CTmin")
MDD_Ttot  <- simulate_MDD(species=species, nM=nM, nF=nF, sigma=sigma_Ttot,  trait_name="Ttot")

CTmax_df <- MDD_CTmax$power_curve %>% mutate(Trait="CTmax", MDD=MDD_CTmax$MDD)
CTmin_df <- MDD_CTmin$power_curve %>% mutate(Trait="CTmin", MDD=MDD_CTmin$MDD)
Ttot_df  <- MDD_Ttot$power_curve  %>% mutate(Trait="Ttot",  MDD=MDD_Ttot$MDD)
power_all <- bind_rows(CTmax_df, CTmin_df, Ttot_df)

# Plot
ggplot(power_all, aes(x=delta, y=power, color=Trait)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=0.8, linetype="dashed", color="black") +
  geom_vline(aes(xintercept=MDD, color=Trait), linetype="dotted") +
  labs(title="Power curves and minimum detectable sex differences",
       x="Sex difference (°C)",
       y="Power",
       color="Trait") +
  theme_minimal() +
  theme(text = element_text(size=12))


males <- rawTempData[rawTempData$Sex == "male",]
females <- rawTempData[rawTempData$Sex == "female",]

(mean(males$CTmax.C))- (mean(females$CTmax.C))
(mean(males$CTmin.C))- (mean(females$CTmin.C))
(mean(males$T.tot))- (mean(females$T.tot))

# OK HERE'S WHAT WE'RE GONNA DO.
## model selection: run lms for interaction, size, and sex. show that model with just species is the most appropriate, but will report the effect + CI for sex
### run DHARMa on lm models to check residuals. they will fail. so then calculate robust standard errors for each selected model (i believe with and without sex gives the same p values for species, so i will include sex to again show no difference)
#### use code above to detect power needed for sex differences. report this and the actual mean differences we had to show that while we might not have found a significant effect of sex in our models, there is a possibility that the small sex differences that do exist could be biologically relevant give a bigger sample size
##### from gpt: The observed mean difference between males and females was 0.12 °C for CTmax. Simulation-based power analysis indicated that, given our sample sizes and residual variation, we could reliably detect sex differences ≥0.7 °C with 80% power. Thus, while we found no statistically significant sex effect, smaller biologically meaningful differences could not be ruled out due to limited sample size.