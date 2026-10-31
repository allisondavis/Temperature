# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────
###### MANUSCRIPT TABLES AND FIGURES ######   NOTE: must run all Rmd files prior to this script
# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────
# NOTE: specific Rmd files needed to run each table/figure are noted in individual headers


# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────
### TABLES ###    NOTE: all tables created using gt() then exported to Word for final editing
# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 1 
# ───────────────────────────────────────────────────────────────────────────────────────────────

#install.packages("tibble")     #v3.3.1
library(tibble)

#Create data set
lit_table <- tribble(
  ~Species, ~Trait, ~Finding, ~Source,
  
  # ── Poecilia latipinna ─────────────────────────────
  "Poecilia latipinna", "Morphology", "Males are smaller and more colorful than females on average", "Scharnweber et al. 2011",
  "Poecilia latipinna", "Feeding and Growth", "Males feed for less time than females due to increased mating behavior", "Scharnweber et al. 2011",
  "Poecilia latipinna", "Feeding and Growth", "Males have shorter intestinal tracts than females, although gut composition is similar", "Scharnweber et al. 2011",
  "Poecilia latipinna", "Feeding and Growth", "Females grow faster than males after reaching sexual maturity", "Snelson 1982",
  "Poecilia latipinna", "Feeding and Growth", "Both sexes mature at larger sizes and later ages in lower salinity and cooler environments", "Trexler et al. 1990",
  "Poecilia latipinna", "Feeding and Growth", "Late-gestation females exhibit higher mass-adjusted routine metabolic rates compared to males", "Timmerman & Chapman 2003",
  "Poecilia latipinna", "Behavior", "Atrazine exposure alters female mate preference and reduces male aggression; boldness and anxiety are affected in both sexes", "MacLaren 2023",
  "Poecilia latipinna", "Behavior", "Females exhibit lower neophobia than males", "Fuss & Witte 2019",
  "Poecilia latipinna", "Behavior", "Males outperform females in color discrimination and reversal learning tasks", "Fuss & Witte 2019",
  "Poecilia latipinna", "Population Dynamics", "Sex ratios are approximately 1:1 among juveniles but female-biased among adults in natural Texas populations", "Baird 1973; Hubbs 1964; Simpson & Gunter 1956; Snelson & Wetherington, 1980",
  "Poecilia latipinna", "Population Dynamics", "Male age and size at maturity vary among populations, while female traits remain consistent", "Trexler et al. 1990",
  "Poecilia latipinna", "Phenotypic Plasticity", "Newborn males are less plastic than females in response to temperature, salinity, and food availability", "Trexler et al. 1990",
  
  # ── Gambusia affinis ─────────────────────────────
  "Gambusia affinis", "Morphology", "Females have higher dry mass, energy reserves, and condition factors but lower burst and critical swimming speeds than males", "Li et al. 2017",
  "Gambusia affinis", "Morphology", "Males exhibit longer peduncle lengths and shorter body heights than females", "Blanchard et al. 2025",
  "Gambusia affinis", "Morphology", "Females have larger caudal-fin aspect ratios and total lengths, enhancing dispersal capacity compared to males", "Blanchard et al. 2025",
  "Gambusia affinis", "Morphology", "Larger females show higher infection rates by the monogenean parasite Salsuginus seculus compared to males", "Renner & Duggan 2024",
  "Gambusia affinis", "Behavior", "Females show increased shoaling and risk sensitivity to mitigate male harassment", "Cummings 2018",
  "Gambusia affinis", "Behavior", "Females exhibit a stronger 'fast-exploratory' cognitive style and higher spatial-temporal learning tasks compared to males", "Wallace et al. 2020",
  "Gambusia affinis", "Behavior", "Females forage more than males, but mixed-sex group foraging is reduced in both sexes", "Arrington et al. 2009",
  "Gambusia affinis", "Behavior", "Males engage in significantly more aggressive attacks than females", "Walsh et al. 2025",
  "Gambusia affinis", "Population Dynamics", "Females are slightly larger than males, increasing vulnerability to predation", "Britton & Moser 1982",
  "Gambusia affinis", "Population Dynamics", "Natural populations are female-biased in adult sex ratios", "Krumholz 1948",
  "Gambusia affinis", "Population Dynamics", "Males experience higher mortality under extreme hypoxia compared to females", "Cech Jr. et al. 1985",
  "Gambusia affinis", "Population Dynamics", "Female fitness declines with increasing female density, independent of male harassment", "Smith 2007",
  "Gambusia affinis", "Population Dynamics", "Environmental contamination reduces male sperm counts and female egg size and number", "Xie et al. 2010",
  "Gambusia affinis", "Phenotypic Plasticity", "Females tolerate colder temperatures better than males", "Wood et al. 2020",
  "Gambusia affinis", "Phenotypic Plasticity", "High predator density selects for larger caudal regions and enhanced burst-swimming performance in both sexes", "Langerhans 2009",
  
  # ── Xiphophorus nigrensis ─────────────────────────────
  "Xiphophorus nigrensis", "Feeding and Growth", "Large males exhibit higher oxygen consumption when females are present, indicating increased energetic costs", "Cummings & Gelineau-Kattner 2009",
  "Xiphophorus nigrensis", "Feeding and Growth", "Females mature earlier than intermediate and large males, while body size continues to increase post-maturity", "Ryan 1989")



# Process and Order Data 
lit_table <- lit_table %>%
  mutate(Trait = factor(Trait, levels = c(
    "Morphology",
    "Feeding and Growth",
    "Behavior",
    "Population Dynamics",
    "Phenotypic Plasticity"))) %>%
  arrange(Species, Trait)

#Final table
gt_table <- lit_table %>%
  gt(groupname_col = "Species") %>%
  tab_header(
    title = md("**Table 1. Sex-specific traits across poeciliid fishes**")) %>%
  cols_label(
    Trait   = "Trait category",
    Finding = "Key finding",
    Source  = "Source") %>%
  fmt_markdown(columns = Source) %>%
  cols_width(
    Trait   ~ px(120),
    Finding ~ px(440), # Width optimized to allow wrapping
    Source  ~ px(150)) %>%
  tab_options(
    table.font.size            = px(10),
    data_row.padding           = px(3),
    row_group.font.weight      = "bold",
    row_group.background.color = "#F5F5F5",
    heading.align              = "left")

gtsave(gt_table, "Table1_sex-specific-traits.rtf")

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 2          NOTE: run trial_rate_calc.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

tt  <- tt  %>% mutate(Test = "CTmax")
tt2 <- tt2 %>% mutate(Test = "CTmin -0.3")
tt3 <- tt3 %>% mutate(Test = "CTmin -0.1")


tt.all <- bind_rows(tt, tt2, tt3)
tt.all <- tt.all %>%
  mutate(
    result = case_when(
      p.value < 0.001 ~ paste0("*p* < 0.001<br>(t = ", statistic,")"),
      p.value < 0.05  ~ paste0(
        "*p* = ", sprintf("%.3f", p.value),
        "<br>(t = ", statistic, ")"),
      TRUE ~ paste0(
        "*p* = ", sprintf("%.3f", p.value),
        "<br>(t = ", statistic, ")")))

tt.wide <- tt.all %>%
  select(Species, Test, result) %>%
  pivot_wider(
    names_from  = Test,
    values_from = result)

spp.trial.rates <- tt.wide %>%
  gt() %>%
  cols_label(
    Species   = "Species") %>%
  fmt_markdown(everything()) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(`CTmax`, `CTmin -0.3`),
      rows = Species == "Swordtail")) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(Species, CTmax, `CTmin -0.3`))) %>%
  cols_align(align = "center", everything()) %>%
  opt_horizontal_padding(scale = 2)
spp.trial.rates

gtsave(spp.trial.rates, "Table2_species-trial-rates.rtf")


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 3        NOTE: run trial_rate_calc.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

rate.sum.table <- trial.rates2 %>%
  gt() %>%
  cols_label(Trial ="Trial type",
             Species = "Species",
             trial.rate = "Rate change (\u00B0C/min)",
             stand.dev = "Standard deviation",
             min = "Minimum rate (\u00B0C)",
             Median = "Median rate (\u00B0C)",
             max = "Maximum rate (\u00B0C)") %>%
  fmt_number(columns = c("trial.rate", "stand.dev", "min", "max", "Median"), decimals = 3) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c("Trial", "Species", "trial.rate", "stand.dev", "min", "Median"))) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "darkgray",
      weight = px(2.5)),
    locations = cells_body(rows = 3)) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())) %>%
  opt_horizontal_padding(scale=2)
rate.sum.table

gtsave(rate.sum.table, "Table3_rate-summaries.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 4        NOTE: run Austin-Daily-Temps.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

annual_avg_rates <- slopes.all %>%
  group_by(Year) %>%
  summarise(
    avg_positive_rate = mean(slope_deg_per_minute[slope_deg_per_minute > 0], na.rm = TRUE),
    avg_negative_rate = mean(slope_deg_per_minute[slope_deg_per_minute < 0], na.rm = TRUE)) %>%
  gt() %>%
  cols_label(
    Year = "Year",
    avg_positive_rate = enc2utf8("Average increase (\u00B0C/min)"),
    avg_negative_rate = enc2utf8("Average decrease (\u00B0C/min)")) %>%
  fmt_number(columns = c(avg_positive_rate, avg_negative_rate), decimals = 6) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(Year, avg_positive_rate))) %>%
  cols_align(align="center", everything()) %>%
  opt_horizontal_padding(scale = 2)


gtsave(annual_avg_rates, "Table4_austin-rates.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 5        NOTE: run Final_temperature-analysis.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

cr2.table <- function(model, trait, label) {
  V <- vcovCR(model, type = "CR2")
  
  as.data.frame(
    conf_int(model, vcov = V, level = 0.95, test = "Satterthwaite", p_values = TRUE)) %>%
    filter(Coef %in% c("Sexmale", "Size.mm")) %>%
    transmute(
      Trait = trait,
      Model = label,
      Predictor = Coef,
      CR2_Estimate = beta,
      CR2_SE = SE,
      CR2_low = CI_L,
      CR2_high = CI_U,
      CR2_p = p_val)
}

cr2.results <- bind_rows(
  cr2.table(mod2b, "CTmin", "B"),
  cr2.table(mod2c, "CTmin", "C"),
  cr2.table(mod2d, "CTmin", "D"),
  cr2.table(mod1b, "CTmax", "B"),
  cr2.table(mod1c, "CTmax", "C"),
  cr2.table(mod1d, "CTmax", "D"))

model.results <- bind_rows(
  ctmin.results %>% mutate(Trait = "CTmin"),
  ctmax.results %>% mutate(Trait = "CTmax")) %>%
  transmute(
    Trait,
    Model,
    Predictor,
    MB_Estimate = Estimate,
    MB_SE = `Std. Error`,
    MB_low = CI_low,
    MB_high = CI_high,
    MB_df = df,
    MB_p = `Pr(>|t|)`)

combined.results <- model.results %>%
  left_join(cr2.results, by = c("Trait", "Model", "Predictor")) %>%
  mutate(
    Trait = factor(Trait, levels = c("CTmin", "CTmax")),
    Model = factor(Model, levels = c("B", "C", "D")), 
    Effect = case_when(
      Predictor == "Sexmale" ~ "Sex",
      Predictor == "Size.mm" ~ "Size", 
      TRUE ~ Predictor)) %>%
  arrange(Trait, Model, match(Predictor, c("Sexmale", "Size.mm")))

combined.results <- combined.results %>%
  mutate(
    p.disp.mb = case_when(
      MB_p < 0.001 ~ "*p* < 0.001",
      MB_p < 0.05  ~ paste0(sprintf("%.3f", MB_p)),
      TRUE            ~ sprintf("%.3f", MB_p))) %>%
  mutate(
    p.disp.cr2 = case_when(
      CR2_p < 0.001 ~ "*p* < 0.001",
      CR2_p < 0.05  ~ paste0(sprintf("%.3f", CR2_p)),
      TRUE            ~ sprintf("%.3f", CR2_p)))

table.results <- combined.results %>%
  mutate(
    Model = dplyr::recode(
      as.character(Model),
      "B" = "Species + Size",
      "C" = "Species + Sex",
      "D" = "Species + Sex \n+Size")) %>%
  transmute(
    Trait,
    Model,
    Effect,
    mb_est = sprintf("%.3f (%.3f)", MB_Estimate, MB_SE),
    mb_ci = sprintf("[%.3f, %.3f]", MB_low, MB_high),
    mb_df = sprintf("%.2f", MB_df),
    mb_p = p.disp.mb,
    
    cr2_est = sprintf("%.3f (%.3f)", CR2_Estimate, CR2_SE),
    cr2_ci = sprintf("[%.3f, %.3f]", CR2_low, CR2_high),
    cr2_p = p.disp.cr2)

coeff.table <- table.results %>%
  gt(rowname_col = "Model", groupname_col = "Trait", row_group_as_column = TRUE) %>%
  tab_stubhead(label = "Model") %>%
  tab_spanner(label = "Model-based inference", columns = c(Effect, mb_est, mb_ci, mb_df, mb_p)) %>%
  tab_spanner(label = "CR2 analysis", columns = c(cr2_est, cr2_ci, cr2_p)) %>%
  cols_label(
    Effect = "Effect",
    mb_est = "Estimate (SE)",
    mb_ci = "95% CI", 
    mb_df = "df",
    mb_p = md("*p*-value"), 
    cr2_est = "Estimate (SE)",
    cr2_ci = "95% CI",
    cr2_p = md("*p*-value")) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(Effect,mb_est, mb_ci,mb_df, cr2_est, cr2_ci))) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(2)),
    locations = cells_body(columns = c(mb_p))) %>%
  cols_align("center", everything()) %>%
  opt_horizontal_padding(scale = 2)

coeff.table

gtsave(coeff.table, "Table5_coefficients.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 6        NOTE: run Final_temperature-analysis.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

sum.stats.tab <- summ.stats.all %>%
  arrange(Species, Trial) %>%
  select(Species, Trial, mean_sd, min, med, max, male, female) %>%
  gt(rowname_col = "Trial", groupname_col = "Species") %>%
  cols_label(
    Species = "Species",
    Trial   = "Trial",
    mean_sd = "Mean (\u00B0C) (SD)",
    min     = "Minimum (\u00B0C)",
    med     = "Median (\u00B0C)",
    max     = "Maximum (\u00B0C)",
    male    = "Males",
    female  = "Females") %>%
  fmt_number(
    columns = c(min, med, max),
    decimals = 2) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(mean_sd, min, med, max, male))) %>%
  cols_align("center", everything()) %>%
  opt_horizontal_padding(scale = 2)


gtsave(sum.stats.tab, "TB6_sum-stats.rtf")




# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 7
# ───────────────────────────────────────────────────────────────────────────────────────────────

df <- data.frame(
  type = c("CTmax", "CTmin"),
  Actual   = c(0.2, 0.1),
  Simulated   = c(0.6, 0.4))

power.an <- df %>%
  gt() %>%
  cols_label(type ="Temperature Category",
             Actual = "Actual difference (\u00B0C)",
             Simulated = "Simulated difference needed (\u00B0C)") %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(type, Actual))) %>%
  cols_align(align = "center", everything())  %>%
  opt_horizontal_padding(scale=2)
power.an

gtsave(power.an, "Table7_power-analysis.rtf")


# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────
### FIGURES ###   NOTE: Figure 1 created in Inkscape
# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────

### Used to imbed fonts
#install.packages("ragg")
library(ragg)
#install.packages("systemfonts")
library(systemfonts)

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 2       NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────
size.pvals <- data.frame(
  Species = c("Swordtail", "Swordtail", "Swordtail", "Sailfin", "Mosquitofish"),
  group1 = c("Female", "Male (large)", "Female", "Female", "Female"),
  group2 = c("Male (large)", "Male (small)", "Male (small)", "Male", "Male"),
  label  = c("p < 0.001*", "p < 0.001*", "p = 0.001*", "ns", "p < 0.001*"),
  y.position = c(45, 46.5, 48, 47, 47))

size.pvals$fontface <- ifelse(size.pvals$label == "ns", "plain", "bold")

size.plot <- ggplot(dat.plot, aes(x = Sex.plot, y = Size.mm)) +
  geom_rect(
    data = subset(datTot, Species == "Sailfin"),
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "grey95") +
  geom_boxplot(linewidth= 0.4) +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  facet_wrap(~ Species, scales = "free_x") +
  scale_y_continuous(breaks = scales::breaks_pretty(n=8)) +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% c("Male (small)", "Male (large)"),
             gsub(" ", "\n", x),
             x)}) +
  stat_pvalue_manual(
    size.pvals,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5,
    fontface = "fontface") +
  theme_classic() + 
  xlab("") +
  ylab("Length (mm)")

#Save as tiff
ggsave(
  filename = "Fig2_sizes.tiff",
  plot = size.plot,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)

#save as eps
ggsave(
  filename = "Fig2_sizes.eps",
  plot = size.plot,
  device = cairo_ps,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 3       NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────

#install.packages("scales")   #v1.4.0

### CTmax - EMM
df.pvals.max <- data.frame(
  group1 = c("Swordtail", "Sailfin", "Swordtail"),
  group2 = c("Sailfin", "Mosquitofish", "Mosquitofish"),
  label  = c("p < 0.001*", "p = 0.013*", "p < 0.001*"),
  y.position = c(39, 39, 39))

df.pvals.max$fontface <- ifelse(df.pvals.max$label == "ns", "plain", "bold")

PW.df.max <- as.data.frame(confint(PW.sp.max))

emmean.max <- ggplot(PW.df.max, aes(x=Species, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.05, 0.25))) +
  stat_pvalue_manual(
    df.pvals.max,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5,
    step.increase = 0.1,
    fontface = "fontface") +
  ylab("Estimated mean CTmax (\u00B0C) \n") +
  xlab("\n Species") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
emmean.max


### CTmin-EMM

df.pvals.min <- data.frame(
  group1 = c("Swordtail", "Sailfin", "Swordtail"),
  group2 = c("Sailfin", "Mosquitofish", "Mosquitofish"),
  label  = c("p = 0.004*", "ns", "p < 0.001*"),
  y.position = c(11.5, 11.5, 11.5))

df.pvals.min$fontface <- ifelse(df.pvals.min$label == "ns", "plain", "bold")

PW.df.min <- as.data.frame(confint(PW.sp.min))

emmean.min <- ggplot(PW.df.min, aes(x=Species, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = c(0.05, 0.25))) +
  ylab("Estimated mean CTmin (\u00B0C) \n") +
  stat_pvalue_manual(
    df.pvals.min,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5,
    step.increase = 0.1,
    fontface = "fontface") +
  xlab("\n Species") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
emmean.min

final.plot_em.spp<- ((emmean.max + emmean.min) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a"))
final.plot_em.spp

#Save as tiff
ggsave(
  filename = "Fig3_em-species.tiff",
  plot = final.plot_em.spp,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 4,          
  units = "in",
  dpi = 600)

#save as eps
ggsave(
  filename = "Fig3_em-species.eps",
  plot = final.plot_em.spp,
  device = cairo_ps,
  width = 6.85, 
  height = 4,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 4       NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────

#install.packages("lme4")

### CTmax- sex and size

# Extract observations used in the full model
df.max <- model.frame(mod1d)

# Generate predictions only within each observed species x sex size range
pred.max <- df.max %>%
  group_by(Species, Sex) %>%
  reframe(
    Size.mm = seq(
      min(Size.mm),
      max(Size.mm),
      length.out = 50))

# Construct fixed-effect design matrix
X <- model.matrix(
  delete.response(
    terms(lme4::nobars(formula(mod1d)))),
  data = pred.max)

# Model coefficients and covariance matrix
beta <- lme4::fixef(mod1d)
V <- as.matrix(vcov(mod1d))

X <- X[, names(beta), drop = FALSE]
V <- V[names(beta), names(beta), drop = FALSE]

# Predicted means and approximate 95% confidence intervals
pred.max$fit <- as.numeric(X %*% beta)

pred.max$SE <- sqrt(
  rowSums((X %*% V) * X))

pred.max$lower <- pred.max$fit - 1.96 * pred.max$SE
pred.max$upper <- pred.max$fit + 1.96 * pred.max$SE

# Plot
sex.size.plot.max <- ggplot(df.max, aes(x = Size.mm, y = CTmax.C, color = Sex)) +
  geom_point(aes(shape= Sex), alpha = 0.6, size = 2) +
  geom_ribbon(data = pred.max, aes(
      x = Size.mm,
      ymin = lower,
      ymax = upper,
      fill = Sex,
      group = Sex),
    inherit.aes = FALSE,
    alpha = 0.15) +
  geom_line(
    data = pred.max, aes(x = Size.mm, y = fit, color = Sex, group = Sex),
    inherit.aes = FALSE,
    linewidth = 1) +
  facet_wrap(~ Species, scales = "free_x") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
  labs(
    x = "Body size (mm)",
    y = "CTmax (°C)",
    color = "Sex",
    fill = "Sex") +
  scale_shape_manual(name = "Sex",
                     values = c("female" = 16, "male" = 17),
                     labels = c("female"="Female", "male"="Male")) +
  scale_color_manual(name="Sex",
                     values = c("female"="#E1AF00", "male"="#78B7C5"),
                     labels = c("female"="Female", "male"="Male")) +
  scale_fill_manual(name="Sex",
                     values = c("female"="#E1AF00", "male"="#78B7C5"),
                     labels = c("female"="Female", "male"="Male")) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))

sex.size.plot.max



### CTmin- sex and size
# Extract observations used in the full model
df.min <- model.frame(mod2d)

# Generate predictions only within each observed species x sex size range
pred.min <- df.min %>%
  group_by(Species, Sex) %>%
  reframe(
    Size.mm = seq(
      min(Size.mm),
      max(Size.mm),
      length.out = 50))

# Construct fixed-effect design matrix
X2 <- model.matrix(
  delete.response(
    terms(lme4::nobars(formula(mod2d)))),
  data = pred.min)

# Model coefficients and covariance matrix
beta2 <- lme4::fixef(mod2d)
V2 <- as.matrix(vcov(mod2d))

X2 <- X2[, names(beta2), drop = FALSE]
V2 <- V2[names(beta2), names(beta2), drop = FALSE]

# Predicted means and approximate 95% confidence intervals
pred.min$fit <- as.numeric(X2 %*% beta2)

pred.min$SE <- sqrt(
  rowSums((X2 %*% V2) * X2))

pred.min$lower <- pred.min$fit - 1.96 * pred.min$SE
pred.min$upper <- pred.min$fit + 1.96 * pred.min$SE

# Plot
sex.size.plot.min <- ggplot(df.min, aes(x = Size.mm, y = CTmin.C, color = Sex)) +
  geom_point(aes(shape= Sex), alpha = 0.6, size = 2) +
  geom_ribbon(data = pred.min, aes(
    x = Size.mm,
    ymin = lower,
    ymax = upper,
    fill = Sex,
    group = Sex),
    inherit.aes = FALSE,
    alpha = 0.15) +
  geom_line(
    data = pred.min, aes(x = Size.mm, y = fit, color = Sex, group = Sex),
    inherit.aes = FALSE,
    linewidth = 1) +
  facet_wrap(~ Species, scales = "free_x") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6)) +
  labs(
    x = "Body size (mm)",
    y = "CTmin (°C)",
    color = "Sex",
    fill = "Sex") +
  scale_shape_manual(name = "Sex",
                     values = c("female" = 16, "male" = 17),
                     labels = c("female"="Female", "male"="Male")) +
  scale_color_manual(name="Sex",
                     values = c("female"="#E1AF00", "male"="#78B7C5"),
                     labels = c("female"="Female", "male"="Male")) +
  scale_fill_manual(name="Sex",
                    values = c("female"="#E1AF00", "male"="#78B7C5"),
                    labels = c("female"="Female", "male"="Male")) +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))

sex.size.plot.min



### Combining

final.plot_sex.size<- ((sex.size.plot.max/sex.size.plot.min) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom"))
final.plot_sex.size

#Save as tiff
ggsave(
  filename = "Fig4_sex.size.tiff",
  plot = final.plot_sex.size,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)

#save as eps
ggsave(
  filename = "Fig4_sex.size.eps",
  plot = final.plot_sex.size,
  device = cairo_ps,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 5       NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────

dat.plot3 <- rawTempData %>%
  mutate(
    Sex.plot3 = case_when(
      Sex == "male" ~ "Male",
      Sex == "female" ~ "Female")) %>%
  mutate(T.tot = CTmax.C - CTmin.C)

range.plot <- ggplot(dat.plot3, aes(x = Sex.plot3, y = T.tot)) +
  geom_rect(
    data = subset(datTot, Species == "Sailfin"),
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "grey95") +
  geom_boxplot() +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  facet_wrap(~ Species, scales = "free_x") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Temperature range (\u00B0C) \n(CTmax-CTmin)")
range.plot

#Save as tiff
ggsave(
  filename = "Fig5_species-range.tiff",
  plot = range.plot,
  device = ragg::agg_tiff,
  width = 3.3, 
  height = 4,          
  units = "in",
  dpi = 600)


#Save as eps
ggsave(
  filename = "Fig5_species-range.eps",
  plot = range.plot,
  device = cairo_ps,
  width = 3.3, 
  height = 4,          
  units = "in",
  dpi = 600)

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 6       NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────

### CTmax-EMM
PW.df.maxSW <- as.data.frame(confint(PW.sw.max))

SW.em.max <- ggplot(PW.df.maxSW, aes(x=Sex.plot, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6)) +
  coord_cartesian(ylim = c(34, 37)) +
  ylab("Estimated mean CTmax (\u00B0C) \n") +
  xlab("\n Sex category") +
  ggtitle("CTmax") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
SW.em.max

### CTmin-EMM
PW.df.minSW <- as.data.frame(confint(PW.sw.min))

SW.em.min <- ggplot(PW.df.minSW, aes(x=Sex.plot, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(9, 11)) +
  ylab("Estimated mean CTmin (\u00B0C) \n") +
  xlab("\n Sex category") +
  ggtitle("CTmin") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
SW.em.min


### Combining

final.plot_em.sw <- ((SW.em.max + SW.em.min) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "a"))
final.plot_em.sw

#Save as tiff
ggsave(
  filename = "Fig6_em-sw.tiff",
  plot = final.plot_em.sw,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 4,          
  units = "in",
  dpi = 600)


#Save as eps
ggsave(
  filename = "Fig6_em-sw.eps",
  plot = final.plot_em.sw,
  device = cairo_ps,
  width = 6.85, 
  height = 4,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
### SUPPLEMENTARY TABLES/FIGURES ###
# ───────────────────────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Online Resource 1 (run Rerun-outlier-removal.Rmd before processing figure below)
# ───────────────────────────────────────────────────────────────────────────────────────────────

out.min <- ggplot(df_flagged, aes(x = CTmin.C, fill = is_outlier_min)) +
  geom_histogram(binwidth = 0.5, color = "white") +
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02"),
                    name = "Data Status",
                    labels = c("Normal", "Outlier"),
                    drop = FALSE) +
  labs(title = "CTmin",
       x = "Critical Thermal Minimum (°C)",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none")

outliers <- ((out.min/out.max) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect"))
outliers

#Save as tiff
ggsave(
  filename = "OnlineResource1_outliers.tiff",
  plot = outliers,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)


#Save as eps
ggsave(
  filename = "OnlineResource1.eps",
  plot = outliers,
  device = cairo_ps,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Online Resource 2 (run Austin-Daily-Temps.Rmd before processing table below)
# ───────────────────────────────────────────────────────────────────────────────────────────────

extreme.table <- extremes.wide %>%
  gt() %>%
  tab_spanner(
    label = "2020",
    columns = c(Year_2020_max_slope, Year_2020_min_slope)) %>%
  tab_spanner(
    label = "2021",
    columns = c(Year_2021_max_slope, Year_2021_min_slope)) %>%
  tab_spanner(
    label = "2022",
    columns = c(Year_2022_max_slope, Year_2022_min_slope)) %>%
  tab_spanner(
    label = "2023",
    columns = c(Year_2023_max_slope, Year_2023_min_slope)) %>%
  tab_spanner(
    label = "2024",
    columns = c(Year_2024_max_slope, Year_2024_min_slope)) %>%
  cols_label(
    Month = "Month") %>%
  cols_label_with(
    fn = function(x) {
      ifelse(grepl("max_slope", x), "Max \nrate",
             ifelse(grepl("min_slope", x), "Min \nrate", x))}) %>%
  fmt_number(
    columns = -Month,
    decimals = 3) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = Month)) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5), style = "dashed"),
    locations = cells_body(columns = contains("max_slope"))) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = contains("min_slope"))) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_column_labels(
      columns = contains("min_slope"))) %>%
  cols_align("center", everything()) %>%
  opt_horizontal_padding(scale=2)
extreme.table

gtsave(extreme.table, "OnlineResource2_extreme-rates-austin.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Online Resource 3 (run Final_temperature-analysis. Rmd before processing tables/figures below)
# ───────────────────────────────────────────────────────────────────────────────────────────────

# ----------------------------------------------------------------------------------
### CTmax & CTmin model selection ###
# ----------------------------------------------------------------------------------
max.sel <- model.sel(mod1, mod1b, mod1c, mod1d, mod1e, mod1f, mod1g)
min.sel <- model.sel(mod2, mod2b, mod2c, mod2d, mod2e, mod2f, mod2g)

prepare.sel <- function(sel, trait) {
  as.data.frame(sel) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(Trait = trait) %>%
    select(Trait, Model, df, AICc, delta, weight)}

sel.all <- bind_rows(
  prepare.sel(min.sel, "CTmin"),
  prepare.sel(max.sel, "CTmax")) %>%
  mutate(Trait = factor(Trait, levels = c("CTmin", "CTmax"))) %>%
  mutate(Model = dplyr::recode(
    as.character(Model),
    mod1 = "Species",
    mod1b = "Species + Size",
    mod1c = "Species + Sex",
    mod1d = "Species + Size \n+Sex",
    mod1e = "Species + \nSize * Sex",
    mod1f = "Species * Sex \n+ Size",
    mod1g = "Species * Size \n+ Sex",
    mod2 = "Species",
    mod2b = "Species + Size",
    mod2c = "Species + Sex",
    mod2d = "Species + Size \n+Sex",
    mod2e = "Species + \nSize * Sex",
    mod2f = "Species * Sex \n+ Size",
    mod2g = "Species * Size \n+ Sex"))

model.table <- sel.all %>%
  gt(rowname_col = "Model", groupname_col = "Trait", row_group_as_column = TRUE) %>%
  #tab_stubhead(label = "Model") %>%
  cols_label(
    Model = "Model",
    df = "df",
    AICc = "AICc",
    delta = "ΔAICc",
    weight = "Weight") %>%
  fmt_number(columns = c(AICc, delta), decimals = 2) %>%
  fmt_number(columns = weight, decimals = 3) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(df, AICc, delta)))

model.table

gtsave(model.table, "OnlineResource3_model-sel.rtf")


# ----------------------------------------------------------------------------------
### CTmax & CTmin main model diagnostics ###
# ----------------------------------------------------------------------------------

library(DHARMa)

# Simulate residuals once
sim.max <- simulateResiduals(mod1, n = 1000)
sim.min <- simulateResiduals(mod2, n = 1000)

build_diagnostic_plot <- function() {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Three rows, two columns
  par(
    mfrow = c(3, 2),
    mar = c(4, 4, 3.5, 1),
    cex = 0.75,
    family = "sans")
  
  # A: CTmax QQ
  plotQQunif(
    sim.max,
    main = "A. CTmax: species-only model checks")
  
  # B: CTmin QQ
  plotQQunif(
    sim.min,
    main = "B. CTmin: species-only model checks")
  
  # C: CTmax residuals vs predicted
  plotResiduals(sim.max)
  
  # D: CTmin residuals vs predicted
  plotResiduals(sim.min)
  
  # E: CTmax residuals vs species
  plotResiduals(sim.max, form = as.factor(model.frame(mod1)$Species))
  
  # F: CTmin residuals vs species
  plotResiduals(sim.min, form = as.factor(model.frame(mod2)$Species))
}

# TIFF
tiff(
  "OR3_mainResidChecks.tiff",
  width = 6.85,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo")
build_diagnostic_plot()
dev.off()

# EPS

cairo_ps(
  "OR3_mainResidChecks.eps",
  width = 6.85,
  height = 7,
  onefile = FALSE,
  fallback_resolution = 600)
build_diagnostic_plot()
dev.off()

# ----------------------------------------------------------------------------------
### CTmax & CTmin alt model diagnostics ###
# ----------------------------------------------------------------------------------

# Simulate residuals 
sim.size.max <- simulateResiduals(mod1b, n = 1000)
sim.sex.max <- simulateResiduals(mod1c, n = 1000)
sim.sss.max <- simulateResiduals(mod1d, n = 1000)
sim.size.min <- simulateResiduals(mod2b, n = 1000)
sim.sex.min <- simulateResiduals(mod2c, n = 1000)
sim.sss.min <- simulateResiduals(mod2d, n = 1000)

# ----------------------------
#### SPECIES + SIZE MODEL ####
# ----------------------------

build_SIZEdiagnostic_plot <- function() {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Four rows, two columns
  par(
    mfrow = c(4, 2),
    mar = c(4, 4, 3.5, 1),
    cex = 0.75,
    family = "sans")
  
  # A: CTmax QQ
  plotQQunif(
    sim.size.max,
    main = "A. CTmax: species+size model checks")
  
  # B: CTmin QQ
  plotQQunif(
    sim.size.min,
    main = "B. CTmin: species+size model checks")
  
  # C: CTmax residuals vs predicted
  plotResiduals(sim.size.max)
  
  # D: CTmin residuals vs predicted
  plotResiduals(sim.size.min)
  
  # E: CTmax residuals vs species
  plotResiduals(sim.size.max, form = as.factor(model.frame(mod1b)$Species))
  
  # F: CTmin residuals vs species
  plotResiduals(sim.size.min, form = as.factor(model.frame(mod2b)$Species))
  
  # G: CTmax residuals vs size
  plotResiduals(sim.size.max, form = model.frame(mod1b)$Size.mm)
  
  # H: CTmin residuals vs size
  plotResiduals(sim.size.min, form = model.frame(mod2b)$Size.mm)
}

# TIFF
tiff(
  "OR3_SizemodResidChecks.tiff",
  width = 6.85,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo")
build_SIZEdiagnostic_plot()
dev.off()

# EPS

cairo_ps(
  "OR3_SizemodResidChecks.eps",
  width = 6.85,
  height = 7,
  onefile = FALSE,
  fallback_resolution = 600)
build_SIZEdiagnostic_plot()
dev.off()

# ---------------------------
#### SPECIES + SEX MODEL ####
# ---------------------------

build_SEXdiagnostic_plot <- function() {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Four rows, two columns
  par(
    mfrow = c(4, 2),
    mar = c(4, 4, 3.5, 1),
    cex = 0.75,
    family = "sans")
  
  # A: CTmax QQ
  plotQQunif(
    sim.sex.max,
    main = "A. CTmax: species+sex model checks")
  
  # B: CTmin QQ
  plotQQunif(
    sim.sex.min,
    main = "B. CTmin: species+sex model checks")
  
  # C: CTmax residuals vs predicted
  plotResiduals(sim.sex.max)
  
  # D: CTmin residuals vs predicted
  plotResiduals(sim.sex.min)
  
  # E: CTmax residuals vs species
  plotResiduals(sim.sex.max, form = as.factor(model.frame(mod1c)$Species))
  
  # F: CTmin residuals vs species
  plotResiduals(sim.sex.min, form = as.factor(model.frame(mod2c)$Species))
  
  # G: CTmax residuals vs sex
  plotResiduals(sim.sex.max, form = as.factor(model.frame(mod1c)$Sex))
  
  # H: CTmin residuals vs sex
  plotResiduals(sim.sex.min, form = as.factor(model.frame(mod2c)$Sex))
}

# TIFF
tiff(
  "OR3_SexmodResidChecks.tiff",
  width = 6.85,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo")
build_SEXdiagnostic_plot()
dev.off()

# EPS

cairo_ps(
  "OR3_SexmodResidChecks.eps",
  width = 6.85,
  height = 7,
  onefile = FALSE,
  fallback_resolution = 600)
build_SEXdiagnostic_plot()
dev.off()

# -----------------------------------
#### SPECIES + SIZE + SEX MODELS ####
# -----------------------------------

build_SSSdiagnostic_plot <- function() {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Five rows, two columns
  par(
    mfrow = c(5, 2),
    mar = c(4, 4, 3.5, 1),
    cex = 0.75,
    family = "sans")
  
  # A: CTmax QQ
  plotQQunif(
    sim.sss.max,
    main = "A. CTmax: species+size+sex model checks")
  
  # B: CTmin QQ
  plotQQunif(
    sim.sss.min,
    main = "B. CTmin: species+size+sex model checks")
  
  # C: CTmax residuals vs predicted
  plotResiduals(sim.sss.max)
  
  # D: CTmin residuals vs predicted
  plotResiduals(sim.sss.min)
  
  # E: CTmax residuals vs species
  plotResiduals(sim.sss.max, form = as.factor(model.frame(mod1d)$Species))
  
  # F: CTmin residuals vs species
  plotResiduals(sim.sss.min, form = as.factor(model.frame(mod2d)$Species))
  
  # G: CTmax residuals vs size
  plotResiduals(sim.sss.max, form = model.frame(mod1d)$Size.mm)
  
  # H: CTmin residuals vs size
  plotResiduals(sim.sss.min, form = model.frame(mod2d)$Size.mm)
  
  # I: CTmax residuals vs sex
  plotResiduals(sim.sss.max, form = as.factor(model.frame(mod1d)$Sex))
  
  # J: CTmin residuals vs sex
  plotResiduals(sim.sss.min, form = as.factor(model.frame(mod2d)$Sex))
}

# TIFF
tiff(
  "OR3_CombomodResidChecks.tiff",
  width = 6.85,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo")
build_SSSdiagnostic_plot()
dev.off()

# EPS

cairo_ps(
  "OR3_CombomodResidChecks.eps",
  width = 6.85,
  height = 7,
  onefile = FALSE,
  fallback_resolution = 600)
build_SSSdiagnostic_plot()
dev.off()

# ----------------------------------------------------------------------------------
### CTmax & CTmin swordtail model selection ###
# ----------------------------------------------------------------------------------
SW.max.sel <- model.sel(SWmod1, SWmod2, SWmod3)
SW.min.sel <- model.sel(SWmod4, SWmod5, SWmod6)

sel.sw <- bind_rows(
  prepare.sel(SW.min.sel, "CTmin"),
  prepare.sel(SW.max.sel, "CTmax")) %>%
  mutate(Trait = factor(Trait, levels = c("CTmin", "CTmax"))) %>%
  mutate(Model = dplyr::recode(
    as.character(Model),
    SWmod1 = "Sex category",
    SWmod2 = "Sex category + Size",
    SWmod3 = "Sex category * Size",
    SWmod4 = "Sex category",
    SWmod5 = "Sex category + Size",
    SWmod6 = "Sex category * Size"))

SW.model.table <- sel.sw %>%
  gt(rowname_col = "Model", groupname_col = "Trait", row_group_as_column = TRUE) %>%
  #tab_stubhead(label = "Model") %>%
  cols_label(
    Model = "Model",
    df = "df",
    AICc = "AICc",
    delta = "ΔAICc",
    weight = "Weight") %>%
  fmt_number(columns = c(AICc, delta), decimals = 2) %>%
  fmt_number(columns = weight, decimals = 3) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(df, AICc, delta)))

SW.model.table

gtsave(SW.model.table, "OnlineResource3_SWmodel-sel.rtf")

# ----------------------------------------------------------------------------------
### Swordtail model diagnostics ###
# ----------------------------------------------------------------------------------


# Simulate residuals once
SW.sim.max <- simulateResiduals(SWmod1, n = 1000)
SW.sim.min <- simulateResiduals(SWmod4, n = 1000)

build_SWdiagnostic_plot <- function() {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Three rows, two columns
  par(
    mfrow = c(3, 2),
    mar = c(4, 4, 3.5, 1),
    cex = 0.75,
    family = "sans")
  
  # A: CTmax QQ
  plotQQunif(
    SW.sim.max,
    main = "A. CTmax: sex category model checks")
  
  # B: CTmin QQ
  plotQQunif(
    SW.sim.min,
    main = "B. CTmin: sex category model checks")
  
  # C: CTmax residuals vs predicted
  plotResiduals(SW.sim.max)
  
  # D: CTmin residuals vs predicted
  plotResiduals(SW.sim.min)
  
  # E: CTmax residuals vs species
  plotResiduals(SW.sim.max, form = as.factor(model.frame(SWmod1)$Sex.plot))
  
  # F: CTmin residuals vs species
  plotResiduals(SW.sim.min, form = as.factor(model.frame(mod2)$Sex.plot))
}

# TIFF
tiff(
  "OR3_SWResidChecks.tiff",
  width = 6.85,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo")
build_SWdiagnostic_plot()
dev.off()

# EPS

cairo_ps(
  "OR3_SWResidChecks.eps",
  width = 6.85,
  height = 7,
  onefile = FALSE,
  fallback_resolution = 600)
build_SWdiagnostic_plot()
dev.off()


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Online Resource 4 
# ───────────────────────────────────────────────────────────────────────────────────────────────

# Figures and tables for the outlier-removed analyses can be produced in Rerun-outlier-removal.Rmd

# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────
### ADDITIONAL TABLES/FIGURES NOT INCLUDED IN MANUSCRIPT ###
# ───────────────────────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Male vs Female Sizes       NOTE: run Final_temperature-analysis.Rmd prior to processing table below
# ───────────────────────────────────────────────────────────────────────────────────────────────

### By species
m.f.sizes <- sizes.fmt %>%
  gt() %>%
  tab_header(title = "Average lengths (mm)") %>%
  cols_label(Species ="Species",
             female = "Female (SD)",
             male = "Male (SD)") %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(Species, female))) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())) %>%
  opt_horizontal_padding(scale=2)
m.f.sizes
### By sizes
m.f.sizes2 <- sizes.fmt2 %>%
  gt() %>%
  tab_header(title = "Average lengths (mm)") %>%
  cols_label(female = "Female (sd)",
             male = "Male (sd)") %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(female))) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())) %>%
  opt_horizontal_padding(scale=2)
m.f.sizes2



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Paired points      NOTE: run Final_temperature-analysis.Rmd prior to processing figure below
# ───────────────────────────────────────────────────────────────────────────────────────────────


datTot2 <- rawTempData[,-c(5:6, 10:11)]


df_wide <- datTot2 %>%
  rename(Min = CTmin.C, Max = CTmax.C)

sex_colors <- c("male" = "#0072B2", "female" = "#D55E00")
plot_species_dynamic <- function(df_sub) {
  df_sub <- df_sub %>% mutate(Sex = factor(Sex))
  # Calculate min/max for this specific species to define the dynamic scale
  min_val <- min(df_sub$Min, na.rm = TRUE)
  max_val <- max(df_sub$Max, na.rm = TRUE)
  range_min <- max(df_sub$Min, na.rm = TRUE) - min(df_sub$Min, na.rm = TRUE)
  range_max <- max(df_sub$Max, na.rm = TRUE) - min(df_sub$Max, na.rm = TRUE)
  
  # Dynamic transformation
  scale_ratio <- range_max / range_min
  min_to_max <- function(x) (x - min_val) * (range_max / range_min) + min(df_sub$Max)
  max_to_min <- function(x) (x - min(df_sub$Max)) * (range_min / range_max) + min_val 
  # Create jitter for X-axis (the same for both points and segments)
  set.seed(123) # Ensures the jitter is the same every time you run it
  df_sub$x_jitter <- runif(nrow(df_sub), -0.03, 0.03)
  
  ggplot(df_sub, aes(group = Name)) +
    # Draw segments: X is 1+jitter, Xend is 1.5+jitter (tightened spacing)
    geom_segment(aes(x = 1 + x_jitter, xend = 1.2 + x_jitter, 
                     y = Max, yend = min_to_max(Min)), 
                 color = "gray24", alpha = 0.5) +
    # Points
    geom_point(aes(x = 1 + x_jitter, y = Max, shape = Sex, color=Sex), size = 2.5, alpha=0.8) +
    geom_point(aes(x = 1.2 + x_jitter, y = min_to_max(Min), shape = Sex, color=Sex), size = 2.5, alpha=0.8) +
    scale_color_manual(values = sex_colors) +
    scale_shape_manual(values = c(16, 17)) + # 16=Circle, 17=Triangle
    scale_y_continuous(
      name = "CTmax (\u00B0C)",
      sec.axis = sec_axis(~max_to_min(.), name = "CTmin (\u00B0C)")
    ) +
    # Set X-axis to match the new 1 and 1.5 positions
    scale_x_continuous(breaks = c(1, 1.2), labels = c("CTmax", "CTmin"), limits = c(0.9, 1.3), expand = c(0.01, 0.01)) +
    labs(title = unique(df_sub$Species)) +
    theme_classic(base_size = 10, base_family = "sans") + 
    theme(
      axis.title.y = element_text(size = 10, family = "sans"), 
      axis.text.y  = element_text(size = 10, family = "sans"),
      axis.text.x  = element_text(size = 10, family = "sans"),
      strip.text    = element_text(size = 10, family = "sans")) +
    theme(axis.title.x = element_blank())
}

# 2. Apply to each species
p_sword <- plot_species_dynamic(df_wide %>% filter(Species == "Swordtail"))
p_sail  <- plot_species_dynamic(df_wide %>% filter(Species == "Sailfin"))
p_mosq  <- plot_species_dynamic(df_wide %>% filter(Species == "Mosquitofish"))

# 3. Stack them
(p_sword / p_sail / p_mosq) + plot_layout(guides = "collect")







