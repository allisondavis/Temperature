# ───────────────────────────────────────────────────────────────────────────────────────────────
###### MANUSCRIPT TABLES AND FIGURES ######   NOTE: must run all Rmd files prior to this script
# ───────────────────────────────────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────────────────────────────────
### TABLES ###    NOTE: all tables created using gt() then exported to Word for final editing
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
  "Xiphophorus nigrensis", "Feeding and Growth", "Females mature earlier than intermediate and large males, while body size continues to increase post-maturity", "Ryan 1988")

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
## Table 2
# ───────────────────────────────────────────────────────────────────────────────────────────────

wt  <- wt  %>% mutate(Test = "CTmax")
wt2 <- wt2 %>% mutate(Test = "CTmin -0.3")
wt3 <- wt3 %>% mutate(Test = "CTmin -0.1")


wt.all <- bind_rows(wt, wt2, wt3)
wt.all <- wt.all %>%
  mutate(
    result = case_when(
      p.value < 0.001 ~ paste0("*p* < 0.001<br>(W = ", statistic,")"),
      p.value < 0.05  ~ paste0(
        "*p* = ", sprintf("%.3f", p.value),
        "<br>(W = ", statistic, ")"),
      TRUE ~ paste0(
        "*p* = ", sprintf("%.3f", p.value),
        "<br>(W = ", statistic, ")")))

wt.wide <- wt.all %>%
  select(Species, Test, result) %>%
  pivot_wider(
    names_from  = Test,
    values_from = result)

spp.trial.rates <- wt.wide %>%
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
## Table 3
# ───────────────────────────────────────────────────────────────────────────────────────────────

rate.pw <- pw_df %>%
  gt() %>%
  cols_label(
    group1 = "Species 1",
    group2 = "Species 2",
    p_adj  = md("*p*-value")) %>%
  fmt_markdown(columns = "p_adj") %>%
  fmt_number(
    columns = p_adj,
    decimals = 3) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(group1, group2))) %>%
  cols_align(align="center", everything()) %>%
  opt_horizontal_padding(scale=2)
rate.pw

gtsave(rate.pw, "Table3_rate-pairwise.rtf")


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 4
# ───────────────────────────────────────────────────────────────────────────────────────────────

rate.sum.table <- trial.rates2 %>%
  gt() %>%
  cols_label(Trial ="Trial type",
             Species = "Species",
             trial.rate = "Temp change (\u00B0C/min)",
             stand.dev = "Standard deviation",
             min = "Minimum (\u00B0C)",
             Median = "Median (\u00B0C)",
             max = "Maximum (\u00B0C)") %>%
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

gtsave(rate.sum.table, "Table4_rate-summaries.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 5
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


gtsave(annual_avg_rates, "Table5_austin-rates.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 6
# ───────────────────────────────────────────────────────────────────────────────────────────────

pw.ctmax <- pairs(PW.sp.max) %>%
  as.data.frame() %>%
  mutate(
    Trial = "CTmax",
    est_se = sprintf("%.2f (%.2f)", estimate, SE)) %>%
  select(Trial, contrast, est_se, p.value)

pw.ctmin <- pairs(PW.sp.min) %>%
  as.data.frame() %>%
  mutate(
    Trial = "CTmin",
    est_se = sprintf("%.2f (%.2f)", estimate, SE)) %>%
  select(Trial, contrast, est_se, p.value)


pw.all <- dplyr::bind_rows(pw.ctmax, pw.ctmin) %>%
  mutate(
    p.display = case_when(
      p.value < 0.001 ~ "*p* < 0.001",
      p.value < 0.05  ~ paste0(sprintf("%.3f", p.value)),
      TRUE            ~ sprintf("%.3f", p.value)))


spp.pw <- pw.all %>%
  gt(groupname_col = "Trial") %>%
  cols_hide(p.value) %>%
  cols_label(
    contrast = "Species contrast",
    est_se  = "Estimate (SE)",
    p.display = md("*p*-value")) %>%
  fmt_markdown(columns = p.display) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(contrast, est_se))) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.display,
      rows = p.value < 0.05)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  cols_align("center", everything()) %>%
  opt_horizontal_padding(scale = 2)
spp.pw

gtsave(spp.pw, "Table6_species-pairwise.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 7
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


gtsave(sum.stats.tab, "TB7_sum-stats.rtf")




# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 8
# ───────────────────────────────────────────────────────────────────────────────────────────────

pw.sex.max <- pairs(PW.sex.max) %>%
  as.data.frame() %>%
  mutate(
    Trial = "CTmax",
    est_se = sprintf("%.2f (%.2f)", estimate, SE)) %>%
  select(Trial, contrast, est_se, p.value)

pw.sex.min <- pairs(PW.sex.min) %>%
  as.data.frame() %>%
  mutate(
    Trial = "CTmin",
    est_se = sprintf("%.2f (%.2f)", estimate, SE)) %>%
  select(Trial, contrast, est_se, p.value)


pw.sex.all <- dplyr::bind_rows(pw.sex.max, pw.sex.min) %>%
  mutate(
    p.display = case_when(
      p.value < 0.001 ~ "*p* < 0.001",
      p.value < 0.05  ~ paste0(sprintf("%.3f", p.value)),
      TRUE            ~ sprintf("%.3f", p.value)))


sex.pw <- pw.sex.all %>%
  gt(groupname_col = "Trial") %>%
  cols_hide(p.value) %>%
  cols_label(
    contrast = "Strategies contrast",
    est_se  = "Estimate (SE)",
    p.display = md("*p*-value")) %>%
  fmt_markdown(columns = p.display) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(contrast, est_se))) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.display,
      rows = p.value < 0.05)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  cols_align("center", everything()) %>%
  opt_horizontal_padding(scale = 2)
sex.pw

gtsave(sex.pw, "Table8_sex-pairwise.rtf")


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 9
# ───────────────────────────────────────────────────────────────────────────────────────────────

df <- data.frame(
  type = c("CTmax", "CTmin"),
  Actual   = c(0.3, 0.2),
  Simulated   = c(0.7, 0.9))

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

gtsave(power.an, "Table9_power-analysis.rtf")


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table 10
# ───────────────────────────────────────────────────────────────────────────────────────────────

stats_tbl <- eff.size.tbl %>%
  left_join(size.test, by = "Species") %>%
  mutate(
    cliffs_delta = sprintf("%.2f", delta),
    p.display = case_when(
      p_value < 0.001 ~ "***p*** **<** **0.001**",
      p_value < 0.05  ~ paste0("*p* = ", sprintf("%.3f", p_value)),
      TRUE            ~ paste0("*p* = ", sprintf("%.3f", p_value))),
    
    wilcox_res = paste0(p.display, " (W = ", sprintf("%.0f", stat), ")")) %>%
  select(Species, cliffs_delta, wilcox_res)

size.table <- sizes.fmt %>%
  left_join(stats_tbl, by="Species")

size.tbl <- size.table %>%
  gt() %>%
  tab_header(title = "Average lengths (mm)") %>%
  cols_label(
    female       = "Female (SD)",
    male         = "Male (SD)",
    cliffs_delta = "Cliff’s δ",
    wilcox_res   = "Wilcoxon") %>%
  fmt_markdown(columns = wilcox_res) %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(Species, female, male, cliffs_delta))) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())) %>%
  opt_horizontal_padding(scale = 2)
size.tbl

gtsave(size.tbl, "Table10_size.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
### FIGURES ###   NOTE: Figure 1 created in Inkscape
# ───────────────────────────────────────────────────────────────────────────────────────────────
### Used to imbed fonts

#install.packages("ragg")
library(ragg)
#install.packages("systemfonts")
library(systemfonts)
# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 2
# ───────────────────────────────────────────────────────────────────────────────────────────────

### CTmax - EMM
df.pvals <- data.frame(
  group1 = c("Swordtail", "Sailfin", "Swordtail"),
  group2 = c("Sailfin", "Mosquitofish", "Mosquitofish"),
  label  = c("p<0.001*", "ns", "p < 0.001*"),
  y.position = c(39, 39, 39))

df.pvals$fontface <- ifelse(df.pvals$label == "ns", "plain", "bold")

robVCV.sp <- vcovHC(mod1f, type= "HC3")
PW.sp.max <- emmeans(mod1f, "Species", vcov. = robVCV.sp)
PW.df <- as.data.frame(PW.sp.max)

emmean.max<-ggplot(PW.df, aes(x=Species, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  stat_pvalue_manual(
    df.pvals,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5,
    step.increase = 0.1,
    fontface = "fontface") +
  ylab("Predicted mean temperature (\u00B0C) \n") +
  xlab("\n Species") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
emmean.max

### CTmin-EMM
robVCV.sp <- vcovHC(mod2f, type= "HC3")
PW.sp.min <- emmeans(mod2f, "Species", vcov. = robVCV.sp)
PW.df <- as.data.frame(PW.sp.min)

df.pvals <- data.frame(
  group1 = c("Swordtail", "Sailfin", "Swordtail"),
  group2 = c("Sailfin", "Mosquitofish", "Mosquitofish"),
  label  = c("p<0.001*", "p = 0.004*", "p<0.001*"),
  y.position = c(10.5, 10.5, 10.5))

df.pvals$fontface <- ifelse(df.pvals$label == "ns", "plain", "bold")

emmean.min <- ggplot(PW.df, aes(x=Species, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  ylab("Predicted mean temperature (\u00B0C) \n") +
  stat_pvalue_manual(
    df.pvals,
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

### CTmax-range x sex
dat.plot <- dat.plot %>%
  mutate(Sex.plot = case_when(
    Sex == "male" ~ "Male",
    Sex == "female" ~ "Female"))

range.max<-ggplot(dat.plot, aes(x = Sex.plot, y = CTmax.C)) +
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
  ylab("Temperature at LOE (\u00B0C)")
range.max

### CTmin-range x sex
range.min <- ggplot(dat.plot, aes(x = Sex.plot, y = CTmin.C)) +
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
  ylab("Temperature at LOE (\u00B0C)")
range.min

### Combining
range.spp <- (range.max + range.min + plot_layout(axis_titles = "collect"))

em.max.spp <- emmean.max + xlab("") + ggtitle("CTmax") + ylim(NA, 41)
em.min.spp <- emmean.min + xlab("") + ggtitle ("CTmin") + ylim(NA, 13)
em.spp <- (em.max.spp + em.min.spp + plot_layout(axis_titles = "collect"))

final.plot_em.spp<- ((em.spp/range.spp) + plot_annotation(tag_levels = "a"))
final.plot_em.spp

#Save as tiff
ggsave(
  filename = "Fig2_em-species.tiff",
  plot = final.plot_em.spp,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 3
# ───────────────────────────────────────────────────────────────────────────────────────────────

dat.plot3 <- rawTempData %>%
  mutate(
    Sex.plot3 = case_when(
      Sex == "male" ~ "Male",
      Sex == "female" ~ "Female"))

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
  filename = "Fig3_species-range.tiff",
  plot = range.plot,
  device = ragg::agg_tiff,
  width = 3.3, 
  height = 4,          
  units = "in",
  dpi = 600)

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 4
# ───────────────────────────────────────────────────────────────────────────────────────────────

### CTmax-EMM
PW.sex.max <- emmeans(SWmod1, "sex.2")
PW.df <- as.data.frame(PW.sex.max)
PW.df <- PW.df %>%
  mutate(
    Sex.plot = case_when(
      sex.2 == "male.coercion" ~ "Male (coercion)",
      sex.2 == "male.courtship" ~ "Male (courtship)",
      sex.2 == "female" ~ "Female")) 

em.max <- ggplot(PW.df, aes(x=Sex.plot, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ylab("Predicted mean temperature (\u00B0C) \n") +
  xlab("\n Mating strategy") +
  ggtitle("\n CTmax") +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% c("Male (coercion)", "Male (courtship)"),
             gsub(" ", "\n", x),
             x)}) +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
em.max

### CTmin-EMM
PW.sex.min <- emmeans(SWmod2, "sex.2")
PW.df <- as.data.frame(PW.sex.min)
PW.df <- PW.df %>%
  mutate(
    Sex.plot = case_when(
      sex.2 == "male.coercion" ~ "Male (coercion)",
      sex.2 == "male.courtship" ~ "Male (courtship)",
      sex.2 == "female" ~ "Female")) 

em.min <- ggplot(PW.df, aes(x=Sex.plot, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ylab("Predicted mean temperature (\u00B0C) \n") +
  xlab("\n Mating strategy") +
  ggtitle("CTmin") +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% c("Male (coercion)", "Male (courtship)"),
             gsub(" ", "\n", x),
             x)}) +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
em.min

### CTmax-range x category
dat.plot.sex <- swordtail_1 %>%
  mutate(
    Sex.plot = case_when(
      sex.2 == "male.coercion" ~ "Male (coercion)",
      sex.2 == "male.courtship" ~ "Male (courtship)",
      sex.2 == "female" ~ "Female")) 

max.plot <- ggplot(dat.plot.sex, aes(x = Sex.plot, y = CTmax.C)) +
  geom_boxplot() +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% c("Male (coercion)", "Male (courtship)"),
             gsub(" ", "\n", x),
             x)}) +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans")) +
  xlab("\n Mating strategy") +
  ylab("Temperature at LOE (\u00B0C)")
max.plot

### CTmin-range x category
min.plot <- ggplot(dat.plot.sex, aes(x = Sex.plot, y = CTmin.C)) +
  geom_boxplot() +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  scale_x_discrete(
    labels = function(x) {
      ifelse(x %in% c("Male (coercion)", "Male (courtship)"),
             gsub(" ", "\n", x),
             x)}) +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans")) +
  xlab("\n Mating strategy") +
  ylab("Temperature at LOE (\u00B0C)")
min.plot

### Combining
range.swd <-  (max.plot + min.plot + plot_layout(axis_titles = "collect"))

em.max <- em.max  + xlab("")
em.min <- em.min  + xlab("")

em.swd <- (em.max + em.min + plot_layout(axis_titles = "collect"))

final.plot_em.sex <- ((em.swd/range.swd) + plot_annotation(tag_levels = "a"))
final.plot_em.sex

#Save as tiff
ggsave(
  filename = "Fig4_em-sex.tiff",
  plot = final.plot_em.sex,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure 5
# ───────────────────────────────────────────────────────────────────────────────────────────────

df.pval <- data.frame(label = c("ns", "ns", "p<0.001*"),
                      Species = c("Swordtail", "Sailfin", "Mosquitofish"),
                      x = c(1.5, 1.5, 1.5),
                      y = c(39, 39, 39))

bracket_data <- data.frame(
  Species = "Mosquitofish", 
  x = c(1, 1, 2, 2),
  y = c(37.5, 38, 38, 37.5))


dat.plot2 <- datTot %>%
  mutate(
    Sex.plot2 = case_when(
      Sex == "male" ~ "Male",
      Sex == "female" ~ "Female"))


size.plot <- ggplot(dat.plot2, aes(x = Sex.plot2, y = Size.mm)) +
  geom_rect(
    data = subset(datTot, Species == "Sailfin"),
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "grey95") +
  geom_boxplot() +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  facet_wrap(~ Species, scales = "free_x") +
  geom_text(df.pval, mapping = aes(x=x, y=y, label=label))+
  geom_path(
    data = bracket_data, 
    aes(x = x, y = y), 
    inherit.aes = FALSE, 
    linewidth = 0.5) +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Length (mm)")
size.plot

#Save as tiff
ggsave(
  filename = "Fig5_size.tiff",
  plot = size.plot,
  device = ragg::agg_tiff,
  width = 3.3, 
  height = 4,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
### SUPPLEMENTARY TABLES/FIGURES ###
# ───────────────────────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Table S1
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
      ifelse(grepl("max_slope", x), "Max",
             ifelse(grepl("min_slope", x), "Min", x))}) %>%
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

gtsave(extreme.table, "TableS1_extreme-rates-austin.rtf")



# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure S1
# ───────────────────────────────────────────────────────────────────────────────────────────────

#visualize CTmin outliers
out.min <- ggplot(df_flagged, aes(x = CTmin.C, fill = is_outlier_min)) +
  geom_histogram(binwidth = 0.5, color = "white") +
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02"),
                    name = "Data Status",
                    labels = c("Normal", "Outlier")) +
  labs(title = "CTmin",
       x = "Critical Thermal Minimum (\u00B0C)",
       y = "Frequency") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
out.min #outliers for sailfin

#visualize CTmax outliers
out.max <- ggplot(df_flagged, aes(x = CTmax.C, fill = is_outlier_max)) +
  geom_histogram(binwidth = 0.5, color = "white") +
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#d95f02"),
                    name = "Data Status",
                    labels = c("Normal", "Outlier")) +
  labs(title = "CTmax",
       x = "Critical Thermal Maximum (\u00B0C)",
       y = "Frequency") +
  theme_classic(base_size = 10, base_family = "sans") + 
  theme(
    axis.title.y = element_text(size = 10, family = "sans"), 
    axis.text.y  = element_text(size = 10, family = "sans"),
    axis.text.x  = element_text(size = 10, family = "sans"),
    strip.text    = element_text(size = 10, family = "sans"))
out.max #outliers for swordtail and mosquitofish

out.plot <- ((out.min/out.max) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect"))
out.plot

#Save as tiff
ggsave(
  filename = "FigS1_outliers.tiff",
  plot = out.plot,
  device = ragg::agg_tiff,
  width = 6.85, 
  height = 7,          
  units = "in",
  dpi = 600)


# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure S2
# ───────────────────────────────────────────────────────────────────────────────────────────────
#install.packages("DHARMa")
library(DHARMa)

simulationOutput <- simulateResiduals(mod1f)
plot(simulationOutput) #get results and type in manually below

#QQ plot test statistics
ks_res   <- "KS test: p= 0.055"
disp_res <- "Dispersion test: p= 0.928"
out_res  <- "Outlier test: p= 0.026"

#Residuals plot test statistics
levene_res <- "Levene Test of homogeneity \nof variance n.s"
group_res <- "Within-group deviation from uniformity n.s."

build_diagnostic_plot <- function() {

# Set up a 1-row, 2-column layout with standard margins
par(family = "sans", lwd = 1, mfrow = c(1, 2), mar = c(5, 4.5, 4, 2))

# ───────────────────
# PANEL 1: QQ PLOT
# ───────────────────
# Turn off DHARMa's internal text printing
plotQQunif(simulationOutput, 
           testUniformity = FALSE, 
           testOutliers = FALSE, 
           testDispersion = FALSE,
           main = "")

# Add custom main title
title("QQ Residuals", cex.main = 0.9, line = 2.5)

# Place manual labels at the top of the plot area
## side = 3 places text in the top margin; line controls vertical distance
mtext(ks_res,   side = 3, line = 1.4, cex = 0.85, adj = 0)
mtext(disp_res, side = 3, line = 0.7, cex = 0.85, adj = 0)
mtext(out_res,  side = 3, line = 0.0, cex = 0.85, adj = 0, font = 2, col = "red")


# ──────────────────────────
# PANEL 2: RESIDUALS PLOT
# ──────────────────────────
#draw boxplots normally
plotResiduals(simulationOutput, main = "")

#draw a white rectangle to cover up DHARMa's text
##force the white box to start at the absolute left edge of the device region
x_left_edge <- grconvertX(0, from = "ndc", to = "user") 
x_limits    <- par("usr")[1:2] 

rect(xleft = x_left_edge,   # Anchors to the absolute left margin barrier
     ybottom = 1.01,        
     xright = x_limits[2] + 2, 
     ytop = 1.45,          
     col = "white", 
     border = NA, 
     xpd = TRUE)

# Add custom main title and stacked text over the clean white space
title("Residual vs. Predicted", cex.main = 0.9, line = 2.5)

mtext(group_res,  side = 3, line = 1.4, cex = 0.75, adj = 0)
mtext(levene_res, side = 3, line = 0.0, cex = 0.75, adj = 0)
} ### END FUNCTION

#save as tiff
tiff("FigS2_DHARMa_CTmax.tiff",
    width = 6.85,
    height = 4.5,
    units = "in",
    res = 600,
    compression = "lzw",
    type = "cairo")

build_diagnostic_plot()
dev.off()




# ───────────────────────────────────────────────────────────────────────────────────────────────
## Figure S3
# ───────────────────────────────────────────────────────────────────────────────────────────────
#install.packages("DHARMa")
library(DHARMa)

simulationOutput <- simulateResiduals(mod2f)
plot(simulationOutput) #get results and type in manually below

#QQ plot test statistics
ks_res   <- "KS test: p= 0.855"
disp_res <- "Dispersion test: p= 0.928"
out_res  <- "Outlier test: p= 1"

#Residuals plot test statistics
levene_res <- "Levene Test of homogeneity \nof variance significant"
group_res <- "Within-group deviation from uniformity n.s."

build_diagnostic_plot <- function() {

# Set up a 1-row, 2-column layout with standard margins
par(family = "sans", lwd = 1, mfrow = c(1, 2), mar = c(5, 4.5, 4, 2))

# ───────────────────
# PANEL 1: QQ PLOT
# ───────────────────
# Turn off DHARMa's internal text printing
plotQQunif(simulationOutput, 
           testUniformity = FALSE, 
           testOutliers = FALSE, 
           testDispersion = FALSE,
           main = "")

# Add custom main title
title("QQ Residuals", cex.main = 0.9, line = 2.5)

# Place manual labels at the top of the plot area
## side = 3 places text in the top margin; line controls vertical distance
mtext(ks_res,   side = 3, line = 1.4, cex = 0.85, adj = 0)
mtext(disp_res, side = 3, line = 0.7, cex = 0.85, adj = 0)
mtext(out_res,  side = 3, line = 0.0, cex = 0.85, adj = 0)


# ──────────────────────────
# PANEL 2: RESIDUALS PLOT
# ──────────────────────────
#draw boxplots normally
plotResiduals(simulationOutput, main = "")

#draw a white rectangle to cover up DHARMa's text
##force the white box to start at the absolute left edge of the device region
x_left_edge <- grconvertX(0, from = "ndc", to = "user") 
x_limits    <- par("usr")[1:2] 

rect(xleft = x_left_edge,   # Anchors to the absolute left margin barrier
     ybottom = 1.01,        
     xright = x_limits[2] + 2, 
     ytop = 1.45,          
     col = "white", 
     border = NA, 
     xpd = TRUE)

# Add custom main title and stacked text over the clean white space
title("Residual vs. Predicted", cex.main = 0.9, line = 2.5)

mtext(group_res,  side = 3, line = 1.4, cex = 0.75, adj = 0)
mtext(levene_res, side = 3, line = 0.0, cex = 0.75, adj = 0, font = 2, col = "red")
} ### END FUNCTION


#save as tiff
tiff("FigS3_DHARMa_CTmin.tiff",
    width = 6.85,
    height = 4.5,
    units = "in",
    res = 600,
    compression = "lzw",
    type = "cairo")

build_diagnostic_plot()
dev.off()


# ───────────────────────────────────────────────────────────────────────────────────────────────
### ADDITIONAL TABLES/FIGURES NOT INCLUDED IN MANUSCRIPT ###
# ───────────────────────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────────────────────
## Male vs Female Sizes
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
## Paired points
# ───────────────────────────────────────────────────────────────────────────────────────────────


datTot2 <- rawTempData[,-c(6:7, 10)]


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







