############################################
###### MANUSCRIPT TABLES AND FIGURES ######   NOTE: run Final_temperature-analysis.Rmd prior to this script
############################################


###############
### TABLES ###    NOTE: all tables created using gt() then exported to Word for final editing
###############

#############
## Table 1 
#############

#############
## Table 2
#############

#############
## Table 3
#############

#############
## Table 4
#############

#############
## Table 5
#############

#############
## Table 6
#############
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


pw.all %>%
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




#############
## Table 7
#############
sum.stats.tab <- summ.stats.all %>%
  arrange(Species, Trial) %>%
  select(Species, Trial, mean_sd, min, med, max, male, female) %>%
  gt(rowname_col = "Trial", groupname_col = "Species") %>%
  cols_label(
    Species = "Species",
    Trial   = "Trial",
    mean_sd = "Mean (°C) (SD)",
    min     = "Minimum (°C)",
    med     = "Median (°C)",
    max     = "Maximum (°C)",
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




#############
## Table 8
#############
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


pw.sex.all %>%
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




#############
## Table 9
#############
df <- data.frame(
  type = c("CTmax", "CTmin"),
  Actual   = c(0.3, 0.2),
  Simulated   = c(0.7, 0.9))

df %>%
  gt() %>%
  cols_label(type ="Temperature Category",
             Actual = "Actual difference (°C)",
             Simulated = "Simulated difference needed (°C)") %>%
  tab_style(
    style = cell_borders(sides = "right", color = "lightgray", weight = px(1.5)),
    locations = cells_body(columns = c(type, Actual))) %>%
  cols_align(align = "center", everything())  %>%
  opt_horizontal_padding(scale=2)




#############
## Table 10
#############
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

size.table %>%
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




################
### FIGURES ###   NOTE: Figure 1 created in Inkscape
################

#############
## Figure 2
#############
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
  ylab("Predicted mean temperature (°C) \n") +
  xlab("\n Species") +
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial"))
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
  ylab("Predicted mean temperature (°C) \n") +
  stat_pvalue_manual(
    df.pvals,
    label = "label",
    tip.length = 0.01,
    bracket.size = 0.5,
    step.increase = 0.1,
    fontface = "fontface") +
  xlab("\n Species") +
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial"))
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
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("Sex") +
  ylab("Temperature at LOE (°C)")
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
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("Sex") +
  ylab("Temperature at LOE (°C)")
range.min

### Combining
range.spp <- (range.max + range.min + plot_layout(axis_titles = "collect"))

em.max.spp <- emmean.max + xlab("") + ggtitle("CTmax") + ylim(NA, 41)
em.min.spp <- emmean.min + xlab("") + ggtitle ("CTmin") + ylim(NA, 13)
em.spp <- (em.max.spp + em.min.spp + plot_layout(axis_titles = "collect"))

((em.spp/range.spp) + plot_annotation(tag_levels = "a"))




#############
## Figure 3
#############
dat.plot3 <- rawTempData %>%
  mutate(
    Sex.plot3 = case_when(
      Sex == "male" ~ "Male",
      Sex == "female" ~ "Female"))

ggplot(dat.plot3, aes(x = Sex.plot3, y = T.tot)) +
  geom_rect(
    data = subset(datTot, Species == "Sailfin"),
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "grey95") +
  geom_boxplot() +
  geom_point(color= "#8D8680", size=2, alpha=0.5) +
  facet_wrap(~ Species, scales = "free_x") +
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("") +
  ylab("Temperature range (°C) \n(CTmax-CTmin)")




#############
## Figure 4
#############
### CTmax-EMM
PW.sex.max <- emmeans(SWmod1, "sex.2")
PW.df <- as.data.frame(PW.sex.max)
em.max <- ggplot(PW.df, aes(x=sex.2, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ylab("Predicted mean temperature (°C) \n") +
  xlab("\n Mating strategy") +
  ggtitle("\n CTmax") +
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial"))
em.max

### CTmin-EMM
PW.sex.min <- emmeans(SWmod2, "sex.2")
PW.df <- as.data.frame(PW.sex.min)
em.min <- ggplot(PW.df, aes(x=sex.2, y=emmean)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ylab("Predicted mean temperature (°C) \n") +
  xlab("\n Mating strategy") +
  ggtitle("CTmin") +
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial"))
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
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("\n Mating strategy") +
  ylab("Temperature at LOE (°C)")
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
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("\n Mating strategy") +
  ylab("Temperature at LOE (°C)")
min.plot

### Combining
range.swd <-  (max.plot + min.plot + plot_layout(axis_titles = "collect"))

my.labels <- c("Female", "Male, \n courtship", "Male, \n coercion")
em.max <- em.max + scale_x_discrete(labels = my.labels) + xlab("")
em.min <- em.min + scale_x_discrete(labels = my.labels) + xlab("")

em.swd <- (em.max + em.min + plot_layout(axis_titles = "collect"))

((em.swd/range.swd) + plot_annotation(tag_levels = "a"))




#############
## Figure 5
#############
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


ggplot(dat.plot2, aes(x = Sex.plot2, y = Size.mm)) +
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
  theme_classic(base_size = 10, base_family = "Arial") + 
  theme(
    axis.title.y = element_text(size = 10, family = "Arial"), 
    axis.text.y  = element_text(size = 10, family = "Arial"),
    axis.text.x  = element_text(size = 10, family = "Arial"),
    strip.text    = element_text(size = 10, family = "Arial")) +
  xlab("") +
  ylab("Length (mm)")




#####################################
### SUPPLEMENTARY TABLES/FIGURES ###
#####################################

#############
## Table S1
#############

#############
## Figure S1
#############

#############
## Figure S2
#############
simulationOutput <- simulateResiduals(mod2f)
#png
png("FigS2_QQ.png",
    width = 3.3,
    height = 5,
    units = "in",
    res = 600)
par(family="sans", lwd = 1)
plotQQunif(simulationOutput)
dev.off()

png("FigS2_residuals.png",
    width = 3.3,
    height = 5,
    units = "in",
    res = 600)
par(family="sans", lwd = 1)
plotResiduals(simulationOutput)
dev.off()

#pdf
cairo_pdf("FigS2_QQ.pdf",
    width = 3.3,
    height = 5,
    family = "sans")
par(family="sans", cex =1, lwd = 1)
plotQQunif(simulationOutput)
dev.off()

cairo_pdf("FigS2_residuals.pdf",
    width = 3.3,
    height = 5,
    family = "sans")
par(family="sans", cex=1, lwd = 1)
plotResiduals(simulationOutput)
dev.off()

#############
## Figure S3
#############
simulationOutput <- simulateResiduals(mod2f)
#png
png("FigS3_QQ.png",
    width = 3.3,
    height = 5,
    units = "in",
    res = 600)
par(family="sans", lwd = 1)
plotQQunif(simulationOutput)
dev.off()

png("FigS3_residuals.png",
    width = 3.3,
    height = 5,
    units = "in",
    res = 600)
par(family="sans", lwd = 1)
plotResiduals(simulationOutput)
dev.off()

#pdf
cairo_pdf("FigS3_QQ.pdf",
          width = 3.3,
          height = 5,
          family = "sans")
par(family="sans", cex =1, lwd = 1)
plotQQunif(simulationOutput)
dev.off()

cairo_pdf("FigS3_residuals.pdf",
          width = 3.3,
          height = 5,
          family = "sans")
par(family="sans", cex=1, lwd = 1)
plotResiduals(simulationOutput)
dev.off()


#############################################################
### ADDITIONAL TABLES/FIGURES NOT INCLUDED IN MANUSCRIPT ###
#############################################################

##########################
## Male vs Female Sizes
##########################
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



#################
## Paired points
#################

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
      name = "CTmax (°C)",
      sec.axis = sec_axis(~max_to_min(.), name = "CTmin (°C)")
    ) +
    # Set X-axis to match the new 1 and 1.5 positions
    scale_x_continuous(breaks = c(1, 1.2), labels = c("CTmax", "CTmin"), limits = c(0.9, 1.3), expand = c(0.01, 0.01)) +
    labs(title = unique(df_sub$Species)) +
    theme_classic(base_size = 10, base_family = "Arial") + 
    theme(
      axis.title.y = element_text(size = 10, family = "Arial"), 
      axis.text.y  = element_text(size = 10, family = "Arial"),
      axis.text.x  = element_text(size = 10, family = "Arial"),
      strip.text    = element_text(size = 10, family = "Arial")) +
    theme(axis.title.x = element_blank())
}

# 2. Apply to each species
p_sword <- plot_species_dynamic(df_wide %>% filter(Species == "Swordtail"))
p_sail  <- plot_species_dynamic(df_wide %>% filter(Species == "Sailfin"))
p_mosq  <- plot_species_dynamic(df_wide %>% filter(Species == "Mosquitofish"))

# 3. Stack them
(p_sword / p_sail / p_mosq) + plot_layout(guides = "collect")







