---
title: "README"
---


## **Out of the frying pan and into the freezer: sex differences in thermal tolerance lacking in livebearing species**

AUTHORS: Sophia E. McKelvey, Cassidy Hawk, Eric N. Iverson, Michael J. Ryan, and Allison C. Davis


***
***

### **EXPERIMENTAL OVERVIEW:**

**Purpose:** Our interest lies in whether sex-biased temperature tolerance exists in 3 livebearing fishes with sexual size dimorphism: sailfin mollies (*Poecilia latipinna*), Western mosquitofish (*Gambusia affinis*), and Panuco swordtails (*Xiphophorous nigrensis*). This stems from an observation of males dying when sampling during hot months but females surviving fine. 

**Main Aims:**

(1) compare the effects of sex, size and their interaction on thermal tolerance across these three species of livebearing fishes (Family: Poeciliidae)

(2) compare thermal tolerance across male morphs within a single livebearing species (*X. nigrensis*)

**Hypothesis:** Preexisting sex differences (e.g. size) decrease resistance to environmental temperature change.

***
***

### **QUICK REPOSITORY NAVIGATION**

* **Austin Daily Temperature/**
  Daily climatology analyses and associated datasets

* **Trial Rates/**
  Trial ramping rate calculations

* **Thermal Tolerance Analysis/**
  Primary statistical analyses of sex and size differences in thermal tolerance

* **Outlier Removed Analysis/**
  Re-analysis excluding major outliers

* **Extras/**
  Raw trial Excel files, climatology source data, and manuscript figure/table scripts

***
***

### **REPRODUCING ANALYSES**

1. Download and unzip all folders into the same parent directory.

2. Open the R project in RStudio.

3. Each analysis folder is independent and can be run separately unless otherwise noted.

   * *Rerun-outlier-removal.Rmd* requires the *Temp ID - Completed Data.csv* in the **Thermal Tolerance Analysis** folder.

4. To recreate manuscript figures and tables, run:

   * *Austin-Daily-Temps.Rmd*
   * *trial_rate_calc.Rmd*
   * *Final_temperature-analysis.Rmd*

   prior to running:

   * *manuscript_figures-tables.R*
   
   Note: one figure (Fig S1) requires *Rerun-outlier-removal.Rmd*; may overwrite data/values from *Final_temperature-analysis.Rmd* in the environment.

5. All file paths are relative and should function provided the original folder structure is maintained.

***
***

### **FILE STRUCTURE & DATA DESCRIPTION**

The four folders correspond to the four analyses used in this study (all Rmd files are annotated):

**1. Austin Daily Temperature:** calculating the daily temperature change rates (°C/min) for 2020-2024 in Austin, TX USA to compare to the temperature rate changes used in our trials.

  \- *Austin-Daily-Temps.Rmd*

  \- *Austin_temp-data.csv*

| VARIABLE | DEFINITION | UNIT/OPTIONS |
| :--- | :--- | :--- | 
| Station | station ID for Camp Mabry, Austin TX USA | USW00013958 |
| Date | recording date | M/D/YYYY |
| Hour | time of recording | H:MM:SS |
| Report_type | type of meteorological event | FM-12/FM-15/FM-16 |
| Source | data collection sub-system | 220/223/343 |
| HourlyDryBulbTemperature | air temperature | °C |


   * Type of meteorological event: 
      - FM-12: global World Meteorological Organization samples generated at synoptic hours (00:00, 06:00, 12:00, 18:00, UTC)
      - FM-15: standard scheduled hourly observations
      - FM-16: unscheduled special reporting triggered by rapid weather changes

   * Data collection sub-system:
      - 220: hourly temperature collection failed, data back filled using Automated Surface Observing System 5-minute Data archive
      - 223: data via United States Air Force Surface Weather Observations
      - 343: data via National Oceanic and Atmospheric Administration Surface Weather Observation
    
    
**2. Trial Rates:** calculating the ramping rate achieved in our trials to compare to the desired goal rate of our automatic Arduino system (± 0.3°C/minute)
    
  \- *trial_rate_calc.Rmd* 
  
  \- *rate_calc.csv*

| VARIABLE | DEFINITION | UNIT/OPTIONS |
| :--- | :--- | :--- | 
| Trial | trial type identifier (CTmax/CTmin) | max/min |
| Species | Poeciliid species of interest | gambusia/sailfin/swordtail |
| Size | size category (for male swordtails) | large/small/NA |
| Date | trial date | MM/DD/YYY |
| Rate | temperature rate calculated from original excel file | °C/min |
| file.id | name of original Excel file | NA |

   * Note: size in this data set refers to a trial's size category, which is pertinent only to swordtails as males were categorized as small or large; all other individuals (including female swordtails) marked as NA

   * Note: temperature rate calculation can be found on the first sheet of the original excel file
   

**3. Thermal Tolerance Analysis:** overall analysis testing for sex and size differences in thermal tolerance between three Poeciliid fish species
    
  \- *Final_temperature-analysis.Rmd*
  
  \- *Temp ID - Completed Data.csv*
  
| VARIABLE | DEFINITION | UNIT/OPTIONS |
| :--- | :--- | :--- | 
| Species | Poeciliid species of interest | X. nigrensis/G. affinis/P. latipinna |
| Size.cat | Size category (for male swordtails) | large/small/reg |
| Sex | Sex of individual | male/female |
| Name | Name chosen for tracking individuals | NA |
| CTmax.num | Number given when fish is removed in CTmax trials | 1-9 |
| CTmin.num | Number given when fish is removed in CTmin trials | 1-9 |
| CTmax.C | Temperature of tank at loss of equilibrium in CTmax trials | °C |
| CTmin.C | Temperature of tank at loss of equilibrium in CTmin trials | °C |
| Size.mm | Standard length of individual | mm |
| file.id.max | file name associated with CTmax trial run | NA |
| file.id.min | file name associated with CTmin trial run | NA |


**4. Outlier Removed Analysis:** a re-run of the overall analysis testing for sex and size differences in thermal tolerance between three Poeciliid fish species without major outliers
    
  \- *Rerun-outlier-removal.Rmd*
  
  \- Uses the *Temp ID - Completed Data.csv* located in the **Thermal Tolerance Analysis** folder


***
    
An additional folder -- **Extras** -- contains the raw excel files corresponding to each trial with additional trial notes, the raw climatology data for the Austin Daily temperature analysis, and an R script for creating the manuscript tables/figures.

  \- Source for climatology data: Kantor D, Casey NW, Menne MJ, Buddenberg A (2023) Local Climatological Data (LCD), Version 2 [Subset used: Austin Camp Mabry, TX US (USW00013958), hourly dry bulb temperature, 2020–2024]. NOAA National Centers for Environmental Information. https://doi.org/10.25921/jp3d-3v19. Accessed 8 April 2025

   * See details at the above DOI for structure of raw climatology data
  
  \- Structure of raw trial excel files ('Trial data sheets'):
   
   * Sheet 1: created post-trial to calculate the temperature rate (°C/min) experienced in that trial. Time column (number of seconds since the start of temperature tracking) and Temp column (°C of tank at associated time) were copied from the 'Data In' sheet. The slope formula (=SLOPE(B/A * 60)) was calculated in a free cell.
   
   * Data In: sheet displaying the Arduino temperature recordings. 'Time' is the time of the day (HH:MM:SS.MS), 'CH1' is the rate calculated every 40s, 'CH2' is the number of seconds since the start of temperature tracking, 'CH3' is the °C of tank at associated time. When the button on the Arduino system was pressed to record an individual, 'CH1' rate is replaced with organism number, 'CH2' is replaced with the temperature, and 'CH3' is replaced with the number of seconds since recording start. 
   
   * Data Out/Settings/Manifest: These sheets were automatically structured with the installation of the Arduino temperature system, and no alterations were made by us during our trials.

  \- *manuscript_figures-tables.R* (code for designing the figures/tables used in the manuscript)


***
***

### **VERSION & PACKAGE INFORMATION**

All Rmd files list required packages at the beginning of the document. Below are all packages required to run all Rmd files.
NOTE: packages required for *Rerun-outlier-removal.Rmd* are equivalent to *Final_temperature-analysis.Rmd*

| PACKAGE/SOFTWARE | VERSION | USED IN |
| :--- | :--- | :--- | 
| Rstudio running R | v4.6.1 | All Rmd files | 
| `curl` | v7.1.0 | All Rmd files | 
| `dplyr` | v1.2.1 | All Rmd files | 
| `tidyr` | v1.3.2 | All Rmd files | 
| `ggplot2` | v4.0.3 | All Rmd files | 
| `gt` | v1.3.0 | All Rmd files |
| `lubridate` | v1.9.5 | *Austin-Daily-Temps.Rmd* |
| `summarytools` | v1.1.5 | *trial_rate_calc.Rmd*/ *Final_temperature-analysis.Rmd* |
| `MuMIn` | v1.48.19 | *Final_temperature-analysis.Rmd* | 
| `lmerTest` | v3.2-1 | *Final_temperature-analysis.Rmd* | 
| `clubSandwich` | v3.1-1 | *Final_temperature-analysis.Rmd* | 
| `ggpubr` | v0.6.3 | *Final_temperature-analysis.Rmd* |
| `DHARMa` | v0.5.0 | *Final_temperature-analysis.Rmd*/ *manuscript_figures-tables.R* | 
| `emmeans` | v2.0.3 | *Final_temperature-analysis.Rmd* | 
| `car` | v3.1-5 | *Final_temperature-analysis.Rmd* | 
| `patchwork` | v1.3.2 | *Final_temperature-analysis.Rmd* | 
| `tibble` | v3.3.1 | *Final_temperature-analysis.Rmd*/ *manuscript_figures-tables.R* | 
| `ragg` | v1.5.2 | *manuscript_figures-tables.R* | 
| `systemfonts` | v1.3.2 | *manuscript_figures-tables.R* | 
| `lme4` | v2.0-1 | *manuscript_figures-tables.R* | 
| `scales` | v1.4.0 | *manuscript_figures-tables.R* | 





***
***


### **CITATION**

If using these data or scripts, please cite:

McKelvey SE, Hawk C, Iverson EN, Ryan MJ, Davis AC. 2026. *Out of the frying pan and into the freezer: sex differences in thermal tolerance lacking in livebearing species.* [Journal name]. DOI: XXXXX

Data repository archived at Zenodo: DOI: 10.5281/zenodo.20466137

Data and scripts are distributed under the Creative Commons Attribution 4.0 International License (CC BY 4.0) unless otherwise noted.


***
***






 
 
 
