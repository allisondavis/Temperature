---
title: "README"
---


## **Out of the frying pan and into the freezer: sex differences in thermal tolerance lacking in livebearing species**

AUTHORS: Sophia E. McKelvey, Cassidy Hawk, Eric N. Iverson, Michael J. Ryan, and Allison C. Davis

***
***

### **EXPERIMENTAL OVERVIEW:**

**Purpose:** Our interest lies in whether sex-biased temperature tolerance exists in 3 livebearing fishes: sailfin mollies (*Poecilia latipinna*), Western mosquitofish (*Gambusia affinis*), and Panuco swordtails (*Xiphophorous nigrensis*). This stems from an observation of males dying when sampling during hot months but females surviving fine. 

**Main Aims:**

(1) compare female and male thermal tolerance across these three species of livebearing fishes (Family: Poeciliidae)

(2) compare thermal tolerance across mating strategies within a single livebearing species (*X. nigrensis*)

**Hypothesis:** Preexisting sex differences (e.g. sexually selected traits, size, behavior etc.) decrease resistance to environmental temperature change.

**Predictions:**

  - Males of fish with elaborate traits have a narrower thermal tolerance--i.e. lower thermal maximums and higher thermal minimums--compared to females of these species.

  - Males of fish without elaborate traits will have a similar tolerance to females of the same species.

  - Males with traits used in courtship (i.e. bright coloration and elongated fins) have a narrower thermal tolerance compared to males that use coercion (within swordtails).

***
***

### **REPOSITORY NAVIGATION**

All four analysis folders are independent and do not require a sequential run order. All necessary data inputs are self-contained within their respective folders, with the exception of the *Rerun-outlier-removal.Rmd* that requires the data within the **Thermal Tolerance Analysis** folder.

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
    
    
**2. Trial Rates:** calculating the temperature rate increase and decrease achieved in our trials to compare to the desired goal rate of our automatic Arduino system (± 0.3°C/minute)
    
  \- *trial_rate_calc.rmd* 
  
  \- *rate_calc.csv*

| VARIABLE | DEFINITION | UNIT/OPTIONS |
| :--- | :--- | :--- | 
| Trial | trial type identifier (CTmax/CTmin) | max/min |
| Species | Poeciliid species of interest | gambusia/sailfin/swordtail |
| Size | size category (for male swordtails) | large/small/NA |
| Date | trial date | MM/DD/YYY |
| Rate | temperature rate calculated from original excel file | °C/min |
| file.id | name of original Excel file | NA |

   * Note: size is pertinent only to swordtails, in which males were categorized as small or large; all other individuals (including female swordtails) marked as NA

   * Note: temperature rate calculation can be found on the first sheet of the original excel file
   

**3. Thermal Tolerance Analysis:** overall analysis testing for sex differences in thermal tolerance between three Poecilidd fish species
    
  \- *Final_temperature-analysis.Rmd*
  
  \- *manuscript_figures-tables.R* (code for designing the figures/tables used in the manuscript)
  
  \- *Temp ID - Completed Data.csv*
  
| VARIABLE | DEFINITION | UNIT/OPTIONS |
| :--- | :--- | :--- | 
| Species | Poeciliid species of interest | X. nigrensis/G. affinis/P. latipinna |
| Mating | Mating style (if present) | courting/coercive/NA |
| Size.cat | Size category (for male swordtails) | large/small/reg |
| Sex | Sex of individual | male/female |
| Name | Name chosen for tracking individuals | NA |
| CTmax.num | Number given when fish is removed in CTmax trials | 1-9 |
| CTmin.num | Number given when fish is removed in CTmin trials | 1-9 |
| CTmax.C | Temperature of tank at loss of equilibrium in CTmax trials | °C |
| CTmin.C | Temperature of tank at loss of equilibrium in CTmin trials | °C |
| T.tot | CTmax.C temperature - CTmin.C temperature (total thermal range) | °C |
| Size.mm | Standard length of individual | mm |


**4. Outlier Removed Analysis:** a re-run of the overall analysis testing for sex differences in thermal tolerance between three Poecilidd fish species without major outliers
    
  \- *Rerun-outlier-removal.Rmd*
  
  \- Uses the *Temp ID - Completed Data.csv* located in the **Thermal Tolerance Analysis** folder


***
    
An additional folder -- **Raw Data** -- contains the raw excel files corresponding to each trial with additional trial notes, and the raw climatology data for the Austin Daily temperature analysis.

  \- Source for climatology data: Kantor D, Casey NW, Menne MJ, Buddenberg A (2023) Local Climatological Data (LCD), Version 2 [Subset used: Austin Camp Mabry, TX US (USW00013958), hourly dry bulb temperature, 2020–2024]. NOAA National Centers for Environmental Information. https://doi.org/10.25921/jp3d-3v19. Accessed 8 April 2025

   * See details at the above DOI for structure of raw climatology data
  
  \- Structure of raw trial excel files:
   
   * Sheet 1: created post-trial to calculate the temperature rate (°C/min) experienced in that trial. Time column (number of seconds since the start of temperature tracking) and Temp column (°C of tank at associated time) were copied from the 'Data In' sheet. The slope formula (=SLOPE(B/A * 60)) was calculated in a free cell.
   
   * Data In: sheet displaying the Arduino temperature recordings. 'Time' is the time of the day (HH:MM:SS.MS), 'CH1' is the rate calculated every 40s, 'CH2' is the number of seconds since the start of temperature tracking, 'CH3' is the °C of tank at associated time. When the button on the Ardunio system was pressed to record an individual, 'CH1' rate is replaced with organism number, 'CH2' is replaced with the temperature, and 'CH3' is replaced with the number of seconds since recording start. 
   
   * Data Out/Settings/Manifest: These sheets were automatically structured with the installation of the Arduino temperature system, and no alterations were made by us during our trials.
   

***
***

### **VERSION & PACKAGE INFORMATION**

All Rmd files list required packages at the beginning of the document. Below are all packages required to run all Rmd files.

| PACKAGE/SOFTWARE | VERSION | USED IN |
| :--- | :--- | :--- | 
| Rstudio running R | v2026.01.0+392/ v4.4.2 | All Rmd files | 
| `curl` | v5.2.1 | All Rmd files | 
| `dplyr` | v1.1.4 | All Rmd files | 
| `tidyr` | v1.3.2 | All Rmd files | 
| `ggplot2` | v3.5.1 | All Rmd files | 
| `gt` | v0.11.1 | All Rmd files |
| `lubridate` | v1.9.3 | *Austin-Daily-Temps.Rmd* |
| `summarytools` | v1.1.4 | *trial_rate_calc.Rmd*/ *Final_temperature-analysis.Rmd* |
| `MuMIn` | v1.48.4 | *Final_temperature-analysis.Rmd* | 
| `lmtest` | v3.1.3 | *Final_temperature-analysis.Rmd* | 
| `sandwich` | v3.1-1 | *Final_temperature-analysis.Rmd* | 
| `ggpubr` | v0.6.2 | *Final_temperature-analysis.Rmd* |
| `DHARMa` | v0.4.7 | *Final_temperature-analysis.Rmd* | 
| `emmeans` | v1.10.5 | *Final_temperature-analysis.Rmd* | 
| `car` | v3.1-2 | *Final_temperature-analysis.Rmd* | 
| `patchwork` | v1.3.2 | *Final_temperature-analysis.Rmd* | 


***
***











 
 
 